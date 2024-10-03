import os
import subprocess
from enum import Enum
from math import cos, pi, sin

from astropy.stats import knuth_bin_width

import h5py

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import mcpl

import numpy as np

import pandas as pd

from uncertainties import ufloat

c2 = 8.98755e16  # m2 / s2
hc = 1.23984186e-12  # MeV m


MCPLColumns = [
    "id",
    "type",
    "E",
    "x",
    "y",
    "z",
    "u",
    "v",
    "w",
    "t",
    "wgt",
    "px",
    "py",
    "pz",
    "userflags",
]
XUNITS = {
    "x": "cm",
    "y": "cm",
    "z": "cm",
    "R": "cm",
    "theta": "$^o$",
    "u": "",
    "v": "",
    "w": "",
    "phi": "$^o$",
    "psi": "$^o$",
    "mu": "",
    "E": "MeV",
    "ln(E0/E)": "",
    "lambda": r"$\AA$",
    "t": "ms",
    "log(t)": "ms",
}
YUNITS = {
    "x": "cm$^{-1}$",
    "y": "cm$^{-1}$",
    "z": "cm$^{-1}$",
    "R": "cm$^{-2}$",
    "theta": "rad$^{-1}$",
    "u": "",
    "v": "",
    "w": "",
    "phi": "rad$^{-1}$",
    "psi": "rad$^{-1}$",
    "mu": "",
    "E": "eV$^{-1}$",
    "ln(E0/E)": "",
    "lambda": r"$\AA^{-1}$",
    "t": "s$^{-1}$",
    "log(t)": "",
}
XLATEX = {
    "x": "x",
    "y": "y",
    "z": "z",
    "R": "R",
    "theta": r"\theta",
    "u": "u",
    "v": "v",
    "w": "w",
    "phi": r"\varphi",
    "psi": r"\vartheta",
    "mu": r"\mu=\cos(\vartheta)",
    "E": "E",
    "lambda": r"\lambda",
    "ln(E0/E)": r"\ln(E_0/E)",
    "t": "t",
    "log(t)": r"\log(t)",
}
YSCALE = {
    "x": 1,
    "y": 1,
    "z": 1,
    "R": 1,
    "theta": 1,
    "u": 1,
    "v": 1,
    "w": 1,
    "phi": 1,
    "psi": 1,
    "mu": 1,
    "E": 1e6,
    "ln(E0/E)": 1,
    "lambda": 1,
    "t": 1e-3,
    "log(t)": 1,
}
PARTMASS = {2112: 939.5653, 22: 0, 11: 0.511, 2212: 938.7829}


def momentum(ptype, pekin):
    mc2 = np.array([PARTMASS[np.abs(i)] for i in ptype])
    pc = pekin * np.sqrt(1.0 + 2 * mc2 / pekin)
    return pc


def velocity(ptype, pekin):
    mc2 = np.array([PARTMASS[np.abs(i)] for i in ptype])
    pc = momentum(ptype, pekin)
    v = np.array(
        [pci / mc2i if mc2i > 0 else c2 ** 0.5 for pci, mc2i in zip(pc, mc2)]
    )
    return v


def wavelength(ptype, pekin):
    return hc / momentum(ptype, pekin) * 1e10


def resample_quadratic(width, center, size):
    u1 = np.random.uniform(-1, 1, size)
    u2 = np.random.uniform(-1, 1, size)
    u3 = np.random.uniform(-1, 1, size)
    mask = (np.abs(u3) >= np.abs(u2)) & (np.abs(u3) >= np.abs(u1))
    newu2 = u2[mask]
    newu3 = u3[~mask]
    samples = np.concatenate([newu2, newu3])
    return center + samples * width / 2


class PDGCode(Enum):
    NEUTRON = 2112
    PHOTON = 22
    ELECTRON = 11
    POSITRON = -11
    PROTON = 2212


class OpenMCCode(Enum):
    NEUTRON = 0
    PHOTON = 1
    ELECTRON = 2
    POSITRON = 3
    PROTON = 4


class VitessCode(Enum):
    NEUTRON = "N"
    PHOTON = "G"
    ELECTRON = "E-"
    POSITRON = "E+"
    PROTON = "P"


class SurfaceSourceFile:
    """Read a surface source file.
    Possible formats: MCPL, HDF5 (OpenMC), SSV (KDSource)

    Parameters
    ----------
    filepath: str
        Path to surface source file.
    S0: float
        Simulation source factor.
        Default: 1.0
    Nparticles: int
        Total source particles simulated.
        Default: 1e9
    dA: float
        Area of the surface [cm$^2$].
        Default: None
    translation: list
        Translation for the position variables [cm].
        Default: [0, 0, 0]
    rotation: list
        Rotation for the position and direction variables [deg].
        Default: [0, 0, 0]
    domain: dict
        Range of the variables: var_min <= 'var' < var_max.
        It must be defined like: {'var':[var_min, var_max]}.
        List of possible variables:
        type [PDG],
        E [MeV], ln(E0/E), lambda [AA],
        x [cm], y [cm], z [cm], R [cm], theta [rad],
        u, v, w, phi [rad], psi [rad],
        t [ms], log(t),
        wgt.
        Default: {}
    set_domain_first: bool
        Define if the variables domains and the normal direction
        must be set before or after the translation and rotation.
        Default: False
    set_rotation_first: bool
        Define if the rotation must be done before or after
        the translation of the coordinates.
        Default: False
    E0: float
        Reference energy [MeV] to calculate the lethargy ln(E0/E).
        Default: 20
    convoluted: bool
        Define if the simulation was performed with a convoluted pulse or not.
        Default: False
    uvw_reference: str
        Versor which is normal to the surface source after the rotation.
        Default: 'w'
    surface: int
        Define the surface id to extract the particles from the tracks file.
        Default: None
    Nmax: int
        Define the maximum number of particles to be readed.
    current: float
        Current, in mA, to normalize the results.
    tpulse: float
        Pulse width, in s, to normalize the results.
    """

    def __init__(
        self,
        filepath,
        S0=1.0,
        Nparticles=1e9,
        dA=None,
        translation=[0, 0, 0],
        rotation=[0, 0, 0],
        domain={},
        set_domain_first=False,
        set_rotation_first=False,
        E0=20,
        convoluted=False,
        uvw_reference="w",
        surface=None,
        Nmax=1e9,
        current=None,
        tpulse=None,
        skip_cloned=[],
        pulse_shape="rectangular",
    ):
        # Initialize SurfaceSource class atributes
        self._filepath = filepath
        self._extension = os.path.splitext(self._filepath)[-1]
        self._translation = translation
        self._rotation = rotation
        self._domain = domain
        self._domain_first = set_domain_first
        self._rotation_first = set_rotation_first
        self._E0 = E0
        self._uvw = uvw_reference
        self._S0 = S0
        self._Np = int(Nparticles)
        self._dA = dA
        self._surface = surface
        self._max_particles = int(Nmax)
        self._tpulse = tpulse
        self._current = current
        self._convoluted = convoluted
        self._shape = pulse_shape
        self._cloned = skip_cloned
        self._df = self.__read__()
        self._df2 = self.get_pandas_dataframe()
        self._brilliance = 1.0

    def __read__(self):
        # OpenMC .h5 format
        if self._extension == ".h5":
            with h5py.File(self._filepath, "r") as fh:
                df = pd.DataFrame(columns=MCPLColumns)
                # Change from OpenMC ParticleType type to MCPL PDGCode
                df["type"] = fh["source_bank"]["particle"]
                for pin, pout in zip((OpenMCCode), (PDGCode)):
                    df.loc[df["type"] == pin.value, "type"] = pout.value
                df["id"] = df.index
                df["E"] = fh["source_bank"]["E"] * 1e-6  # eV to MeV
                df["x"] = fh["source_bank"]["r"]["x"]  # cm
                df["y"] = fh["source_bank"]["r"]["y"]  # cm
                df["z"] = fh["source_bank"]["r"]["z"]  # cm
                df["u"] = fh["source_bank"]["u"]["x"]
                df["v"] = fh["source_bank"]["u"]["y"]
                df["w"] = fh["source_bank"]["u"]["z"]
                df["wgt"] = fh["source_bank"]["wgt"]
                df["t"] = fh["source_bank"]["time"] * 1e3  # s to ms
                df["px"] = 0.0
                df["py"] = 0.0
                df["pz"] = 0.0
                df["userflags"] = 0x00000000
                df = df[df["id"] < self._max_particles]
        # MCPL format
        elif self._extension == ".gz" or self._extension == ".mcpl":
            with mcpl.MCPLFile(
                self._filepath, blocklength=self._max_particles
            ) as mcpl_file:
                pb = mcpl_file.read_block()
                plist = np.array(
                    [
                        pb.pdgcode,
                        pb.ekin,
                        pb.x,
                        pb.y,
                        pb.z,
                        pb.ux,
                        pb.uy,
                        pb.uz,
                        pb.time,
                        pb.weight,
                        pb.polx,
                        pb.poly,
                        pb.polz,
                        pb.userflags,
                    ]
                )
                df = pd.DataFrame(plist.T, columns=MCPLColumns[1:])
                if self._surface is not None:
                    df = df[df["userflags"] == int(self._surface)]
                df["id"] = np.arange(0, len(df))
                df = df[MCPLColumns]
        # VITESS format (for detectors, not for dump files)
        elif self._extension == ".vitess":
            df = pd.DataFrame(
                np.loadtxt(self._filepath, usecols=list(np.arange(3, 15))),
                columns=[
                    "t",
                    "E",
                    "wgt",
                    "z",
                    "x",
                    "y",
                    "w",
                    "u",
                    "v",
                    "px",
                    "py",
                    "pz",
                ],
            )
            df["E"] = (72.3 / (252.8 * df["E"].to_numpy())) ** 2 * 1e-6
            df["id"] = [
                int(i[2:]) - 1
                for i in np.loadtxt(self._filepath, dtype="str", usecols=(0))
            ]
            df["type"] = np.loadtxt(self._filepath, dtype="str", usecols=(1))
            df["userflags"] = 0x00000000
            for pin, pout in zip((VitessCode), (PDGCode)):
                df.loc[df["type"] == pin.value, "type"] = pout.value
            df = df[MCPLColumns]
        # ASCII format, for all the other cases
        else:
            df = pd.DataFrame(
                np.loadtxt(self._filepath, skiprows=5), columns=MCPLColumns
            )

        if len(self._cloned) > 0:
            df = self.search_duplicated_particles(df, self._cloned)

        # Print some info regarding the particles in the file
        print(
            "Number of particles in file {:s}: {:d}".format(
                self._filepath, len(df)
            )
        )
        for p, idx in zip((PDGCode), ["n", "g", "e-", "e+", "p"]):
            print(
                "{:.0f}% {:s}, ".format(
                    len(df[df["type"] == p.value]) / len(df["type"]) * 100, idx
                ),
                end="",
            )
        print("")
        return df

    def get_pandas_dataframe(self):
        """Obtain a DataFrame with the particles that were saved in the file"""
        df = self._df.copy()
        # Define the coordinates system
        if self._uvw == "u":
            X, Y, Z, U, V, W = "y", "z", "x", "v", "w", "u"
        elif self._uvw == "v":
            X, Y, Z, U, V, W = "z", "x", "y", "w", "u", "v"
        elif self._uvw == "w":
            X, Y, Z, U, V, W = "x", "y", "z", "u", "v", "w"
        self._X = X
        self._Y = Y
        self._Z = Z
        # Check domain
        if self._domain_first:
            df = self._set_domain(df)
        # Translate and rotate the coordinates
        if self._rotation_first:
            df = self._set_rotation(df)
            df = self._set_translation(df)
        else:
            df = self._set_translation(df)
            df = self._set_rotation(df)
        # Normalize the direction
        self._normalize_direction(df)
        # Convolute time considering the pulse width
        if self._tpulse is not None and self._convoluted is False:
            self._convolute(
                df, {"t": self._tpulse * 1e3}, self._shape
            )  # from s to ms
            self._convoluted = True
        # Calculate extra coordinates
        df["R"] = (df[X].to_numpy() ** 2 + df[Y].to_numpy() ** 2) ** 0.5
        df["theta"] = np.arctan2(df[Y].to_numpy(), df[X].to_numpy())
        df["mu"] = df[W].to_numpy()
        df["psi"] = np.arccos(df[W].to_numpy())
        df["phi"] = np.arctan2(df[V].to_numpy(), df[U].to_numpy())
        df["ln(E0/E)"] = np.log(self._E0 / df["E"].to_numpy())
        df["log(t)"] = np.log10(df["t"].to_numpy())
        # df['lambda'] = 72.3 / 252.8 * (df['E'] * 1e6)**-0.5
        df["lambda"] = wavelength(df["type"].to_numpy(), df["E"].to_numpy())
        # Check domain
        if not self._domain_first:
            df = self._set_domain(df)
        return df

    def search_duplicated_particles(self, df, vars):
        newwgt = (
            df.groupby(
                vars,
                sort=False,
                as_index=False,
            )
            .sum(numeric_only=True)["wgt"]
            .values
        )
        newdf = df.drop_duplicates(vars, ignore_index=True).copy()
        newdf["wgt"] = newwgt
        return newdf

    def _set_translation(self, df):
        # Translate the coordinates (only position)
        df["x"] = df["x"].to_numpy() + self._translation[0]
        df["y"] = df["y"].to_numpy() + self._translation[1]
        df["z"] = df["z"].to_numpy() + self._translation[2]
        return df

    def _set_rotation(self, df):
        # Rotate the coordinates (posisition and direction)
        phi, theta, psi = np.array(self._rotation) * (pi / 180)
        c3, s3 = cos(phi), sin(phi)
        c2, s2 = cos(theta), sin(theta)
        c1, s1 = cos(psi), sin(psi)
        rotation_matrix = np.array(
            [
                [c1 * c2, c1 * s2 * s3 - c3 * s1, s1 * s3 + c1 * c3 * s2],
                [c2 * s1, c1 * c3 + s1 * s2 * s3, c3 * s1 * s2 - c1 * s3],
                [-s2, c2 * s3, c2 * c3],
            ]
        )
        df["x"], df["y"], df["z"] = rotation_matrix @ np.array(
            [df["x"].to_numpy(), df["y"].to_numpy(), df["z"].to_numpy()]
        )
        df["u"], df["v"], df["w"] = rotation_matrix @ np.array(
            [df["u"].to_numpy(), df["v"].to_numpy(), df["w"].to_numpy()]
        )
        return df

    def _normalize_direction(self, df):
        # Normalize the direction versor
        unorm = (
            df["u"].to_numpy() ** 2
            + df["v"].to_numpy() ** 2
            + df["w"].to_numpy() ** 2
        )

        df["u"], df["v"], df["w"] = (
            df["u"].to_numpy() / unorm ** 0.5,
            df["v"].to_numpy() / unorm ** 0.5,
            df["w"].to_numpy() / unorm ** 0.5,
        )
        return df

    def _set_domain(self, df):
        # Check the filters coordinates
        for pvar, (pmin, pmax) in self._domain.items():
            if pvar == "psi" or pvar == "phi" or pvar == "theta":
                pmin, pmax = np.deg2rad([pmin, pmax])
            if pmin is not None:
                df = df[df[pvar] >= pmin]
            if pmax is not None:
                df = df[df[pvar] < pmax]
        return df

    def _convolute(self, df, convolution, shape):
        # Convolute the desired variables
        # Set the seed to reproduce the results
        np.random.seed(len(df))
        if shape == "rectangular":
            for cvar, cwidth in convolution.items():
                df[cvar] = df[cvar].to_numpy() + np.random.uniform(
                    low=0.0, high=cwidth, size=len(df[cvar])
                )
        elif shape == "triangular":
            for cvar, cwidth in convolution.items():
                df[cvar] = df[cvar].to_numpy() + np.random.triangular(
                    left=0.0, mode=cwidth / 2, right=cwidth, size=len(df[cvar])
                )
        elif shape == "quadratic":
            for cvar, cwidth in convolution.items():
                df[cvar] = df[cvar].to_numpy() + resample_quadratic(
                    width=cwidth, center=cwidth / 2, size=len(df[cvar])
                )

        return df

    def _get_brilliance(self, df, norm_vars):
        self._brilliance = 1.0
        Sunit = self._get_units(norm_vars)

        if self._dA is not None:
            self._brilliance /= self._dA

        if "mAs" in norm_vars:
            self._brilliance /= 1.6022e-19 * 1e3
            if self._current is not None:
                self._brilliance *= self._current
            if self._tpulse is not None:
                self._brilliance *= self._tpulse
            norm_vars = [x for x in norm_vars if x not in ["mAs"]]

        for nvar in norm_vars:
            pmin = min(df[nvar]) * YSCALE[nvar]
            pmax = max(df[nvar]) * YSCALE[nvar]
            if nvar == "R":
                dp = 0.5 * (pmax ** 2 - pmin ** 2)
            elif nvar == "psi":
                dp = -(np.cos(pmax) - np.cos(pmin))
            else:
                dp = pmax - pmin

            self._brilliance /= dp

        return Sunit

    def _get_units(self, norm_vars):
        Sunit = ""
        dA = False
        dOmega = False

        if "mAs" in norm_vars:
            if self._current is None and self._tpulse is None:
                Sunit += "(mA s)$^{-1}$ "
            elif self._current is None:
                Sunit += "mA$^{-1}$ "
            elif self._tpulse is None:
                Sunit += "s$^{-1}$ "
            else:
                if "t" not in norm_vars:
                    Sunit += "pulse$^{-1}$ "
            norm_vars = [x for x in norm_vars if x not in ["mAs"]]

        if self._dA is not None:
            Sunit += "cm$^{-2}$ "
            dA = True

        if "x" in norm_vars and "y" in norm_vars and dA is False:
            Sunit += "cm$^{-2}$ "
            norm_vars = [x for x in norm_vars if x not in ["x", "y"]]

        if "x" in norm_vars and "z" in norm_vars and dA is False:
            Sunit += "cm$^{-2}$ "
            norm_vars = [x for x in norm_vars if x not in ["x", "z"]]

        if "y" in norm_vars and "z" in norm_vars and dA is False:
            Sunit += "cm$^{-2}$ "
            norm_vars = [x for x in norm_vars if x not in ["y", "z"]]

        if "R" in norm_vars and "theta" in norm_vars and dA is False:
            Sunit += "cm$^{-2}$ "
            norm_vars = [x for x in norm_vars if x not in ["R", "theta"]]

        if "mu" in norm_vars and "phi" in norm_vars and dOmega is False:
            Sunit += "sr$^{-1}$ "
            norm_vars = [x for x in norm_vars if x not in ["mu", "phi"]]
            dOmega = True

        if "psi" in norm_vars and "phi" in norm_vars and dOmega is False:
            Sunit += "sr$^{-1}$ "
            norm_vars = [x for x in norm_vars if x not in ["psi", "phi"]]
            dOmega = True

        for nvar in norm_vars:
            Sunit += "{:s} ".format(YUNITS[nvar])

        return Sunit

    def _get_info(self, df, vars):
        info = "Surface: {:d}\n\n".format(
            self._surface if self._surface is not None else 0
        )
        info += "NSIM: {:.2e}\n\n".format(self._Np)
        info += "NMCPL: {:d}\n\n".format(len(df))
        info += "S0: {:.2e}\n\n".format(self._S0)
        info += "Area: {:.2e} cm$^2$\n\n".format(
            self._dA if self._dA is not None else 0.0
        )
        info += "Translation: [{:.1f}, {:.1f}, {:.1f}]\n\n".format(
            self._translation[0], self._translation[1], self._translation[2]
        )
        info += "Rotation: [{:.1f}, {:.1f}, {:.1f}]\n\n".format(
            self._rotation[0], self._rotation[1], self._rotation[2]
        )
        for pvar in vars:
            if pvar != "mAs":
                pmin = min(df[pvar])
                pmax = max(df[pvar])
                if pvar == "theta" or pvar == "psi" or pvar == "phi":
                    pmin, pmax = np.rad2deg([pmin, pmax])
                info += (
                    "${:s}$: [{:.2e}, {:.2e}] {:s}".format(
                        XLATEX[pvar], pmin, pmax, XUNITS[pvar]
                    )
                    + "\n\n"
                )
            else:
                if self._current is not None:
                    info += "Current: {:.2f} mA\n\n".format(self._current)
                if self._tpulse is not None:
                    info += "Pulse width: {:.2f} us\n\n".format(
                        self._tpulse * 1e6
                    )
        return info

    def get_distribution(
        self,
        vars,
        bins,
        scales=["linear"],
        factor=1.0,
        filters={},
        norm_vars=[],
        convolution=None,
        total=False,
    ):
        """Obtain a N-D distribution for given variables

        Parameters
        ----------
        vars: list of str
            Variables to calculate the distribution
        bins: list of int or list of iterable
            Number or array of bins to calculate the distribution.
            If int, then the distributions are going to be calculated between
            min(var), max(var) for int bins.
            If 0, the number of bins will be calculated optimally with
            the knuth_bin_width function from astropy.
        scales: list of str
            Scales of the axis. The options are linear or log.
            Default: ['linear']
        factor: float
            Multiply the distribution by a specific numerical factor.
            Default: 1.0
        filters: dict
            Specify the other coordinates domain for the distribution
            It must be defined like: {'var':[var_min, var_max]}
            Default: {}
        norm_vars: list of str
            List of the variables to calculate the brilliance, i.e. to
            normalize the distribution.
            List of possible variables:
            mAs [mAs], E [MeV], ln(E0/E), lambda [AA],
            x [cm], y [cm], z [cm], R [cm], theta [rad],
            u, v, w, phi [rad], psi [rad],
            t [ms], log(t),
        convolution: dict
            Specify the coordinates to perform a convolution.
            It must be defined like: {'var':[var_convolution]}
            Default: {}
        total: bool
            Define if the integral value is required.
            Default: False

        Returns
        -------
        p_df: DataFrame
            Pandas DataFrame with the bins and the distribution.
        bins: list of arrays
            List of the bins where the distributions were calculated.
        p_info: str
            String with the information of the distributions calculated.
        p_total: ufloat
            If total is True, the integral value of the distribution for the
            required normalization variables.
        p_units: str
            If total is True, the units for the integral value of the
            distribution.
        """
        p_info = ""
        old_filters = self._domain
        self._domain = filters
        # df = self.get_pandas_dataframe()
        df = self._df2.copy()
        df = self._set_domain(df)
        self._domain = old_filters

        # Normalize the distribution by the brilliance and convolute
        if convolution is not None:
            df = self._convolute(df, convolution)
        p_units = self._get_brilliance(
            df, [x for x in norm_vars if x not in vars]
        )

        if len(bins) != len(vars):
            print("Both lists of variables and bins should have the same size")
            return

        xmed = [0] * len(bins)
        xleft = [0] * len(bins)
        xright = [0] * len(bins)
        dx = [0] * len(bins)

        for i, (bin, var, scale) in enumerate(zip(bins, vars, scales)):
            # If bins is int, create a mesh from var-min to var-max
            if type(bin) is int:
                if bin == 0:
                    if scale == "log":
                        bins[i] = knuth_bin_width(
                            np.log10(df[var].to_numpy()), return_bins=True
                        )[1]
                        bins[i] = 10 ** bins[i]
                    else:
                        bins[i] = knuth_bin_width(
                            df[var].to_numpy(), return_bins=True
                        )[1]

                else:
                    if scale == "log":
                        bins[i] = np.logspace(
                            np.log10(df[var].min()),
                            np.log10(df[var].max()),
                            bin,
                        )
                    else:
                        bins[i] = np.linspace(
                            df[var].min(), df[var].max(), bin
                        )
            else:
                # If var is angle, convert degrees to radians
                if var == "psi" or var == "phi" or var == "theta":
                    bins[i] = np.deg2rad(bins[i])

        # Calculate the histograms and the total current
        for i, bin in enumerate(bins):
            xleft[i] = bin[:-1]
            xright[i] = bin[1:]
            xmed[i] = 0.5 * (xleft[i] + xright[i])
            dx[i] = xright[i] - xleft[i]

        xleft = np.meshgrid(*xleft, indexing="ij")
        xright = np.meshgrid(*xright, indexing="ij")
        xmed = np.meshgrid(*xmed, indexing="ij")
        dx = np.meshgrid(*dx, indexing="ij")

        p_mean = np.histogramdd(
            sample=df[vars].to_numpy(), bins=bins, weights=df["wgt"]
        )[0]
        p_stdv = (
            np.histogramdd(
                sample=df[vars].to_numpy(), bins=bins, weights=df["wgt"] ** 2
            )[0]
            ** 0.5
        )
        p_total = ufloat(sum(df["wgt"]), sum(df["wgt"] ** 2) ** 0.5)

        # p_mean, p_stdv = p_mean.T, p_stdv.T

        # Normalize the distribution considering the source factor,
        # the number of source particles, the brilliance, and the extra factor
        p_mean *= factor * self._brilliance * self._S0 / self._Np
        p_stdv *= factor * self._brilliance * self._S0 / self._Np
        p_total *= factor * self._brilliance * self._S0 / self._Np

        # Calculate the relative error
        with np.errstate(divide="ignore", invalid="ignore"):
            p_erel = np.nan_to_num(p_stdv / p_mean)

        for i, (bin, var) in enumerate(zip(bins, vars)):
            if var in norm_vars:
                p_mean = p_mean / YSCALE[var]
                p_stdv = p_stdv / YSCALE[var]
                if var == "psi":
                    p_mean = p_mean / (-np.cos(xright[i]) + np.cos(xleft[i]))
                    p_stdv = p_stdv / (-np.cos(xright[i]) + np.cos(xleft[i]))
                elif var == "R":
                    p_mean = p_mean / (xright[i] ** 2 / 2 - xleft[i] ** 2 / 2)
                    p_stdv = p_stdv / (xright[i] ** 2 / 2 - xleft[i] ** 2 / 2)
                else:
                    p_mean = p_mean / dx[i]
                    p_stdv = p_stdv / dx[i]

            # Convert the angles from radians to degrees
            if var == "psi" or var == "phi" or var == "theta":
                bins[i] = np.rad2deg(bins[i])
                xleft[i] = np.rad2deg(xleft[i])
                xright[i] = np.rad2deg(xright[i])
                xmed[i] = np.rad2deg(xmed[i])

        p_info += self._get_info(df, norm_vars)
        if len(bins) == 1:
            p_info += "Max: {:.3e} at {:.3e} {:s}\n\n".format(
                max(p_mean), xmed[0][p_mean == max(p_mean)][0], XUNITS[vars[0]]
            )
        if type(factor) is int or type(factor) is float:
            p_info += "Integral: {:.3e} [{:s}]\n".format(p_total, p_units)

        # Save the results into a Dataframe
        p_df = pd.DataFrame()
        for i, (var, xl, xr, xm) in enumerate(zip(vars, xleft, xright, xmed)):
            p_df["{:s}-min".format(var)] = xl.ravel()
            p_df["{:s}-max".format(var)] = xr.ravel()
            p_df["{:s}".format(var)] = xm.ravel()

        p_df["mean"] = p_mean.ravel()
        p_df["stdv"] = p_stdv.ravel()
        p_df["erel"] = p_erel.ravel()

        if total:
            return p_total, p_units
        else:
            return p_df, bins, p_info

    def get_domain(
        self,
        vars,
        filters={},
    ):
        old_filters = self._domain
        self._domain = filters
        df = self._df2.copy()
        df = self._set_domain(df)
        self._domain = old_filters

        dm = {}
        for var in vars:
            dm[var] = df[var].min(), df[var].max()
            if var == "theta" or var == "psi" or var == "phi":
                dm[var] = np.rad2deg([dm[var][0], dm[var][1]])
        return dm

    def plot_distribution(
        self,
        vars,
        bins,
        scales=["linear", "linear"],
        factor=1.0,
        filters={},
        norm_vars=[],
        errors=True,
        convolution=None,
        zlevels=1,
        info=False,
        ylabel=None,
        zlabel=None,
        zscale="linear",
        tolerance=1.0,
        vmin=None,
        vmax=None,
        peak_brilliance=False,
        **kwargs,
    ):
        """Plot 1-D or 2-D distribution for given variables

        Parameters
        ----------
        vars: list of str
            Variables to calculate the distribution
        bins: list of int or list of iterable
            Number or array of bins to calculate the distribution.
            If int, then the distributions are going to be calculated between
            min(var), max(var) for int bins.
            If 0, the number of bins will be calculated optimally with
            the knuth_bin_width function from astropy.
        scales: list of str
            Scales of the axis. The options are linear or log.
            Default: ['linear', 'linear']
        factor: float
            Multiply the distribution by a specific numerical factor.
            Default: 1.0
        filters: dict
            Specify the other coordinates domain for the distribution
            It must be defined like: {'var':[var_min, var_max]}
            Default: {}
        norm_vars: list of str
            List of the variables to calculate the brilliance, i.e. to
            normalize the distribution.
        errors: bool
            Specify if the error bars should be plotted or not.
            Default: True
        convolution: dict
            Specify the coordinates to perform a convolution.
            It must be defined like: {'var':[var_convolution]}
            Default: {}
        zlevels: int
            Number of countour surfaces to plot with the colormap (if more
            than one variable is required). Deprecated when zscale is 'log'
            Default: 1
        info: bool
            Specify if the information of the plot should be included or not.
            Default: False
        ylabel: str
            Label for the y-axis.
            Default: 'Intensity' (for 1-D plots)
        zlabel: str
            Label for the z-axis (2-D plots).
            Default: 'Intensity'
        zscale: str
            Scale for the colorbar (2-D plots).
            Default: 'linear'
        tolerance: float
            When the realive error is greater or equal than the tolerance,
            the mean values of the distribution are not shown.
            Default: 1.0
        **kwargs
            Extra arguments that will be passed to matplotlib
            (label, color, and so on)

        Returns
        -------
        matplotlib 1-D or 2-D plot for the required variables with the proper
        normalization.
        """
        xscale, yscale = scales
        scales = [scales[i] for i in range(0, len(bins))]
        df, bins, pinfo = self.get_distribution(
            vars, bins, scales, factor, filters, norm_vars, convolution
        )
        if peak_brilliance and len(vars) == 2:
            df = (
                df.groupby([var for var in vars if var != "t"])
                .max()
                .reset_index()
            )
        if xscale == "log":
            norm_vars = [var for var in norm_vars if var != vars[0]]
        if yscale == "log" and len(vars) > 1:
            norm_vars = [var for var in norm_vars if var != vars[1]]
        Jlabel = self._get_units(norm_vars)

        if len(bins) == 1:
            Jbrill = "Intensity" if ylabel is None else "{:s}".format(ylabel)
            Jlabel = (
                "{:s}".format(Jbrill)
                if Jlabel == ""
                else "{:s} [ {:s} ]".format(Jbrill, Jlabel)
            )
            if xscale == "log":
                if info:
                    print("Plotting x*f(x) instead of f(x) (xscale='log')")
                df["mean"] = (
                    df["mean"] * df["{:s}".format(vars[0])] * YSCALE[vars[0]]
                )
                df["stdv"] = df["erel"] * df["mean"]
            df.loc[df["erel"] >= tolerance, "mean"] = 0
            df.loc[df["erel"] >= tolerance, "stdv"] = 0
            if errors:
                plt.errorbar(
                    df[vars], df["mean"], df["stdv"], ds="steps-mid", **kwargs
                )
            else:
                plt.plot(
                    df[vars].to_numpy(),
                    df["mean"].to_numpy(),
                    ds="steps-mid",
                    **kwargs,
                )
            plt.xscale(xscale)
            plt.yscale(yscale)
            if XUNITS[vars[0]] != "":
                plt.xlabel(
                    r"${:s}$ [{:s}]".format(XLATEX[vars[0]], XUNITS[vars[0]])
                )
            else:
                plt.xlabel(r"${:s}$".format(XLATEX[vars[0]]))
            plt.ylabel(r"{:s}".format(Jlabel))

        elif len(bins) == 2:
            Jbrill = "Intensity" if zlabel is None else zlabel
            Jlabel = (
                "{:s}".format(Jbrill)
                if Jlabel == ""
                else "{:s} [ {:s} ]".format(Jbrill, Jlabel)
            )

            len1 = len(bins[0]) - 1
            len2 = len(bins[1]) - 1
            bins1 = np.array(bins[0])
            bins2 = np.array(bins[1])

            if xscale == "log" and yscale == "log":
                if info:
                    print(
                        "Plotting x*y*f(x,y) instead of f(x,y)\
                        (xscale='log', yscale='log')"
                    )
                df["mean"] = (
                    df["mean"] * df["{:s}".format(vars[0])] * YSCALE[vars[0]]
                )
                df["mean"] = (
                    df["mean"] * df["{:s}".format(vars[1])] * YSCALE[vars[1]]
                )
                df["stdv"] = df["erel"] * df["mean"]
            elif xscale == "log":
                if info:
                    print(
                        "Plotting x*f(x,y) instead of f(x,y)\
                        (xscale='log')"
                    )
                df["mean"] = (
                    df["mean"] * df["{:s}".format(vars[0])] * YSCALE[vars[0]]
                )
                df["stdv"] = df["erel"] * df["mean"]

            elif yscale == "log":
                if info:
                    print(
                        "Plotting y*f(x,y) instead of f(x,y)\
                        (yscale='log')"
                    )
                df["mean"] = (
                    df["mean"] * df["{:s}".format(vars[1])] * YSCALE[vars[1]]
                )
                df["stdv"] = df["erel"] * df["mean"]

            df.loc[df["erel"] >= tolerance, "mean"] = 0
            df.loc[df["erel"] >= tolerance, "stdv"] = 0
            df.loc[df["mean"] == 0, "mean"] = np.nan

            z_mean = df["mean"].to_numpy().reshape(len1, len2).T
            # z_stdv = df['stdv'].to_numpy().reshape(len2, len1)
            if zscale == "linear":
                plt.pcolor(
                    0.5 * (bins1[1:] + bins1[:-1]),
                    0.5 * (bins2[1:] + bins2[:-1]),
                    z_mean,
                    shading="auto",
                    vmin=vmin,
                    vmax=vmax,
                    **kwargs,
                )
                cbar = plt.colorbar(label=r"{:s}".format(Jlabel))
                if zlevels > 0:
                    lvls = np.linspace(cbar.vmin, cbar.vmax, zlevels)
            elif zscale == "log":
                plt.pcolor(
                    0.5 * (bins1[1:] + bins1[:-1]),
                    0.5 * (bins2[1:] + bins2[:-1]),
                    z_mean,
                    shading="auto",
                    norm=LogNorm(vmin, vmax),
                    **kwargs,
                )
                cbar = plt.colorbar(label=r"{:s}".format(Jlabel))
                if zlevels > 0:
                    lvls = np.logspace(
                        np.floor(np.log10(cbar.vmin)),
                        np.ceil(np.log10(cbar.vmax)),
                        zlevels
                        * int(
                            np.ceil(np.log10(cbar.vmax))
                            - np.floor(np.log10(cbar.vmin))
                        )
                        + 1,
                    )

            plt.xscale(xscale)
            plt.yscale(yscale)
            if XUNITS[vars[0]] != "":
                plt.xlabel(
                    r"${:s}$ [{:s}]".format(XLATEX[vars[0]], XUNITS[vars[0]])
                )
            else:
                plt.xlabel(r"${:s}$".format(XLATEX[vars[0]]))
            if XUNITS[vars[1]] != "":
                plt.ylabel(
                    r"${:s}$ [{:s}]".format(XLATEX[vars[1]], XUNITS[vars[1]])
                )
            else:
                plt.ylabel(r"${:s}$".format(XLATEX[vars[1]]))

            if zlevels > 0:
                cntrs = plt.contour(
                    0.5 * (bins[0][1:] + bins[0][:-1]),
                    0.5 * (bins[1][1:] + bins[1][:-1]),
                    z_mean,
                    colors="grey",
                    levels=lvls,
                    vmin=cbar.vmin,
                    vmax=cbar.vmax,
                )
                cbar.add_lines(cntrs)

        if info:
            plt.figtext(
                x=1.0, y=0.1, s=r"{:s}".format(pinfo), fontdict={"size": 8}
            )

    def save_source_file(self, filepath, **kwargs):
        """Save the particles list.
        Possible formats: MCPL, HDF5 (OpenMC), SSV (KDSource)

        Parameters
        ----------
        filepath: str
            Path to the surface source file.
        **kwargs:
            Extra arguments for the .h5 format

        """
        create_source_file(self._df2.copy(), filepath, **kwargs)


def create_source_file(df, filepath, **kwargs):
    """Generate a source file from a Pandas DataFrame.
    Possible formats: MCPL, HDF5 (OpenMC), SSV (KDSource)

    Parameters
    ----------
    df: DataFrame
        Pandas DataFrame that contains the particles list.
    filepath: str
        Path to surface source file.
    """

    new_extension = os.path.splitext(filepath)[-1]
    df = df[MCPLColumns]

    # Write in the OpenMC .h5 format
    if new_extension == ".h5":
        print("Saving into OpenMC format (HDF5)")
        for pin, pout in zip((PDGCode), (OpenMCCode)):
            df.loc[df["type"] == pin.value, "type"] = pout.value

        pos_dtype = np.dtype([("x", "<f8"), ("y", "<f8"), ("z", "<f8")])
        source_dtype = np.dtype(
            [
                ("r", pos_dtype),
                ("u", pos_dtype),
                ("E", "<f8"),
                ("time", "<f8"),
                ("wgt", "<f8"),
                ("delayed_group", "<i4"),
                ("surf_id", "<i4"),
                ("particle", "<i4"),
            ]
        )
        arr = np.array(
            [
                (
                    (s[3], s[4], s[5]),
                    (s[6], s[7], s[8]),
                    (s[2] * 1e6),
                    (s[9] * 1e-3),
                    (s[10]),
                    (0),
                    (0),
                    (s[1]),
                )
                for s in df.values
            ],
            dtype=source_dtype,
        )
        kwargs.setdefault("mode", "w")
        with h5py.File(filepath, **kwargs) as fh:
            fh.attrs["filetype"] = np.string_("source")
            fh.create_dataset("source_bank", data=arr, dtype=source_dtype)

    # Write in MCPL format
    elif new_extension == ".mcpl" or new_extension == ".gz":
        new_extension = ".mcpl"
        print("Saving into MCPL format")
        create_source_file(df, "temp.txt")
        subprocess.call(["ssv2mcpl", "temp.txt", filepath])
        subprocess.call(["rm", "temp.txt"])

    # Write the ASCII-based format file in the other cases
    else:
        print("Saving into SSV format (ASCII)")
        with open(filepath, "w") as fo:
            sheader = ""
            fmtstr = ""

            sheader += "#MCPL-ASCII\n"
            sheader += "#GENERATED FROM KDSOURCE\n"
            sheader += "#NPARTICLES: {:d}\n".format(len(df))
            sheader += "#END-HEADER\n"
            sheader += "index     "
            sheader += "pdgcode               "
            sheader += "ekin[MeV]                   "
            sheader += "x[cm]                   "
            sheader += "y[cm]                   "
            sheader += "z[cm]                      "
            sheader += "ux                      "
            sheader += "uy                      "
            sheader += "uz                "
            sheader += "time[ms]                  "
            sheader += "weight                   "
            sheader += "pol-x                   "
            sheader += "pol-y                   "
            sheader += "pol-z  "
            sheader += "userflags\n"
            fo.write(sheader)

            fmtstr += "%5i %11i %23.18g %23.18g %23.18g %23.18g "
            fmtstr += "%23.18g %23.18g %23.18g %23.18g %23.18g "
            fmtstr += "%23.18g %23.18g %23.18g 0x%08x\n"
            for s in df.values:
                fo.write(
                    fmtstr
                    % (
                        s[0],
                        s[1],
                        s[2],
                        s[3],
                        s[4],
                        s[5],
                        s[6],
                        s[7],
                        s[8],
                        s[9],
                        s[10],
                        s[11],
                        s[12],
                        s[13],
                        int(s[14]),
                    )
                )
    print("Done, saved into {:s} file".format(filepath))
