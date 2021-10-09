#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for PList object
"""

from xml.etree import ElementTree as ET
import os, subprocess

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st
import mcpl

from .utils import pt2pdg,pdg2pt


def convert2mcpl(filename, readformat):
    """
    Convert particle list with MCPL-compatible format to MCPL format.

    If there is a file with same name as the file to convert, and MCPL
    format, it will be assumed that particle list has already been
    converted.

    Conversion is executed with subprocess.run, calling the corresponding
    conversion executable. For this to work, MCPL binaries must be in PATH.

    Parameters
    ----------
    filename: str
        Name of file with particle list to convert.
    readformat: str
        Particle list format. Valid formats are:
        - 'ssw': MCNP binary surface source.
        - 'phits': PHITS particle list.
        - 'ptrac': MCNP text surface source.
        - 'stock': TRIPOLI-4 stored particles.
        - 'ssv': Text file with space-separated values, with format
                 resulting from 'mcpltool --text' command.

    Returns
    -------
    mcplname: str
        Name of the created or existing MCPL file.
    """
    # Search MCPL file with same name
    already_converted = False
    if readformat == "mcpl":
        already_converted = True
    else:
        filename_base = filename.split('.')[0]
        if os.path.isfile(filename_base+".mcpl.gz"):
            already_converted = True
            filename = filename_base+".mcpl.gz"
        elif os.path.isfile(filename_base+".mcpl"):
            already_converted = True
            filename = filename_base+".mcpl"
    if already_converted:
        print("Using existing file {}".format(filename))
        return filename
    # Must convert
    if not readformat in ["ssw", "phits", "ptrac", "stock", "ssv"]:
        raise Exception("Invalid {} format".format(readformat))
    print("Converting {} file to MCPL...".format(readformat))
    mcplname = filename.split(".")[0]+".mcpl"
    result = subprocess.run(["kstool", readformat+"2mcpl", filename, mcplname],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=True) # if result.returncode != 0: raise exception
    stdout = result.stdout.split()
    mcplname = stdout[stdout.index(b'Created')+1].decode("utf-8")
    print("Done. Created {}".format(mcplname))
    return mcplname

def join2mcpl(filenames, readformat):
    """
    Merge particle lists with MCPL-compatible format into MCPL file.

    Each file is converted to MCPL with convert2mcpl method, and then merged
    with mcpltool command.

    Merged MCPL file name is constructed as the intersection of filenames,
    if any, or set as 'merged.mcpl.gz' otherwise.

    Parameters
    ----------
    filenames: list
        Names of files with particle lists to convert.
    readformat: str
        Particle list format for all files. See convert2mcpl for valid
        formats.

    Returns
    -------
    mergedname: str
        Name of the created MCPL file.
    """
    if np.isscalar(filenames):
        filenames = [filenames]
    mcplnames = []
    for name in filenames:
        mcplnames.append(convert2mcpl(name, readformat))
    if len(mcplnames) == 1:
        return mcplnames[0]
    joinedname = []
    for chars in zip(*mcplnames):
        if len(set(chars))==1:
            joinedname.append(set(chars).pop())
    joinedname = ''.join(joinedname)
    if len(joinedname) == 0:
        joinedname = "joined.mcpl.gz"
    subprocess.run(["kstool", "mcpltool", joinedname, *mcplnames],
                   check=True) # if result.returncode != 0: raise exception
    print("Joined MCPL files into {}".format(joinedname))
    return joinedname

ptmap = {'n':2112, 'p':22, 'e':11, 'e+':-11}

def savessv(pt, parts, ws, outfile): # Equivalent to convert2ascii (in mcpl.py)
    """
    Save particle list to Space-Separated Values file.

    This function is equivalent to convert2ascii from mcpl, and results in
    the same format as 'mcpltool --text' command. Generated SSV file can be
    converted to MCPL format with convert2mcpl.

    Parameters
    ----------
    pt: str
        Particle type. See pt2pdg for available particle types.
    parts: array-like
        Array of particles. Must have shape (N, 7), with columns ordered as
        in varnames list.
    ws: array-like
        Particle weights.
    outfile: str
        Name of SSV file. If it exist, its content will be overwritten.
    """
    print("Writing particles into SSV file...")
    with open(outfile,'w') as fout:
        fout.write("#MCPL-ASCII\n#ASCII-FORMAT: v1\n#NPARTICLES: %i\n"%len(parts))
        fout.write("#END-HEADER\n")
        fout.write("index     pdgcode               ekin[MeV]                   x[cm]          "
                   +"         y[cm]                   z[cm]                      ux                  "
                   +"    uy                      uz                time[ms]                  weight  "
                   +"                 pol-x                   pol-y                   pol-z  userflags\n")
        pdgcode = ptmap[pt]
        fmtstr="%5i %11i %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g 0x%08x\n"
        for idx,p in enumerate(parts):
            fout.write(fmtstr%(idx,pdgcode,p[0],p[1],p[2],p[3],p[4],p[5],p[6],0,ws[idx],0,0,0,0))
    print("Done. All particles written into {}".format(outfile))

def appendssv(pt, parts, ws, outfile):
    """
    Append particle list to Space-Separated Values file.

    This function uses the same format as savessv, but appends particles
    without writing a header.

    Parameters
    ----------
    pt: str
        Particle type. See pt2pdg for available particle types.
    parts: array-like
        Array of particles. Must have shape (N, 7), with columns ordered as
        in varnames list.
    ws: array-like
        Particle weights.
    outfile: str
        Name of SSV file. If it exist, its content will be overwritten.
    """
    print("Appending particles into SSV file...")
    with open(outfile,'a') as fout:
        pdgcode = ptmap[pt]
        fmtstr="%5i %11i %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g 0x%08x\n"
        for idx,p in enumerate(parts):
            fout.write(fmtstr%(idx,pdgcode,p[0],p[1],p[2],p[3],p[4],p[5],p[6],0,ws[idx],0,0,0,0))
    print("Done. All particles appended into {}".format(outfile))

class PList:
    def __init__(self, filename, readformat="mcpl", pt='n', trasl=None, rot=None, switch_x2z=False, set_params=True):
        """
        Object defining particle list. It is a wrapper for a MCPL file.

        It also allows performing a spatial transformation right after reading
        particles (during get method), which is useful for using particle lists
        in simulations with different reference systems.

        Parameters
        ----------
        filename: str or list
            File or list of files containing the particle list, in any
            MCPL-compatible format.
        readformat: str
            Particle list format. Can be 'mcpl' or any valid value of
            convert2mcpl readformat parameter.
        pt: str
            Particle type. See pt2pdg for available particle types.
        trasl: array-like, optional
            Spatial translation. Default is no translation.
        rot: numpy.ndarray or scipy.spatial.transform.Rotation, optional
            Spatial rotation. Can be a scipy Rotation object, or any array-like
            which can be used to generate a scipy Rotation object (rotation
            vector, quaternion or rotation matrix). Rotation is applied after
            translation.
        switch_x2z: bool
            If true, the following permutation is applied after translation and
            rotation:
                (x, y, z) -> (y, z, x)
        set_params: bool
            Whether to set I (sum of weights) and p2 (sum of squared weights)
            after construction.
        """
        if np.isscalar(filename):
            self.filename = convert2mcpl(filename, readformat) # Convert format to MCPL
        else:
            self.filename = join2mcpl(filename, readformat) # Convert format to MCPL and join
        if trasl is not None:
            trasl = np.array(trasl).reshape(-1)
            if trasl.shape != (3,):
                raise ValueError("Invalid trasl")
        if rot is not None:
            if not isinstance(rot, st.Rotation):
                rot = np.array(rot)
                if rot.shape == (4,):
                    rot = st.Rotation.from_quat(rot)
                elif rot.shape == (3,3):
                    rot = st.Rotation.from_matrix(rot)
                elif rot.shape == (3,):
                    rot = st.Rotation.from_rotvec(rot)
                else:
                    raise ValueError("Invalid rot")
        self.pt = pt
        self.trasl = trasl
        self.rot = rot
        self.x2z = switch_x2z
        self.params_set = False
        if set_params:
            self.set_params()

    def set_params(self):
        """Set parameters I (sum of weights) and p2 (sum of squared weights)."""
        pl = mcpl.MCPLFile(self.filename)
        I = p2 = N = 0
        for pb in pl.particle_blocks:
            mask = np.logical_and(pb.weight>0, pb.pdgcode==pt2pdg(self.pt))
            ws = pb.weight[mask]
            N += len(ws)
            I += ws.sum()
            p2 += (ws**2).sum()
        if N == 0:
            raise Exception("Empty particle list")
        N_eff = I**2 / p2
        print("I = {}\np2 = {}\nN = {}\nN_eff = {}".format(I, p2, N, N_eff))
        self.I = I
        self.p2 = p2
        self.N = N
        self.N_eff = N_eff
        self.params_set = True

    def get(self, N=-1, skip=0):
        """
        Get particles from MCPL file.

        After loading particles, translation, rotation and switch_x2z
        transformations will be applied (if any).

        Parameters
        ----------
        N: int
            Number of particles to read. The real number N' of particles read
            may be lower if end of particle list is reached or there are
            particles with zero weight. -1 to read all particles.
        skip: int
            Number of particles to skip in the list before starting to read.

        Returns
        -------
        [parts,ws]: list
            Array of particles and weights. Particles array will have shape
            (N', 7), with columns ordered as in varnames list.
        """
        if N < 0:
            N = mcpl.MCPLFile(self.filename).nparticles
        pl = mcpl.MCPLFile(self.filename, blocklength=N)
        pl.skip_forward(skip)
        pb = pl.read_block()
        if pb is None:
            return np.zeros(0),np.zeros((0,7))
        parts = np.stack((pb.ekin,pb.x,pb.y,pb.z,pb.ux,pb.uy,pb.uz), axis=1)
        ws = pb.weight
        mask = np.logical_and(ws>0, pb.pdgcode==pt2pdg(self.pt))
        parts = parts[mask]
        ws = ws[mask]
        if self.trasl is not None: # Apply translation
            parts[:,1:4] += self.trasl
        if self.rot is not None: # Apply rotation
            parts[:,1:4] = self.rot.apply(parts[:,1:4])
            parts[:,4:7] = self.rot.apply(parts[:,4:7])
        if self.x2z: # Apply permutation (x,y,z) -> (y,z,x)
            ekin,x,y,z,ux,uy,uz = parts.T
            parts = np.array([ekin,y,z,x,uy,uz,ux]).T
        return parts,ws

    def save(self, pltree):
        """Save PList parameters into XML tree."""
        ET.SubElement(pltree, "pt").text = self.pt
        ET.SubElement(pltree, "mcplname").text = os.path.abspath(self.filename)
        trasl = np.array_str(self.trasl)[1:-1] if self.trasl is not None else ""
        ET.SubElement(pltree, "trasl").text = trasl
        rot = np.array_str(self.rot.as_rotvec())[1:-1] if self.rot is not None else ""
        ET.SubElement(pltree, "rot").text = rot
        ET.SubElement(pltree, "x2z").text = str(int(self.x2z))

    @staticmethod
    def load(pltree):
        """Load parameters from XML tree and build PList."""
        pt = pltree[0].text
        filename = pltree[1].text
        if pltree[2].text: trasl = np.array(pltree[2].text.split(), dtype="float64")
        else: trasl = None
        if pltree[3].text: rot = np.array(pltree[3].text.split(), dtype="float64")
        else: rot = None
        switch_x2z = bool(int(pltree[4].text))
        return PList(filename, pt=pt, trasl=trasl, rot=rot, switch_x2z=switch_x2z)
