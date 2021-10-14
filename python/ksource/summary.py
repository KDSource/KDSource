#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for collecting McStas and TRIPOLI-4 simulation outputs
"""

import numpy as np
import matplotlib.pyplot as plt

from .plist import PList
from .tally import T4Tally

def read_bashoutput(bashoutput, mccode):
    """
    Process bash output of McStas or TRIPOLI-4 simulation.

    Searches for number of produced particles (with KSource) and
    simulation times.

    Parameters
    ----------
    bashoutput: str
        Name of file containing bash output of simulation.
    mccode: "McStas" or "TRIPOLI"
        MC code used for simulation.
    """
    t_simul = 0
    I_source = 1
    if mccode == "McStas":
        with open(bashoutput, "r") as file:
            for line in file:
                if "KSource" in line:
                    I_source = np.double(line.split()[-3])
                if "Finally" in line:
                    line = line.split()
                    units = line[-1][1:-1]
                    t_simul = np.double(line[-2])
                    if units == "min": t_simul *= 60
                    elif units == "h": t_simul *= 3600
    elif mccode == "TRIPOLI":
        with open(bashoutput, "r") as file:
            for line in file:
                if "Produced particles" in line:
                    I_source = np.double(line.split()[-3])
                if "simulation time" in line:
                    t_simul = np.double(line.split()[-1])
    else:
        raise ValueError("Invalid mccode.")
    return [t_simul, I_source]

class Summary:
    def __init__(self, mccode, folder, bashoutput="bash.out", n_detectors=[], p_detectors=[], t4output=None, tallies=[]):
        """
        Object representing summary of MC simulation.

        After a Monte Carlo simulation, this class helps collecting the
        results, such as recorded particle lists or tallies, and
        simulation times. It also computes integral magnitudes of said
        results.

        Parameters
        ----------
        mccode: "McStas" or "TRIPOLI"
            MC code used for simulation.
        folder: str
            Directory containing simulation output files.
        bashoutput: str
            Name of file containing bash output of simulation.
        n_detectors: list, optional
            Neutron lists recorded.
        p_detectors: list, optional
            Photon lists recorded.
        t4outup: str, optional
            Name of TRIPOLI output file. Ignore if mccode is not
            TRIPOLI.
        tallies: list, optional
            Tally names recorded in TRIPOLI. Ignore if mccode is not
            TRIPOLI.
        """
        if any(mccode == mc for mc in ["McStas", "TRIPOLI"]):
            self.mccode = mccode
        else:
            raise ValueError("Invalid mccode.")
        self.folder = folder
        self.bashoutput = bashoutput
        self.t4output = t4output
        self.n_detectors = n_detectors
        self.p_detectors = p_detectors
        self.tallies = tallies
        self.initialized = False
    def compute(self):
        """
        Load results and compute integral magnitudes.
        """
        if self.bashoutput is not None:
            self.t_simul, self.I_source = read_bashoutput(self.folder+"/"+self.bashoutput, self.mccode)
        else:
            self.t_simul = 0
            self.I_source = 1
        self.n_det_scores = []
        for detector in self.n_detectors:
            if self.mccode == "McStas": readformat = "mcpl"
            if self.mccode == "TRIPOLI": readformat = "stock"
            plist = PList(self.folder+"/"+detector, readformat)
            self.n_det_scores.append([plist.I, np.sqrt(plist.p2)])
        self.p_det_scores = []
        for detector in self.p_detectors:
            if self.mccode == "McStas": readformat = "mcpl"
            if self.mccode == "TRIPOLI": readformat = "stock"
            plist = PList(self.folder+"/"+detector, readformat)
            self.p_det_scores.append([plist.I, np.sqrt(plist.p2)])
        self.tally_scores = []
        if self.t4output is not None:
            for tallyname in self.tallies:
                tally = T4Tally(self.folder+"/"+self.t4output, tallyname)
                self.tally_scores.append([self.I_source*np.sum(tally.I), self.I_source*np.sqrt(np.sum(tally.err**2))])
        self.initialized = True
    def save(self, filename):
        """
        Save results in text file, able to be imported to spreadsheet.

        Parameters
        ----------
        filename: str
            Name of file to save results.
        """
        if not self.initialized:
            print("Must compute results before saving.")
            return
        with open(self.folder+"/"+filename, "w") as file:
            file.write("t_simul\t{}\n".format(self.t_simul))
            file.write("I_source\t{}\n".format(self.I_source))
            file.write("n_detectors:\n")
            for det in self.n_detectors: file.write(det.split(sep='.')[0]+"\t\t")
            file.write("\n")
            np.savetxt(file, np.reshape(self.n_det_scores, (1,-1)))
            file.write("p_detectors:\n")
            for det in self.p_detectors: file.write(det.split(sep='.')[0]+"\t\t")
            file.write("\n")
            np.savetxt(file, np.reshape(self.p_det_scores, (1,-1)))
            file.write("tallies:\n")
            for tallyname in self.tallies: file.write(tallyname+"\t\t")
            file.write("\n")
            np.savetxt(file, np.reshape(self.tally_scores, (1,-1)))
        print("Summary successfully saved.")