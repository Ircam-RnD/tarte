import numpy as np
import h5py
import os.path as path
import subprocess

''' 
This class is made to work with the WebsterFDTDreadWrite.cpp file. The main function fills an hdf5 file the relevant
simulation parameters, launch the simulation and returns the results as a python dictionnary.
'''


class WebsterFDTDReadWrite:
    def __init__(self, program_directory):
        self.program_directory = program_directory
        # Initializes with default values for all the necessary parameters

        self.sr = 44100
        self.duration = 0.1
        self.l0 = 17e-2

        self.yielding_walls = True
        self.radiation = True
        self.compute_powers = False
        self.store_spatial_distributions = False

        # Geometry articulation parameters
        self.articulation_mode = 0  # 0 for constant section, 1 for formants setting
        self.constant_section = np.pi * 1e-4  # Only used if articulation_mode == 0
        self.F1 = 300  # Only used if articulation_mode == 1
        self.F2 = 1200  # Only used if articulation_mode == 1

    def Qin(self, t):
        # Overwrite this function to change the input flow rate
        return np.zeros_like(t)

    def run_simulation(self, fname: str, mute=False) -> dict:
        """Write inputs, run the C++ solver, and return the relevant arrays. if mute==True, the C++ code stdout output is ignored."""
        self.N_samples = int(self.sr * self.duration)
        self.t = np.linspace(0, self.duration, self.N_samples)
        self.QinVec = self.Qin(self.t)

        with h5py.File(fname, "w") as f:
            f.attrs["sr"] = self.sr
            f.attrs["duration"] = self.duration
            f["Qin"] = self.QinVec
            f.attrs["l0"] = self.l0

            f.attrs["yieldingWalls"] = self.yielding_walls
            f.attrs["radiation"] = self.radiation
            f.attrs["computePowers"] = self.compute_powers
            f.attrs["storeSpatialDistributions"] = self.store_spatial_distributions

            f.attrs["articulationMode"] = self.articulation_mode
            f.attrs["constantSection"] = self.constant_section
            f.attrs["F1"] = self.F1
            f.attrs["F2"] = self.F2

        if (mute == False):
            subprocess.run(
                [self.program_directory, f"{path.realpath(fname)}"],
                check=True,
                # stdout=subprocess.DEVNULL,
                # stderr=subprocess.DEVNULL,
            )
        else:
            subprocess.run(
                [self.program_directory, f"{path.realpath(fname)}"],
                check=True,
                stdout=subprocess.DEVNULL,
                # stderr=subprocess.DEVNULL,
            )

        with h5py.File(fname, "r") as f:
            results = {}
            results["radiatedPressure"] = f["radiatedPressure"][:]
            if (self.store_spatial_distributions):
                results["densityDistribution"] = f["densityDistribution"][:]
                results["velocityDistribution"] = f["velocityDistribution"][:]
            if (self.compute_powers):
                results["PStoredFluid"] = f["PStoredFluid"][:]
                results["PStoredFluidKinetic"] = f["PStoredFluidKinetic"][:]
                results["PStoredFluidPotential"] = f["PStoredFluidPotential"][:]
                results["PStoredWalls"] = f["PStoredWalls"][:]
                results["PStoredRadiation"] = f["PStoredRadiation"][:]
                results["PDissWalls"] = f["PDissWalls"][:]
                results["PDissRadiation"] = f["PDissRadiation"][:]
                results["PExchanged"] = f["PExchanged"][:]
                results["Ptot"] = f["Ptot"][:]
            return results


if __name__ == "__main__":
    runner = WebsterFDTDReadWrite(
        "../build/examples/WebsterFDTDReadWrite/tarte-websterFDTDReadWrite")
    mute = False

    Amp = 1e-4
    width = 1e-4

    def Qin(t):
        return Amp * np.sin(np.pi * t / width) * (t < width)
    runner.duration = 1
    runner.store_spatial_distributions = False
    runner.compute_powers = False
    runner.Qin = Qin
    results = runner.run_simulation("ReadWriteExample.hdf5", mute)
