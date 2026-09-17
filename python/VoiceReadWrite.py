import numpy as np
import h5py
import os.path as path
import subprocess

''' 
This class is made to work with the VoiceReadWrite.cpp file. The main function fills an hdf5 file the relevant
simulation parameters, launch the simulation and returns the results as a python dictionnary.
'''


class VoiceReadWrite:
    def __init__(self, program_directory):
        self.program_directory = program_directory
        # Initializes with default values for all the necessary parameters

        self.sr = 44100
        self.duration = 0.1
        self.l0 = 17e-2

        self.yielding_walls = True
        self.radiation = True
        self.compute_powers = False
        # Storing spatial data makes the files big: caution
        self.store_spatial_distributions_vt = False

        # Geometry articulation parameters
        self.articulation_mode = 0  # 0 for constant section, 1 for formants setting
        self.constant_section = np.pi * 0.25e-4  # Only used if articulation_mode == 0
        self.F1 = 300  # Only used if articulation_mode == 1
        self.F2 = 1200  # Only used if articulation_mode == 1

        # Larynx parameters
        self.a_ct = 0.
        self.a_ta = 0.
        self.a_lc = 0.49

        self.noise_ratio = 0
        self.kt = 1.3
        self.epsilon_smooth = 1e-5
        self.lambda_sav = 0

        self.store_larynx_state = False
        self.store_pressure_drop = False
        self.store_glottal_flow = False
        self.store_epsilon_sav = False

    def Psub(self, t):
        # Overwrite this function to change the subglottal pressure
        return np.zeros_like(t)

    def run_simulation(self, fname: str, mute=False) -> dict:
        """Write inputs, run the C++ solver, and return the relevant arrays. if mute==True, the C++ code stdout output is ignored."""
        self.N_samples = int(self.sr * self.duration)
        self.t = np.linspace(0, self.duration, self.N_samples)
        self.PinVec = self.Psub(self.t)

        with h5py.File(fname, "w") as f:
            f.attrs["sr"] = self.sr
            f.attrs["duration"] = self.duration
            f["Psub"] = self.PinVec
            f["t"] = self.t

            f.attrs["l0"] = self.l0

            f.attrs["yieldingWalls"] = self.yielding_walls
            f.attrs["radiation"] = self.radiation
            f.attrs["computePowers"] = self.compute_powers
            f.attrs["storeSpatialDistributionsVT"] = self.store_spatial_distributions_vt

            f.attrs["articulationMode"] = self.articulation_mode
            f.attrs["constantSection"] = self.constant_section
            f.attrs["F1"] = self.F1
            f.attrs["F2"] = self.F2

            # Larynx parameters
            f.attrs["a_ct"] = self.a_ct
            f.attrs["a_ta"] = self.a_ta
            f.attrs["a_lc"] = self.a_lc

            f.attrs["noiseRatio"] = self.noise_ratio
            f.attrs["kt"] = self.kt
            f.attrs["epsilonSmooth"] = self.epsilon_smooth
            f.attrs["lambdaSav"] = self.lambda_sav

            f.attrs["storeLarynxState"] = self.store_larynx_state
            f.attrs["storePressureDrop"] = self.store_pressure_drop
            f.attrs["storeGlottalFlow"] = self.store_glottal_flow
            f.attrs["storeEpsilonSav"] = self.store_epsilon_sav

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
            results["rho0"] = f.attrs["rho0"]
            results["c0"] = f.attrs["c0"]
            results["Psub"] = self.PinVec
            results["radiatedPressure"] = f["radiatedPressure"][:]
            results["inputPressure"] = f["inputPressure"][:]
            if (self.store_spatial_distributions_vt):
                results["densityDistribution"] = f["densityDistribution"][:]
                results["velocityDistribution"] = f["velocityDistribution"][:]
            if (self.store_larynx_state):
                results["foldsDisplacement"] = f["foldsDisplacement"][:]
                results["effectiveOpenings"] = f["effectiveOpenings"][:]
            if (self.store_glottal_flow):
                results["glottalFlow"] = f["glottalFlow"][:]
            if (self.store_pressure_drop):
                results["pressureDrop"] = f["pressureDrop"][:]

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
    runner = VoiceReadWrite(
        "../build/examples/VoiceReadWrite/tarte-VoiceReadWrite")
    mute = False
    runner.lambda_sav = 1000

    Pmax = 400
    trise = 0.05

    def Psub(t):
        return Pmax * (t >= trise) + Pmax * (t/trise) * (t < trise)
    runner.duration = 1
    runner.store_spatial_distributions = True
    runner.compute_powers = True
    runner.store_pressure_drop = True
    runner.store_glottal_flow = True
    runner.store_larynx_state = True
    runner.store_epsilon_sav = True
    runner.Psub = Psub
    results = runner.run_simulation("VoiceReadWriteExample2.hdf5", mute)
