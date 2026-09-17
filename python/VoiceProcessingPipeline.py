import pathlib
import sys
import copy
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from functools import partial
import pandas as pd
import numpy as np
from librosa import pyin
import scipy.io.wavfile
import json
from VoiceReadWrite import VoiceReadWrite

"""
This file provides an example of data generation pipeline, with
direct signal processing preventing storage of large files. Multiprocessing
is used.
The processed results are stored as a single hdf5 file exported form a panda dataframe,
enabling easy re-open. Common parameters and general settings used in the data generation are stored
separately in a json file.
"""


class VoiceProcessingPipeline():
    def __init__(self):
        # ----------------- General settings --------------------- #
        self.output_dir = "./results_voice_processed_data"
        self.cpp_runner_directory = "../build/examples/VoiceReadWrite/tarte-VoiceReadWrite"
        # If store_data is set to True, simulation data is deleted along the way, only analysis is kept
        self.store_data = False
        # If store_sounds is set to True, the radiated pressures are exported as wav files
        self.store_sounds = True
        self.N_pools = cpu_count()  # Number of processes in parralel
        self.chunksize = 1  # Multiprocessing chunksize
        self.mute = True  # Mute output of cpp code

        # ----------------- Experimental configuration description --------------------- #
        # This dict should include what you want to be kept alongside the processed data.
        self.experimental_config = {
            "N_sims": 100,
            "Detail":
                "Only subglottal pressure varying, straight tube, and default muscle activation. Fundamental frequency is estimated using Yin.",
            "Varying_Parameters": ["P_sub"],
            "Output_Descriptors": ["Oscillating", "Fundamental_Frequency"]
        }

        # ----------------- Base runner class --------------------- #
        self.base_runner = VoiceReadWrite(self.cpp_runner_directory)
        # Change here all of the default parameters common to all simulations
        self.base_runner.duration = 1

        # ----------------- Varying parameters --------------------- #
        # Change here the values of changing parameters and related functions
        self.P_maxs = np.linspace(0, 2000, self.experimental_config["N_sims"])

    def Psub(self, Pmax, t):
        return Pmax * (t >= 0.05) + Pmax * (t/0.05) * (t < 0.05)

    # ----------------- Processing function --------------------- #
    # Overwrite this function to change the type of processing done on the data and written to
    # the file.

    def process_results(self, results):
        Psub = results["Psub"]
        Pmax = np.max(Psub)
        radiated_pressure = results["radiatedPressure"]
        f0, _, _ = pyin(radiated_pressure[-4096:], fmin=80, fmax=1000,
                        sr=self.base_runner.sr, frame_length=4096)
        f0 = f0[-1]
        amp = np.sum(np.abs(radiated_pressure[-4096:])) / 4096
        oscillates = amp > 0.01
        return {"Pmax": Pmax, "oscillates": oscillates, "f0": f0}, radiated_pressure

    def run_and_process(self, index):
        runner = copy.deepcopy(self.base_runner)
        # Replace functions parameters with those of the current simulation
        runner.Psub = partial(self.Psub, self.P_maxs[index])
        results = runner.run_simulation(
            self.fnames[index], self.mute)
        processed, radiated_pressure = self.process_results(results)
        if self.store_sounds:
            scipy.io.wavfile.write(self.wav_dir.joinpath(
                f"{index}.wav"), rate=self.base_runner.sr, data=radiated_pressure / (np.max(np.abs(radiated_pressure)) + 1e-16))
        if not self.store_data:
            pathlib.Path.unlink(self.fnames[index])
        return processed

    # ----------------- Actual computation --------------------- #
    # It should'nt be needed to mess too much with the following function if you are okay
    # with the way the data is stored.
    def run_pipeline(self):
        # ----------------- Storage setup --------------------- #
        if (pathlib.Path.exists(pathlib.Path(self.output_dir))):
            print(
                "Output directory already exists. Continuing will overwrite data. Do you want to proceed (y/n)?")
            yes = {'yes', 'y', 'ye', ''}
            no = {'no', 'n'}
            choice = input().lower()
            while True:
                if choice in yes:
                    print("Continuing")
                    break
                elif choice in no:
                    sys.exit()
                else:
                    sys.stdout.write("Please respond with 'yes' or 'no'")
        pathlib.Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        # Write the experimental config file
        with open(pathlib.Path(self.output_dir).joinpath("experimental_config.json"), 'w') as f:
            json.dump(self.experimental_config, f, sort_keys=True,
                      indent=4, ensure_ascii=False)
        # Create a folder for storing raw simulation data
        data_dir = pathlib.Path(self.output_dir).joinpath("data")
        pathlib.Path(data_dir).mkdir(parents=True, exist_ok=True)
        # Create a folder for storing wav files if needed
        if self.store_sounds:
            self.wav_dir = pathlib.Path(self.output_dir).joinpath("sounds")
            pathlib.Path(self.wav_dir).mkdir(parents=True, exist_ok=True)

        # ----------------- Run the pipeline --------------------- #
        self.fnames = [pathlib.Path(data_dir).joinpath(
            f"{i}.hdf5") for i in range(self.experimental_config["N_sims"])]

        processed_data_list = []
        with Pool(self.N_pools) as p:
            processed_data_list.append(p.map(self.run_and_process, tqdm([
                i for i in range(self.experimental_config["N_sims"])], total=self.experimental_config["N_sims"]), self.chunksize))

        if not self.store_data:
            pathlib.Path.rmdir(data_dir)
        processed_data = pd.DataFrame(*processed_data_list)
        processed_data.to_csv(pathlib.Path(self.output_dir).joinpath(
            "processed_data.csv"), mode="w")


if __name__ == "__main__":
    pipeline = VoiceProcessingPipeline()
    pipeline.run_pipeline()
