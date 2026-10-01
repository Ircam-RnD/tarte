from smt.surrogate_models import KPLS
import pandas as pd
import pathlib
import numpy as np
import matplotlib.pyplot as plt


def train(data, param_names, descriptor="f0"):
    ''' Trains a surrogate model to predict the descriptor as a function of the parameters. 
    Data must be a panda dataframe. Data masking must be done beforehand if needed.'''
    sm = KPLS(theta0=[1], eval_noise=True, hyper_opt="Cobyla")
    sm.options["print_global"] = False
    sm.options["print_training"] = False
    sm.options["print_prediction"] = False

    # Mask non-oscillating points before training
    masked_results = data.mask(
        data['oscillates'] == False).dropna(axis="rows")
    xt = masked_results[param_names].to_numpy()
    yt = masked_results[descriptor].to_numpy()
    # masked_results.np
    sm.set_training_values(xt, yt)
    sm.train()
    return sm


if __name__ == "__main__":
    result_dir = "./results_voice_processed_data"

    df = pd.read_csv(pathlib.Path(
        result_dir).joinpath("processed_data.csv"))

    masked_results = df["oscillates"]
    # df.mask(
    #     df['oscillates'] == False).dropna(axis="rows")

    sm = train(df, ["Pmax"], "f0")

    plt.Figure()
    Ptests = np.linspace(0, 2000, 500)
    plt.plot(Ptests, sm.predict_values(Ptests))
    plt.scatter(df["Pmax"], df["f0"])
    plt.xlabel(r"max($P_{sub}$) [Pa]")
    plt.ylabel("f0 [Hz]")
    plt.grid()
    plt.show()
    # sm.save(join(result_dir, "kriging.bin"))
