#define MINIAUDIO_IMPLEMENTATION
#include "utility/audiowrite.h"
#include "utility/maths.h"

#include <iostream>
#include <string>
#include <vector>
#include <websterFDTD.h>

int main(int, char*[])
{
    std::string path = "websterFDTD.wav";
    float samplerate = 44100;
    float duration = 1;
    std::size_t num_sample = static_cast<int>(samplerate * duration);
    tarte::WebsterFDTD<double> resonator(samplerate, 1.0f);
    tarte::Articulation art;
    art.SetFromVowel(tarte::vowels::u);
    resonator.set_yielding_walls(true);
    resonator.set_radiation(true);
    resonator.set_l0(17e-2);
    resonator.set_time_varying_geometry(false);
    resonator.set_compute_powers(true);
    resonator.SetTargetGeometryFromArticulation(art);

    std::vector<float> samples;
    samples.resize(num_sample);

    // Run a simulation with the default parameters and a dirac impulse as input
    float max_power_error = 0;
    float max_power_exchanged = 0;
    for (int i = 0; i < samples.size(); i++) {
        if (i == 0) {
            resonator.Process(1e-3);
        } else {
            resonator.Process(0);
        }
        samples[i] = resonator.ReadRadiatedPressure();

        if (abs(resonator.ReadPowerTotal()) > max_power_error) {
            max_power_error = abs(resonator.ReadPowerTotal());
        }

        if (abs(resonator.ReadPowerFluidStored()) > max_power_exchanged) {
            max_power_exchanged = abs(resonator.ReadPowerFluidStored());
        }
    }
    std::cout << "Maximum relative error on power balance: " << max_power_error / max_power_exchanged << std::endl;
    // Write an example wav file
    tarte::WriteWav(path, samples, samplerate);
    return 0;
}