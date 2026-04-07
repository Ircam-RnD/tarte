#define MINIAUDIO_IMPLEMENTATION
#include "utility/audiowrite.h"
#include "utility/vectormaths.h"

#include <iostream>
#include <string>
#include <vector>
#include <websterFDTD.h>

int main(int, char*[])
{
    std::string path = "websterFDTD.wav";
    float samplerate = 44100;
    float duration = 10;
    std::size_t num_sample = static_cast<int>(samplerate * duration);
    tarte::WebsterFDTD<float> resonator(samplerate, 1.0f);
    resonator.set_yielding_walls(false);
    resonator.SetConstantSection(1e-4);

    std::vector<float> samples;
    samples.resize(num_sample);

    // Run a simulation with the default parameters and a dirac impulse as input
    for (int i = 0; i < samples.size(); i++) {
        resonator.Process((((double)rand() / (RAND_MAX)) - 0.5) * 1e-4);
        samples[i] = resonator.ReadRadiatedPressure();
    }
    utility::normalize_vector(samples);
    // Write an example wav file
    utility::WriteWav(path, samples, samplerate);
    return 0;
}