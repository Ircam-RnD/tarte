#define MINIAUDIO_IMPLEMENTATION
#include "utility/audiowrite.h"
#include "utility/maths.h"

#include <duffing.h>
#include <iostream>
#include <string>
#include <vector>

int main(int, char*[])
{
    std::string path = "duffing.wav";
    float samplerate = 44100;
    float duration = 1;
    std::size_t num_sample = static_cast<int>(samplerate * duration);
    tarte::Duffing<float> duffing_oscillator(samplerate);
    duffing_oscillator.set_nonlinearity(1e6);

    std::vector<float> samples;
    samples.resize(num_sample);

    // Run a simulation with the default parameters and a dirac impulse as input
    for (int i = 0; i < samples.size(); i++) {
        if (i == 0)
            duffing_oscillator.Process(1.0f);
        else
            duffing_oscillator.Process(0.0f);
        samples[i] = duffing_oscillator.ReadDisplacement();
    }
    tarte::normalize_vector(samples);
    // Write an example wav file
    tarte::WriteWav(path, samples, samplerate);
    return 0;
}