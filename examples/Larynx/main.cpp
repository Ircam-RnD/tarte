#define MINIAUDIO_IMPLEMENTATION
#include "utility/audiowrite.h"
#include "utility/maths.h"

#include <iostream>
#include <larynx.h>
#include <string>
#include <vector>

int main(int, char*[])
{
    std::string path = "larynx.wav";
    float samplerate = 44100;
    float duration = 10;
    float subglottal_pressure = 300;
    std::size_t num_sample = static_cast<int>(samplerate * duration);
    tarte::Larynx<double> proc(samplerate, true);
    proc.get_resonator()->set_l0(17e-2);
    proc.get_resonator()->SetConstantSection(25e-4);

    std::vector<float> samples;
    samples.resize(num_sample);

    // Run a simulation with the default parameters and a dirac impulse as input
    for (int i = 0; i < samples.size(); i++) {
        proc.Process(subglottal_pressure);
        samples[i] = proc.ReadRadiatedPressure();
    }
    tarte::normalize_vector(samples);
    // Write an example wav file
    tarte::WriteWav(path, samples, samplerate);
    return 0;
}