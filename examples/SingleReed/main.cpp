#define MINIAUDIO_IMPLEMENTATION
#include "articulation.h"
#include "single_reed.h"
#include "utility/audiowrite.h"
#include "utility/maths.h"

#include <iostream>
#include <string>
#include <vector>

int main(int, char*[])
{
    std::string path = "single_reed.wav";
    float samplerate = 44100;
    float duration = 10;
    float mouth_pressure = 2000;
    std::size_t num_sample = static_cast<int>(samplerate * duration);
    tarte::SingleReed<double> proc(samplerate);
    tarte::Articulation art;
    proc.get_resonator()->set_l0(0.66);
    proc.get_resonator()->set_time_varying_geometry(false);
    proc.get_resonator()->SetConstantSection(M_PI * 1e-4);

    std::vector<float> samples;
    samples.resize(num_sample);

    // Run a simulation with the default parameters and a dirac impulse as input
    for (int i = 0; i < samples.size(); i++) {
        proc.Process(mouth_pressure);
        samples[i] = proc.ReadRadiatedPressure();
    }
    tarte::normalize_vector(samples);
    // Write an example wav file
    tarte::WriteWav(path, samples, samplerate);
    return 0;
}