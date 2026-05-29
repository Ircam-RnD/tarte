#define MINIAUDIO_IMPLEMENTATION
#include "articulation.h"
#include "larynx.h"
#include "utility/audiowrite.h"
#include "utility/maths.h"
#include "vowels.h"

#include <iostream>
#include <string>
#include <vector>

int main(int, char*[])
{
    std::string path = "larynx.wav";
    float samplerate = 44100;
    float duration = 10;
    float subglottal_pressure = 800;
    std::size_t num_sample = static_cast<int>(samplerate * duration);
    tarte::Larynx<double> proc(samplerate, true);
    tarte::Articulation art;
    art.SetFromVowel(tarte::vowels::i);
    proc.get_resonator()->set_l0(17e-2);
    proc.get_resonator()->set_time_varying_geometry(false);
    proc.get_resonator()->SetTargetGeometryFromArticulation(art);

    proc.set_lambda_sav(10);
    proc.set_contact_stiffness(15);
    // proc.set_noise_ratio(0.3);

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