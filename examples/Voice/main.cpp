#define MINIAUDIO_IMPLEMENTATION
#include "articulation.h"
#include "utility/audiowrite.h"
#include "utility/maths.h"
#include "vocal_folds/body_cover_pair.h"
#include "voice.h"
#include "vowels.h"

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

int main(int, char*[])
{
    std::string path = "voice.wav";
    float samplerate = 96000;
    float duration = 1.0;
    float subglottal_pressure = 800;
    std::size_t num_sample = static_cast<int>(samplerate * duration);
    tarte::Voice<tarte::BodyCoverPair<double>, double> proc(samplerate, true);
    tarte::Articulation art;
    art.SetFromVowel(tarte::vowels::a);
    proc.get_resonator()->set_l0(17e-2);
    proc.get_resonator()->set_time_varying_geometry(false);
    proc.get_resonator()->SetTargetGeometryFromArticulation(art);

    proc.set_lambda_sav(1000);
    // proc.set_contact_stiffness(15);
    // proc.set_eta_stiffness(0);
    // proc.set_noise_ratio(0.3);

    // tarte::VFPair<tarte::BodyCoverVF<double>, double> vf_pair;

    std::vector<float> samples;
    samples.resize(num_sample);

    auto t1 = std::chrono::high_resolution_clock::now();
    // Run a simulation with the default parameters and a dirac impulse as input
    for (int i = 0; i < samples.size(); i++) {
        proc.Process(subglottal_pressure);
        samples[i] = proc.ReadRadiatedPressure();
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Constant geometry, yielding walls : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / duration
              << " milliseconds of computation per second of output\n";
    tarte::normalize_vector(samples);
    // Write an example wav file
    tarte::WriteWav(path, samples, samplerate);
    return 0;
}