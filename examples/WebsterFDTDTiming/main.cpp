#define MINIAUDIO_IMPLEMENTATION
#include "utility/audiowrite.h"
#include "utility/maths.h"

#include <chrono>
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
    tarte::WebsterFDTD<float, 300> resonator(samplerate, 1.0f);
    resonator.set_yielding_walls(false);
    resonator.SetConstantSection(1e-4);
    resonator.set_l0(1.0);

    std::cout << "Number of segment is " << resonator.get_N() << std::endl;

    std::vector<float> samples;
    samples.resize(num_sample);

    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < samples.size(); i++) {
        resonator.Process((((double)rand() / (RAND_MAX)) - 0.5) * 1e-4);
        samples[i] = resonator.ReadRadiatedPressure();
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Constant geometry, no yielding walls : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / duration
              << " milliseconds of computation per second of output\n";

    resonator.set_yielding_walls(true);
    t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < samples.size(); i++) {
        resonator.Process((((double)rand() / (RAND_MAX)) - 0.5) * 1e-4);
        samples[i] = resonator.ReadRadiatedPressure();
    }
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Constant geometry, yielding walls : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / duration
              << " milliseconds of computation per second of output\n";

    resonator.set_time_varying_geometry(true);
    t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < samples.size(); i++) {
        resonator.Process((((double)rand() / (RAND_MAX)) - 0.5) * 1e-4);
        samples[i] = resonator.ReadRadiatedPressure();
    }
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Varying geometry, yielding walls : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / duration
              << " milliseconds of computation per second of output\n";
    return 0;
}