#define MINIAUDIO_IMPLEMENTATION
#include "utility/audiowrite.h"
#include "utility/maths.h"

#include <chrono>
#include <iostream>
#include <string>
#include <utility/ResultsStorage.h>
#include <vector>
#include <websterFDTD.h>

int main(int argc, char const* argv[])
{
    auto t1 = std::chrono::high_resolution_clock::now();
    std::string path = "websterFDTDReadWrite.hdf5";
    if (argc > 1) {
        path = argv[1];
    }
    std::cout << "Running with file path:  " << path << std::endl;

    /*
        1. Read simulation parameters from file, initialize model
    */

    // Opens the file
    ResultsStorage storage(path);
    // Samplerate, duration
    float sr;
    float duration;
    storage.readAttribute("sr", sr);
    storage.readAttribute("duration", duration);
    std::size_t N_samples = static_cast<int>(sr * duration);

    // Input flow vector
    Eigen::VectorXd Qin = Eigen::VectorXd::Zero(N_samples);
    storage.readVector("Qin", Qin);

    // length
    float l0;
    storage.readAttribute("l0", l0);

    // Flags
    bool yielding_walls, radiation, compute_powers;
    storage.readAttribute("yieldingWalls", yielding_walls);
    storage.readAttribute("radiation", radiation);
    storage.readAttribute("computePowers", compute_powers);

    // Initialize model
    tarte::WebsterFDTD<double> proc(sr);
    proc.set_l0(l0);
    proc.set_yielding_walls(yielding_walls);
    proc.set_radiation(radiation);
    proc.set_compute_powers(compute_powers);
    proc.set_time_varying_geometry(false);

    // Set geometry
    int articulation_mode;
    storage.readAttribute("articulationMode", articulation_mode);
    tarte::Articulation art;
    float constant_section;
    switch (articulation_mode) {
    case 0: // Constant section (straight tube)
        storage.readAttribute("constantSection", constant_section);
        proc.SetConstantSection(constant_section);
        std::cout << "Using constant section geometry with S = " << constant_section << " m^2." << std::endl;
        break;

    case 1: // Set from formants (interpolation done in the geometrical space, not in the formant space)
        float F1, F2;
        storage.readAttribute("F1", F1);
        storage.readAttribute("F2", F2);
        art.SetFromFormants(F1, F2);
        proc.SetTargetGeometryFromArticulation(art);
        std::cout << "Using geometry computed from target F1 = " << F1 << ", F2 = " << F2 << " Hz." << std::endl;
        break;

    default: // Constant section (straight tube)
        storage.readAttribute("constantSection", constant_section);
        proc.SetConstantSection(constant_section);
        std::cout << "Using constant section geometry with S=" << constant_section << " m^2." << std::endl;
        break;
    }

    /*
        2. Prepare storage
    */

    std::vector<float> radiated_pressure;
    radiated_pressure.resize(N_samples);

    /*
        3. Run the simulation
    */

    // Run a simulation with the default parameters and a dirac impulse as input
    for (int i = 0; i < N_samples; i++) {
        proc.Process(Qin[i]);
        radiated_pressure[i] = proc.ReadRadiatedPressure();
    }

    /*
        4. Write results to the file
    */

    storage.writeVector("radiatedPressure", radiated_pressure);

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration<double>(t2 - t1).count();
    storage.writeAttribute("elapsedTime", elapsed_time);
    storage.close();
    return 0;
}