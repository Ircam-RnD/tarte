
#include <chrono>
#include <iostream>
#include <string>
#include <utility/ResultsStorage.h>
#include <voice.h>
#include <websterFDTD.h>

int main(int argc, char const* argv[])
{
    auto t1 = std::chrono::high_resolution_clock::now();
    std::string path = "VoiceReadWrite.hdf5";
    if (argc > 1) {
        path = argv[1];
    }
    std::cout << "Running with file path:  " << path << std::endl;

    /*
        1. Read simulation parameters from file, initialize model
    */

    // Opens the file
    ResultsStorage storage(path);
    // ------------- General settings --------------- //
    float sr;
    float duration;
    storage.readAttribute("sr", sr);
    storage.readAttribute("duration", duration);
    std::size_t N_samples = static_cast<int>(sr * duration);

    // Input pressure vector
    Eigen::VectorXd Psub = Eigen::VectorXd::Zero(N_samples);
    storage.readVector("Psub", Psub);

    // General flag
    bool compute_powers;
    storage.readAttribute("computePowers", compute_powers);

    // ------------- Resonator settings --------------- //
    // length
    float l0;
    storage.readAttribute("l0", l0);

    // Flags
    bool yielding_walls, radiation, store_spatial_distributions_vt;
    storage.readAttribute("yieldingWalls", yielding_walls);
    storage.readAttribute("radiation", radiation);
    storage.readAttribute("storeSpatialDistributionsVT", store_spatial_distributions_vt);

    // ------------- Larynx settings ----------------- //
    // Mucle activations
    float a_ct, a_ta, a_lc;
    storage.readAttribute("a_ct", a_ct);
    storage.readAttribute("a_ta", a_ta);
    storage.readAttribute("a_lc", a_lc);

    // Noise ratio
    float noise_ratio;
    storage.readAttribute("noiseRatio", noise_ratio);

    // Recovery coefficient
    float kt;
    storage.readAttribute("kt", kt);

    // Area smoothing
    float epsilon_smooth;
    storage.readAttribute("epsilonSmooth", epsilon_smooth);

    // Sav control parameter
    float lambda_sav;
    storage.readAttribute("lambdaSav", lambda_sav);

    // Storage flags
    bool store_larynx_state, store_pressure_drop, store_glottal_flow, store_epsilon_sav;
    storage.readAttribute("storeLarynxState", store_larynx_state);
    storage.readAttribute("storePressureDrop", store_pressure_drop);
    storage.readAttribute("storeGlottalFlow", store_glottal_flow);
    storage.readAttribute("storeEpsilonSav", store_epsilon_sav);

    // ------------- Model initialization --------------- //
    tarte::Voice<tarte::BodyCoverPair<double>, double> proc(sr, true);
    // Vocal tract parameters
    proc.get_resonator()->set_l0(l0);
    proc.get_resonator()->set_yielding_walls(yielding_walls);
    proc.get_resonator()->set_radiation(radiation);
    proc.get_resonator()->set_compute_powers(compute_powers);
    proc.get_resonator()->set_time_varying_geometry(false);

    // Larynx parameters
    proc.get_vocal_folds()->set_muscles_activation(a_ct, a_ta, a_lc);
    proc.get_vocal_folds()->set_noise_ratio(noise_ratio);
    proc.get_vocal_folds()->set_kt(kt);
    proc.get_vocal_folds()->set_epsilon_smooth(epsilon_smooth);
    proc.set_lambda_sav(lambda_sav);

    // Set geometry
    int articulation_mode;
    storage.readAttribute("articulationMode", articulation_mode);
    tarte::Articulation art;
    float constant_section;
    switch (articulation_mode) {
    case 0: // Constant section (straight tube)
        storage.readAttribute("constantSection", constant_section);
        proc.get_resonator()->SetConstantSection(constant_section);
        std::cout << "Using constant section geometry with S = " << constant_section << " m^2." << std::endl;
        break;

    case 1: // Set from formants (interpolation done in the geometrical space, not in the formant space)
        float F1, F2;
        storage.readAttribute("F1", F1);
        storage.readAttribute("F2", F2);
        art.SetFromFormants(F1, F2);
        proc.get_resonator()->SetTargetGeometryFromArticulation(art);
        std::cout << "Using geometry computed from target F1 = " << F1 << ", F2 = " << F2 << " Hz." << std::endl;
        break;

    default: // Constant section (straight tube)
        storage.readAttribute("constantSection", constant_section);
        proc.get_resonator()->SetConstantSection(constant_section);
        std::cout << "Using constant section geometry with S=" << constant_section << " m^2." << std::endl;
        break;
    }

    /*
        2. Prepare storage
    */
    using Vector = Eigen::VectorXd;
    using Matrix = Eigen::MatrixXd;

    Vector radiated_pressure, input_pressure;
    radiated_pressure = Vector::Zero(N_samples);
    input_pressure = Vector::Zero(N_samples);

    // Resonator
    Matrix density_distribution, velocity_distribution;
    if (store_spatial_distributions_vt) {
        density_distribution = Matrix::Zero(proc.get_resonator()->get_N(), N_samples);
        velocity_distribution = Matrix::Zero(proc.get_resonator()->get_N() - 1, N_samples);
    }

    // Larynx
    Matrix folds_displacement, effective_openings;
    Vector pressure_drop, glottal_flow, epsilon_sav;
    if (store_larynx_state) {
        folds_displacement = Matrix::Zero(proc.get_vocal_folds()->get_N(), N_samples);
        effective_openings = Matrix::Zero(proc.get_vocal_folds()->get_half_N(), N_samples);
    }
    if (store_pressure_drop) {
        pressure_drop = Vector::Zero(N_samples);
    }
    if (store_glottal_flow) {
        glottal_flow = Vector::Zero(N_samples);
    }
    if (store_epsilon_sav) {
        epsilon_sav = Vector::Zero(N_samples);
    }

    // Powers (only for resonator)
    Vector P_stored_fluid, P_stored_fluid_kinetic, P_stored_fluid_potential;
    Vector P_stored_walls;
    Vector P_stored_radiation;
    Vector P_diss_walls, P_diss_radiation;
    Vector P_in;
    Vector P_tot;
    if (compute_powers) {
        P_stored_fluid = Vector::Zero(N_samples);
        P_stored_fluid_kinetic = Vector::Zero(N_samples);
        P_stored_fluid_potential = Vector::Zero(N_samples);
        P_stored_walls = Vector::Zero(N_samples);
        P_stored_radiation = Vector::Zero(N_samples);
        P_diss_walls = Vector::Zero(N_samples);
        P_diss_radiation = Vector::Zero(N_samples);
        P_in = Vector::Zero(N_samples);
        P_tot = Vector::Zero(N_samples);
    }

    /*
        3. Run the simulation
    */

    // Run a simulation with the default parameters and a dirac impulse as input
    for (int i = 0; i < N_samples; i++) {
        proc.Process(Psub[i]);
        radiated_pressure(i) = proc.ReadRadiatedPressure();
        input_pressure(i) = proc.get_resonator()->ReadInputPressure();

        if (store_spatial_distributions_vt) {
            density_distribution.col(i) = proc.get_resonator()->ReadCurrentDensityDistribution();
            velocity_distribution.col(i) = proc.get_resonator()->ReadCurrentVelocityDistribution();
        }

        if (store_larynx_state) {
            folds_displacement.col(i) = proc.ReadFoldsDisplacement();
            effective_openings.col(i) = proc.ReadEffectiveOpenings();
        }
        if (store_pressure_drop) {
            pressure_drop(i) = proc.ReadPressureDrop();
        }
        if (store_glottal_flow) {
            glottal_flow(i) = proc.ReadMeanGlottalFlow();
        }
        if (store_epsilon_sav) {
            epsilon_sav(i) = proc.ReadEpsilonSav();
        }

        if (compute_powers) {
            P_stored_fluid(i) = proc.get_resonator()->ReadPowerStoredFluid();
            P_stored_fluid_kinetic(i) = proc.get_resonator()->ReadPowerStoredFluidKinetic();
            P_stored_fluid_potential(i) = proc.get_resonator()->ReadPowerStoredFluidPotential();
            P_stored_walls(i) = proc.get_resonator()->ReadPowerStoredWalls();
            P_stored_radiation(i) = proc.get_resonator()->ReadPowerStoredRadiation();
            P_diss_walls(i) = proc.get_resonator()->ReadPowerDissipatedWalls();
            P_diss_radiation(i) = proc.get_resonator()->ReadPowerDissipatedRadiation();
            P_in(i) = proc.get_resonator()->ReadPowerExchanged();
            P_tot(i) = proc.get_resonator()->ReadPowerTotal();
        }
    }

    /*
        4. Write results to the file (as well as some attributes)
    */
    storage.writeAttribute("rho0", proc.get_resonator()->get_rho0());
    storage.writeAttribute("c0", proc.get_resonator()->get_c0());

    storage.writeVector("radiatedPressure", radiated_pressure);
    storage.writeVector("inputPressure", input_pressure);

    if (store_spatial_distributions_vt) {
        storage.writeMatrix("densityDistribution", density_distribution);
        storage.writeMatrix("velocityDistribution", velocity_distribution);
    }

    if (store_larynx_state) {
        storage.writeMatrix("foldsDisplacement", folds_displacement);
        storage.writeMatrix("effectiveOpenings", effective_openings);
    }
    if (store_pressure_drop) {
        storage.writeVector("pressureDrop", pressure_drop);
    }
    if (store_glottal_flow) {
        storage.writeVector("glottalFlow", glottal_flow);
    }
    if (store_epsilon_sav) {
        storage.writeVector("epsilonSav", epsilon_sav);
    }

    if (compute_powers) {
        storage.writeVector("PStoredFluid", P_stored_fluid);
        storage.writeVector("PStoredFluidKinetic", P_stored_fluid_kinetic);
        storage.writeVector("PStoredFluidPotential", P_stored_fluid_potential);
        storage.writeVector("PStoredWalls", P_stored_walls);
        storage.writeVector("PStoredRadiation", P_stored_radiation);
        storage.writeVector("PDissWalls", P_diss_walls);
        storage.writeVector("PDissRadiation", P_diss_radiation);
        storage.writeVector("PExchanged", P_in);
        storage.writeVector("Ptot", P_tot);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration<double>(t2 - t1).count();
    storage.writeAttribute("elapsedTime", elapsed_time);
    storage.close();
    return 0;
}