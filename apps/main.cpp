#include <iostream>
#include <iomanip>
#include <string>

#include "gem/core/types.hpp"
#include "gem/core/constants.hpp"
#include "gem/core/config.hpp"
#include "gem/hadrons/meson.hpp"
#include "gem/hadrons/baryon.hpp"

using namespace gem;

void printHeader() {
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "  GEM Quark Model Calculator\n";
    std::cout << "  Gaussian Expansion Method for\n";
    std::cout << "  Multiquark State Masses\n";
    std::cout << "========================================\n\n";
}

void printMesonResult(const std::string& name, Meson& meson, Real exp_mass_mev) {
    Real mass_gev = meson.calculateMass();
    Real mass_mev = mass_gev * constants::GEV_TO_MEV;  // Convert to MeV for display
    Real diff = mass_mev - exp_mass_mev;
    Real percent = (diff / exp_mass_mev) * 100.0;

    std::cout << std::setw(12) << name << ": "
              << std::fixed << std::setprecision(1)
              << std::setw(8) << mass_mev << " MeV "
              << "(exp: " << std::setw(6) << exp_mass_mev << " MeV, "
              << "diff: " << std::showpos << std::setw(6) << diff << " MeV, "
              << std::setw(5) << percent << "%)" << std::noshowpos;

    if (meson.isValid()) {
        std::cout << "  [OK]";
    } else {
        std::cout << "  [FAILED]";
    }
    std::cout << "\n";
}

void calculateMesons() {
    std::cout << "=== Meson Masses ===\n\n";
    std::cout << "Light Mesons:\n";

    // Pion (π): ud̄, S=0
    Meson pion = MesonFactory::pion();
    printMesonResult("pion (π)", pion, 140.0);

    // Rho (ρ): ud̄, S=1
    Meson rho = MesonFactory::rho();
    printMesonResult("rho (ρ)", rho, 775.0);

    // Kaon (K): us̄, S=0
    Meson kaon = MesonFactory::kaon();
    printMesonResult("kaon (K)", kaon, 494.0);

    // K* (K*): us̄, S=1
    Meson kstar = MesonFactory::kstar();
    printMesonResult("K* (K*)", kstar, 892.0);

    std::cout << "\nHeavy-Light Mesons:\n";

    // D meson: cd̄, S=0
    Meson dmeson = MesonFactory::dmeson();
    printMesonResult("D meson", dmeson, 1870.0);

    // D* meson: cd̄, S=1
    Meson dstar = MesonFactory::dstar();
    printMesonResult("D* meson", dstar, 2010.0);

    std::cout << "\nHeavy Quarkonia:\n";

    // η_c: cc̄, S=0
    Meson etac = MesonFactory::etac();
    printMesonResult("eta_c (ηc)", etac, 2984.0);

    // J/ψ: cc̄, S=1
    Meson jpsi = MesonFactory::jpsi();
    printMesonResult("J/psi (J/ψ)", jpsi, 3097.0);

    // η_b: bb̄, S=0
    Meson etab = MesonFactory::etab();
    printMesonResult("eta_b (ηb)", etab, 9399.0);

    // Υ: bb̄, S=1
    Meson upsilon = MesonFactory::upsilon();
    printMesonResult("Upsilon (Υ)", upsilon, 9460.0);

    std::cout << "\n";
}

void calculateBaryons() {
    std::cout << "=== Baryon Masses ===\n\n";

    // Proton: uud, S=1/2
    Baryon proton = BaryonFactory::proton();
    Real proton_mass = proton.calculateMass() * constants::GEV_TO_MEV;
    std::cout << std::setw(12) << "proton (p)" << ": "
              << std::fixed << std::setprecision(1)
              << std::setw(8) << proton_mass << " MeV "
              << "(exp: 938 MeV)";
    if (proton.isValid()) {
        std::cout << "  [OK]";
    } else {
        std::cout << "  [FAILED]";
    }
    std::cout << "\n";

    // Neutron: udd, S=1/2
    Baryon neutron = BaryonFactory::neutron();
    Real neutron_mass = neutron.calculateMass() * constants::GEV_TO_MEV;
    std::cout << std::setw(12) << "neutron (n)" << ": "
              << std::fixed << std::setprecision(1)
              << std::setw(8) << neutron_mass << " MeV "
              << "(exp: 940 MeV)";
    if (neutron.isValid()) {
        std::cout << "  [OK]";
    } else {
        std::cout << "  [FAILED]";
    }
    std::cout << "\n";

    // Delta: uuu, S=3/2
    Baryon delta = BaryonFactory::delta_pp();
    Real delta_mass = delta.calculateMass() * constants::GEV_TO_MEV;
    std::cout << std::setw(12) << "Delta (Δ++)" << ": "
              << std::fixed << std::setprecision(1)
              << std::setw(8) << delta_mass << " MeV "
              << "(exp: 1232 MeV)";
    if (delta.isValid()) {
        std::cout << "  [OK]";
    } else {
        std::cout << "  [FAILED]";
    }
    std::cout << "\n";

    // Lambda: uds, S=1/2
    Baryon lambda = BaryonFactory::lambda();
    Real lambda_mass = lambda.calculateMass() * constants::GEV_TO_MEV;
    std::cout << std::setw(12) << "Lambda (Λ)" << ": "
              << std::fixed << std::setprecision(1)
              << std::setw(8) << lambda_mass << " MeV "
              << "(exp: 1116 MeV)";
    if (lambda.isValid()) {
        std::cout << "  [OK]";
    } else {
        std::cout << "  [FAILED]";
    }
    std::cout << "\n\n";
}

void printUsage(const char* program) {
    std::cout << "Usage: " << program << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  --meson       Calculate meson masses only\n";
    std::cout << "  --baryon      Calculate baryon masses only\n";
    std::cout << "  --all         Calculate all hadron masses (default)\n";
    std::cout << "  --help        Show this help message\n";
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    bool calc_mesons = true;
    bool calc_baryons = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "--meson") {
            calc_mesons = true;
            calc_baryons = false;
        } else if (arg == "--baryon") {
            calc_mesons = false;
            calc_baryons = true;
        } else if (arg == "--all") {
            calc_mesons = true;
            calc_baryons = true;
        }
    }

    printHeader();

    std::cout << "Using default model parameters:\n";
    std::cout << "  alpha_s (Coulomb) = " << constants::default_potential::ALPHA_S_COULOMB << "\n";
    std::cout << "  alpha_s (hyperfine) = " << constants::default_potential::ALPHA_S_HYPERFINE << "\n";
    std::cout << "  b (string tension) = " << constants::default_potential::B_STRING << " GeV^2\n";
    std::cout << "  V0 = " << constants::default_potential::V0 << " GeV\n";
    std::cout << "  Basis size = " << constants::default_basis::N_BASIS << "\n";
    std::cout << "  r range = [" << constants::default_basis::R_MIN << ", "
              << constants::default_basis::R_MAX << "] fm\n\n";

    if (calc_mesons) {
        calculateMesons();
    }

    if (calc_baryons) {
        calculateBaryons();
    }

    std::cout << "========================================\n";
    std::cout << "  Calculation complete\n";
    std::cout << "========================================\n\n";

    return 0;
}
