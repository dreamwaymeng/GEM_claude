#ifndef GEM_CORE_CONSTANTS_HPP
#define GEM_CORE_CONSTANTS_HPP

#include "types.hpp"
#include <cmath>

namespace gem {
namespace constants {

// Mathematical constants
constexpr Real PI = 3.14159265358979323846;
constexpr Real SQRT_PI = 1.7724538509055160273;
constexpr Real PI_SQUARED = 9.8696044010893586188;

// Physical constants - ALL IN GeV UNITS
constexpr Real HBAR_C = 0.1973269804;  // GeV·fm (ℏc)
constexpr Real HBAR_C_SQUARED = 0.038937929;  // (GeV·fm)^2

// Unit conversions
constexpr Real GEV_TO_MEV = 1000.0;
constexpr Real MEV_TO_GEV = 0.001;
constexpr Real FM_TO_GEV_INV = 1.0 / 0.1973269804;  // fm to GeV^-1
constexpr Real GEV_INV_TO_FM = 0.1973269804;        // GeV^-1 to fm

// Quark masses (in GeV) - constituent masses from AL1 model
namespace quark_mass {
    constexpr Real UP = 0.315;      // GeV
    constexpr Real DOWN = 0.315;    // GeV (same as up)
    constexpr Real STRANGE = 0.577; // GeV
    constexpr Real CHARM = 1.836;   // GeV
    constexpr Real BOTTOM = 5.227;  // GeV
    constexpr Real TOP = 173.0;     // GeV - Not used in hadron spectroscopy
}

// Current quark masses (in GeV) - for reference
namespace current_quark_mass {
    constexpr Real UP = 0.0022;
    constexpr Real DOWN = 0.0047;
    constexpr Real STRANGE = 0.096;
    constexpr Real CHARM = 1.28;
    constexpr Real BOTTOM = 4.18;
    constexpr Real TOP = 173.1;
}

// Color factors
namespace color_factor {
    // <λ_i · λ_j> / 4 for different configurations
    constexpr Real MESON_QQ_BAR = -4.0 / 3.0;     // q-q̄ in color singlet
    constexpr Real BARYON_QQ = -2.0 / 3.0;        // q-q in antisymmetric 3̄

    // Casimir operators
    constexpr Real C_F = 4.0 / 3.0;  // Fundamental representation
    constexpr Real C_A = 3.0;        // Adjoint representation
    constexpr Real T_F = 0.5;        // Index of fundamental

    // Number of colors
    constexpr int N_C = 3;
}

// Spin factors
namespace spin_factor {
    // <s_i · s_j> for different spin states
    constexpr Real SPIN_0 = -3.0 / 4.0;  // S=0: (0-3/4)/1 = -3/4
    constexpr Real SPIN_1 = 1.0 / 4.0;   // S=1: (2-3/4)/3 = 1/4

    // For S(S+1)/2 - 3/4
    inline Real spin_spin_expectation(int twice_S) {
        Real S = static_cast<Real>(twice_S) / 2.0;
        return S * (S + 1) / 2.0 - 0.75;
    }
}

// Default potential parameters (AL1 model) - ALL IN GeV
namespace default_potential {
    // Different alpha_s for Coulomb and hyperfine
    constexpr Real ALPHA_S_COULOMB = 0.380175;    // Coulomb coupling (dimensionless)
    constexpr Real ALPHA_S_HYPERFINE = 1.39568;   // Hyperfine coupling (dimensionless)

    // Hyperfine smearing: tau = mu^b_exp / a, where mu is reduced mass in GeV
    constexpr Real HYPERFINE_A = 1.6553;          // Smearing parameter a
    constexpr Real HYPERFINE_B_EXP = 0.2204;      // Mass exponent b

    // Confinement
    constexpr Real B_STRING = 0.1653;             // String tension (GeV^2)
    constexpr Real V0 = 0.624075;                 // Constant shift (GeV)

    // Legacy compatibility
    constexpr Real ALPHA_S = ALPHA_S_COULOMB;
    constexpr Real SIGMA_HYPERFINE = 1.0;         // Not used directly anymore
}

// Basis parameters defaults - using r_min/r_max in fm
namespace default_basis {
    constexpr int N_BASIS = 20;
    constexpr Real R_MIN = 0.0001; // fm (smallest radius -> largest nu)
    constexpr Real R_MAX = 3.0;    // fm (largest radius -> smallest nu)
    // nu = 1/r^2 in fm^-2
    // nu_min = 1/(3.0)^2 = 0.111 fm^-2
    // nu_max = 1/(0.0001)^2 = 10^8 fm^-2
}

// Helper function to get quark mass (returns GeV)
inline Real getQuarkMass(QuarkFlavor flavor) {
    switch (flavor) {
        case QuarkFlavor::Up:      return quark_mass::UP;
        case QuarkFlavor::Down:    return quark_mass::DOWN;
        case QuarkFlavor::Strange: return quark_mass::STRANGE;
        case QuarkFlavor::Charm:   return quark_mass::CHARM;
        case QuarkFlavor::Bottom:  return quark_mass::BOTTOM;
        case QuarkFlavor::Top:     return quark_mass::TOP;
        default: return quark_mass::UP;
    }
}

// Get quark flavor name
inline const char* getQuarkName(QuarkFlavor flavor) {
    switch (flavor) {
        case QuarkFlavor::Up:      return "u";
        case QuarkFlavor::Down:    return "d";
        case QuarkFlavor::Strange: return "s";
        case QuarkFlavor::Charm:   return "c";
        case QuarkFlavor::Bottom:  return "b";
        case QuarkFlavor::Top:     return "t";
        default: return "?";
    }
}

// Convert r (fm) to nu (GeV^2): nu = (hbar_c / r)^2
inline Real r_to_nu(Real r_fm) {
    return (HBAR_C / r_fm) * (HBAR_C / r_fm);
}

// Convert nu (GeV^2) to r (fm): r = hbar_c / sqrt(nu)
inline Real nu_to_r(Real nu_gev2) {
    return HBAR_C / std::sqrt(nu_gev2);
}

} // namespace constants
} // namespace gem

#endif // GEM_CORE_CONSTANTS_HPP
