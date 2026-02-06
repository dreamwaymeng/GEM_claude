#ifndef GEM_CORE_TYPES_HPP
#define GEM_CORE_TYPES_HPP

#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <array>
#include <string>

namespace gem {

// Basic types
using Real = double;
using Complex = std::complex<double>;

// Matrix types
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;

// Fixed-size matrices for SU(3) and SU(2)
using Matrix3c = Eigen::Matrix3cd;  // 3x3 complex for SU(3) color
using Matrix2c = Eigen::Matrix2cd;  // 2x2 complex for SU(2) spin

// Quark flavors
enum class QuarkFlavor {
    Up,
    Down,
    Strange,
    Charm,
    Bottom,
    Top
};

// Quark type (quark or antiquark)
enum class QuarkType {
    Quark,
    Antiquark
};

// Spin quantum number (S = 0, 1/2, 1, 3/2, ...)
struct Spin {
    int twice_S;  // 2*S to handle half-integers

    Spin() : twice_S(0) {}
    explicit Spin(int twice_s) : twice_S(twice_s) {}

    Real value() const { return static_cast<Real>(twice_S) / 2.0; }
    bool isHalfInteger() const { return twice_S % 2 != 0; }
};

// Quark properties
struct Quark {
    QuarkFlavor flavor;
    QuarkType type;
    Real mass;  // in MeV

    Quark() : flavor(QuarkFlavor::Up), type(QuarkType::Quark), mass(0.0) {}
    Quark(QuarkFlavor f, QuarkType t, Real m) : flavor(f), type(t), mass(m) {}
};

// Hadron types
enum class HadronType {
    Meson,      // q q̄
    Baryon,     // qqq
    Tetraquark, // qq q̄q̄
    Pentaquark  // qqqq q̄
};

// Color representation
enum class ColorRepresentation {
    Singlet,    // 1
    Triplet,    // 3
    Antitriplet,// 3̄
    Sextet,     // 6
    Antisextet, // 6̄
    Octet       // 8
};

// Quantum numbers for a hadron state
struct QuantumNumbers {
    int J;      // Total angular momentum (times 2 for half-integer)
    int P;      // Parity (+1 or -1)
    int C;      // C-parity (+1, -1, or 0 if not defined)
    int I;      // Isospin (times 2 for half-integer)
    int I3;     // Isospin z-component (times 2)
    int S_flavor; // Strangeness
    int charm;  // Charm quantum number
    int bottom; // Bottom quantum number

    QuantumNumbers() : J(0), P(1), C(0), I(0), I3(0), S_flavor(0), charm(0), bottom(0) {}
};

// Result of eigenvalue calculation
struct EigenResult {
    Vector eigenvalues;
    Matrix eigenvectors;
    bool converged;
    int iterations;

    EigenResult() : converged(false), iterations(0) {}
};

// Gaussian basis parameter set
struct BasisParameters {
    int n_basis;      // Number of basis functions
    Real r_min;       // Minimum radius in fm (gives largest nu)
    Real r_max;       // Maximum radius in fm (gives smallest nu)

    BasisParameters() : n_basis(15), r_min(0.1), r_max(3.0) {}
    BasisParameters(int n, Real rmin, Real rmax) : n_basis(n), r_min(rmin), r_max(rmax) {}

    // Compute nu_min and nu_max in fm^-2 from r_min/r_max
    // nu = 1/r^2, so larger r -> smaller nu
    Real nu_min() const { return 1.0 / (r_max * r_max); }  // fm^-2
    Real nu_max() const { return 1.0 / (r_min * r_min); }  // fm^-2
};

// Potential parameters
struct PotentialParameters {
    Real alpha_s_coul;   // Coulomb coupling constant
    Real alpha_s_hyp;    // Hyperfine coupling constant
    Real b;              // String tension (GeV^2)
    Real V0;             // Constant potential shift (GeV)

    // Hyperfine smearing: tau = mu^b_exp / a (mu in GeV)
    Real hyp_a;          // Smearing parameter a
    Real hyp_b_exp;      // Mass exponent b

    PotentialParameters()
        : alpha_s_coul(0.380175), alpha_s_hyp(1.39568),
          b(0.1653), V0(0.624075),
          hyp_a(1.6553), hyp_b_exp(0.2204) {}

    // Legacy getter for backward compatibility
    Real alpha_s() const { return alpha_s_coul; }
};

// Calculation configuration
struct CalculationConfig {
    BasisParameters basis;
    PotentialParameters potential;
    int max_iterations;
    Real convergence_threshold;
    bool verbose;

    CalculationConfig()
        : max_iterations(1000), convergence_threshold(1e-10), verbose(false) {}
};

} // namespace gem

#endif // GEM_CORE_TYPES_HPP
