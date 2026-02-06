#include <iostream>
#include <iomanip>
#include <cmath>

#include "gem/core/types.hpp"
#include "gem/core/constants.hpp"
#include "gem/basis/gaussian_basis.hpp"
#include "gem/potential/potentials.hpp"
#include "gem/matrix/hamiltonian.hpp"
#include "gem/solver/eigenvalue_solver.hpp"
#include "gem/hadrons/meson.hpp"

using namespace gem;

void testCharmonium() {
    std::cout << "=== Charmonium (cc̄) Test with AL1 Parameters ===\n\n";

    Real m_c = constants::quark_mass::CHARM;  // 1.836 GeV
    std::cout << "Charm quark mass: " << m_c << " GeV (" << m_c * 1000 << " MeV)\n";
    std::cout << "Constituent mass sum: " << 2 * m_c << " GeV\n\n";

    // Basis parameters - using r_min/r_max in fm (matching Mathematica)
    BasisParameters basis_params;
    basis_params.n_basis = 20;
    basis_params.r_min = 0.0001;  // fm (gives largest nu)
    basis_params.r_max = 3.0;     // fm (gives smallest nu)

    GaussianBasis basis(basis_params);
    std::cout << "Basis size: " << basis.size() << "\n";
    std::cout << "r range: [" << basis_params.r_min << ", " << basis_params.r_max << "] fm\n";
    std::cout << "nu range: [" << basis.nu(0) << ", " << basis.nu(basis.size()-1) << "] fm^-2\n\n";

    // Reduced mass (in GeV)
    Real mu = m_c / 2.0;  // Equal mass case: μ = m/2
    std::cout << "Reduced mass μ = " << mu << " GeV\n\n";

    // AL1 model potential parameters
    PotentialParameters pot_params;  // Uses defaults from AL1 model
    std::cout << "AL1 Model Parameters:\n";
    std::cout << "  alpha_s (Coulomb) = " << pot_params.alpha_s_coul << "\n";
    std::cout << "  alpha_s (hyperfine) = " << pot_params.alpha_s_hyp << "\n";
    std::cout << "  b (string tension) = " << pot_params.b << " GeV^2\n";
    std::cout << "  V0 = " << pot_params.V0 << " GeV\n";
    std::cout << "  hyp_a = " << pot_params.hyp_a << "\n";
    std::cout << "  hyp_b_exp = " << pot_params.hyp_b_exp << "\n\n";

    Real color_factor = constants::color_factor::MESON_QQ_BAR;  // -4/3
    std::cout << "Color factor: " << color_factor << "\n\n";

    Real mass_s0 = 0, mass_s1 = 0;

    // Test both spin states
    for (int S = 0; S <= 1; ++S) {
        Real spin_factor = (S == 0) ? constants::spin_factor::SPIN_0
                                    : constants::spin_factor::SPIN_1;
        std::cout << "--- S = " << S << " ---\n";
        std::cout << "Spin factor <s1·s2> = " << spin_factor << "\n";

        // Create potential
        QuarkPotential potential(pot_params, m_c, m_c, color_factor, spin_factor);

        // Print computed tau
        std::cout << "Computed tau = " << potential.hyperfine().tau() << " GeV\n";

        // Build Hamiltonian
        TwoBodyHamiltonian hamiltonian(basis, mu, potential);
        hamiltonian.build(color_factor, spin_factor);

        // Solve with SVD-stabilized method
        EigenvalueSolver solver(EigenvalueSolver::Method::SvdStabilize);
        solver.setTolerance(1e-10);
        EigenResult result = solver.solve(hamiltonian);

        if (result.converged) {
            Real binding = result.eigenvalues(0);  // in GeV
            Real total_mass = 2 * m_c + binding;   // in GeV
            Real total_mass_mev = total_mass * constants::GEV_TO_MEV;  // in MeV
            std::cout << "Binding energy: " << binding << " GeV (" << binding * 1000 << " MeV)\n";
            std::cout << "Total mass: " << total_mass << " GeV (" << total_mass_mev << " MeV)\n";

            if (S == 0) mass_s0 = total_mass_mev;
            else mass_s1 = total_mass_mev;

            // Check condition number
            Matrix S_mat = basis.overlapMatrix();
            Eigen::JacobiSVD<Matrix> svd(S_mat);
            Real cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
            std::cout << "Effective basis size: " << result.eigenvalues.size() << "\n";
            std::cout << "Overlap condition number: " << cond << "\n";
        } else {
            std::cout << "Solver did not converge!\n";
        }
        std::cout << "\n";
    }

    std::cout << "=== Hyperfine Splitting ===\n";
    std::cout << "M(J/ψ) - M(η_c) = " << mass_s1 - mass_s0 << " MeV\n";
    std::cout << "Expected: ~113 MeV\n\n";

    // Compare with Meson class
    std::cout << "=== Using Meson class ===\n\n";

    Meson etac = MesonFactory::etac();
    Real mass_etac = etac.calculateMass() * constants::GEV_TO_MEV;  // Convert to MeV
    std::cout << "eta_c (S=0): " << mass_etac << " MeV\n";

    Meson jpsi = MesonFactory::jpsi();
    Real mass_jpsi = jpsi.calculateMass() * constants::GEV_TO_MEV;  // Convert to MeV
    std::cout << "J/psi (S=1): " << mass_jpsi << " MeV\n";

    std::cout << "\nHyperfine splitting: " << mass_jpsi - mass_etac << " MeV\n";
    std::cout << "Expected (experimental): 113 MeV\n";
    std::cout << "η_c experimental mass: 2984 MeV\n";
    std::cout << "J/ψ experimental mass: 3097 MeV\n";
}

int main() {
    testCharmonium();
    return 0;
}
