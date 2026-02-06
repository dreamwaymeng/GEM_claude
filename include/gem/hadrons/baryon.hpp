#ifndef GEM_HADRONS_BARYON_HPP
#define GEM_HADRONS_BARYON_HPP

#include "../core/types.hpp"
#include "../core/config.hpp"
#include "../basis/gaussian_basis.hpp"
#include "../basis/jacobi_coordinates.hpp"
#include "../potential/potentials.hpp"
#include "../matrix/hamiltonian.hpp"
#include "../solver/eigenvalue_solver.hpp"

namespace gem {

/**
 * @brief Baryon (qqq) system
 *
 * A baryon consists of three quarks in a color singlet state.
 * The color wave function is antisymmetric: |1> = (1/√6) ε_{ijk} |i,j,k>
 * The color factor for each pair is <λ_i·λ_j>/4 = -2/3.
 *
 * Uses (ρ, λ) Jacobi coordinates:
 * ρ = r_1 - r_2
 * λ = (m_1 r_1 + m_2 r_2)/(m_1+m_2) - r_3
 *
 * For spin-1/2 baryons (proton, neutron, etc.): mixed symmetry spin
 * For spin-3/2 baryons (Δ, Ω, etc.): symmetric spin
 */
class Baryon {
public:
    /**
     * @brief Construct a baryon with given quark content
     * @param flavor1 Flavor of quark 1
     * @param flavor2 Flavor of quark 2
     * @param flavor3 Flavor of quark 3
     * @param twice_total_spin Total spin times 2 (1 for S=1/2, 3 for S=3/2)
     */
    Baryon(QuarkFlavor flavor1, QuarkFlavor flavor2, QuarkFlavor flavor3,
           int twice_total_spin = 1);

    /**
     * @brief Construct with custom masses
     */
    Baryon(Real m1, Real m2, Real m3, int twice_total_spin = 1);

    /**
     * @brief Set the basis parameters
     */
    void setBasis(const BasisParameters& params_rho, const BasisParameters& params_lambda);
    void setBasis(int n_basis, Real r_min, Real r_max);

    /**
     * @brief Set the potential parameters
     */
    void setPotential(const PotentialParameters& params);

    /**
     * @brief Set configuration
     */
    void setConfig(const Config& config);

    /**
     * @brief Calculate the baryon mass
     * @return The total mass in GeV
     */
    Real calculateMass();

    /**
     * @brief Get the ground state binding energy
     */
    Real bindingEnergy() const { return binding_energy_; }

    /**
     * @brief Get the constituent mass sum
     */
    Real constituentMass() const { return m1_ + m2_ + m3_; }

    /**
     * @brief Get the total mass
     */
    Real totalMass() const { return m1_ + m2_ + m3_ + binding_energy_; }

    /**
     * @brief Get color factor for baryon pairs
     */
    static constexpr Real colorFactor() { return constants::color_factor::BARYON_QQ; }

    /**
     * @brief Get the spin configuration
     */
    int twiceTotalSpin() const { return twice_total_spin_; }

    /**
     * @brief Check if calculation was successful
     */
    bool isValid() const { return valid_; }

    /**
     * @brief Get quark masses
     */
    Real mass1() const { return m1_; }
    Real mass2() const { return m2_; }
    Real mass3() const { return m3_; }

    /**
     * @brief Get reduced masses
     */
    Real mu_rho() const { return jacobi_.mu_rho(); }
    Real mu_lambda() const { return jacobi_.mu_lambda(); }

private:
    Real m1_, m2_, m3_;              // Quark masses
    int twice_total_spin_;           // 1 for S=1/2, 3 for S=3/2
    bool pair12_identical_;          // True if quarks 1,2 are same flavor (Σ-type)

    BasisParameters basis_params_rho_;
    BasisParameters basis_params_lambda_;
    PotentialParameters pot_params_;

    JacobiCoordinates3Body jacobi_;
    std::unique_ptr<GaussianBasis> basis_rho_;
    std::unique_ptr<GaussianBasis> basis_lambda_;

    Real binding_energy_;
    Vector wave_function_;
    Vector eigenvalues_;
    bool valid_;

    void initialize();
};

/**
 * @brief Factory for creating common baryons
 */
class BaryonFactory {
public:
    // Light baryons (S = 1/2)
    static Baryon proton() {
        return Baryon(QuarkFlavor::Up, QuarkFlavor::Up, QuarkFlavor::Down, 1);
    }
    static Baryon neutron() {
        // Put identical quarks (dd) in positions 1,2 for correct Σ-type spin structure
        return Baryon(QuarkFlavor::Down, QuarkFlavor::Down, QuarkFlavor::Up, 1);
    }
    static Baryon lambda() {
        return Baryon(QuarkFlavor::Up, QuarkFlavor::Down, QuarkFlavor::Strange, 1);
    }
    static Baryon sigma_plus() {
        return Baryon(QuarkFlavor::Up, QuarkFlavor::Up, QuarkFlavor::Strange, 1);
    }
    static Baryon xi() {
        // Put identical quarks (ss) in positions 1,2 for correct Σ-type spin structure
        return Baryon(QuarkFlavor::Strange, QuarkFlavor::Strange, QuarkFlavor::Up, 1);
    }

    // Decuplet baryons (S = 3/2)
    static Baryon delta_pp() {
        return Baryon(QuarkFlavor::Up, QuarkFlavor::Up, QuarkFlavor::Up, 3);
    }
    static Baryon sigma_star() {
        return Baryon(QuarkFlavor::Up, QuarkFlavor::Up, QuarkFlavor::Strange, 3);
    }
    static Baryon xi_star() {
        return Baryon(QuarkFlavor::Up, QuarkFlavor::Strange, QuarkFlavor::Strange, 3);
    }
    static Baryon omega() {
        return Baryon(QuarkFlavor::Strange, QuarkFlavor::Strange, QuarkFlavor::Strange, 3);
    }

    // Charmed baryons
    static Baryon lambda_c() {
        return Baryon(QuarkFlavor::Up, QuarkFlavor::Down, QuarkFlavor::Charm, 1);
    }
    static Baryon sigma_c() {
        return Baryon(QuarkFlavor::Up, QuarkFlavor::Up, QuarkFlavor::Charm, 1);
    }
    static Baryon omega_c() {
        return Baryon(QuarkFlavor::Strange, QuarkFlavor::Strange, QuarkFlavor::Charm, 1);
    }

    // Triple heavy baryons (S = 3/2)
    static Baryon omega_ccc() {
        return Baryon(QuarkFlavor::Charm, QuarkFlavor::Charm, QuarkFlavor::Charm, 3);
    }
    static Baryon omega_bbb() {
        return Baryon(QuarkFlavor::Bottom, QuarkFlavor::Bottom, QuarkFlavor::Bottom, 3);
    }
};

} // namespace gem

#endif // GEM_HADRONS_BARYON_HPP
