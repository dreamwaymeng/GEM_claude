#ifndef GEM_HADRONS_MESON_HPP
#define GEM_HADRONS_MESON_HPP

#include "../core/types.hpp"
#include "../core/config.hpp"
#include "../basis/gaussian_basis.hpp"
#include "../basis/jacobi_coordinates.hpp"
#include "../quantum/wave_function.hpp"
#include "../potential/potentials.hpp"
#include "../matrix/hamiltonian.hpp"
#include "../solver/eigenvalue_solver.hpp"

namespace gem {

/**
 * @brief Meson (qq̄) system
 *
 * A meson consists of a quark and an antiquark.
 * The color wave function is the singlet: |1> = (1/√3)(|rr̄> + |gḡ> + |bb̄>)
 * The color factor is <λ·λ>/4 = -4/3.
 *
 * For pseudoscalar mesons (π, K, etc.): J^PC = 0^-+, S = 0
 * For vector mesons (ρ, ω, etc.): J^PC = 1^--, S = 1
 */
class Meson {
public:
    /**
     * @brief Construct a meson with given quark content
     * @param flavor1 Flavor of the quark
     * @param flavor2 Flavor of the antiquark
     * @param total_spin Total spin (0 for pseudoscalar, 1 for vector)
     */
    Meson(QuarkFlavor flavor1, QuarkFlavor flavor2, int total_spin = 0);

    /**
     * @brief Construct with custom masses
     */
    Meson(Real m1, Real m2, int total_spin = 0);

    /**
     * @brief Set the basis parameters
     */
    void setBasis(const BasisParameters& params);
    void setBasis(int n_basis, Real r_min, Real r_max);

    /**
     * @brief Set the potential parameters
     */
    void setPotential(const PotentialParameters& params);

    /**
     * @brief Set the configuration
     */
    void setConfig(const Config& config);

    /**
     * @brief Calculate the meson mass
     * @return The total mass (constituent masses + binding energy) in GeV
     */
    Real calculateMass();

    /**
     * @brief Get the ground state binding energy
     */
    Real bindingEnergy() const { return binding_energy_; }

    /**
     * @brief Get the constituent mass sum
     */
    Real constituentMass() const { return m1_ + m2_; }

    /**
     * @brief Get the total mass (constituent + binding)
     */
    Real totalMass() const { return m1_ + m2_ + binding_energy_; }

    /**
     * @brief Get the RMS radius of the meson
     */
    Real rmsRadius() const { return rms_radius_; }

    /**
     * @brief Get the reduced mass
     */
    Real reducedMass() const { return jacobi_.reducedMass(); }

    /**
     * @brief Get color factor
     */
    static constexpr Real colorFactor() { return constants::color_factor::MESON_QQ_BAR; }

    /**
     * @brief Get spin factor for the meson
     */
    Real spinFactor() const;

    /**
     * @brief Get the wave function coefficients
     */
    const Vector& waveFunction() const { return wave_function_; }

    /**
     * @brief Get all eigenvalues (excited states)
     */
    const Vector& eigenvalues() const { return eigenvalues_; }

    /**
     * @brief Get the n-th excited state energy
     */
    Real excitedStateEnergy(int n) const;

    /**
     * @brief Check if calculation was successful
     */
    bool isValid() const { return valid_; }

    /**
     * @brief Get quark masses
     */
    Real mass1() const { return m1_; }
    Real mass2() const { return m2_; }

    /**
     * @brief Get total spin
     */
    int totalSpin() const { return total_spin_; }

private:
    Real m1_, m2_;                    // Quark masses
    int total_spin_;                  // 0 (pseudoscalar) or 1 (vector)

    BasisParameters basis_params_;
    PotentialParameters pot_params_;

    JacobiCoordinates2Body jacobi_;
    std::unique_ptr<GaussianBasis> basis_;
    std::unique_ptr<QuarkPotential> potential_;

    Real binding_energy_;
    Real rms_radius_;
    Vector wave_function_;
    Vector eigenvalues_;
    bool valid_;

    void initialize();
};

/**
 * @brief Factory for creating common mesons
 */
class MesonFactory {
public:
    // Light mesons
    static Meson pion() { return Meson(QuarkFlavor::Up, QuarkFlavor::Down, 0); }
    static Meson rho() { return Meson(QuarkFlavor::Up, QuarkFlavor::Down, 1); }
    static Meson kaon() { return Meson(QuarkFlavor::Up, QuarkFlavor::Strange, 0); }
    static Meson kstar() { return Meson(QuarkFlavor::Up, QuarkFlavor::Strange, 1); }

    // Heavy-light mesons
    static Meson dmeson() { return Meson(QuarkFlavor::Charm, QuarkFlavor::Down, 0); }
    static Meson dstar() { return Meson(QuarkFlavor::Charm, QuarkFlavor::Down, 1); }
    static Meson bmeson() { return Meson(QuarkFlavor::Bottom, QuarkFlavor::Up, 0); }

    // Heavy quarkonia
    static Meson jpsi() { return Meson(QuarkFlavor::Charm, QuarkFlavor::Charm, 1); }
    static Meson etac() { return Meson(QuarkFlavor::Charm, QuarkFlavor::Charm, 0); }
    static Meson upsilon() { return Meson(QuarkFlavor::Bottom, QuarkFlavor::Bottom, 1); }
    static Meson etab() { return Meson(QuarkFlavor::Bottom, QuarkFlavor::Bottom, 0); }
};

} // namespace gem

#endif // GEM_HADRONS_MESON_HPP
