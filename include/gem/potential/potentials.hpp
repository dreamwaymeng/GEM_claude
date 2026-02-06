#ifndef GEM_POTENTIAL_POTENTIALS_HPP
#define GEM_POTENTIAL_POTENTIALS_HPP

#include "../core/types.hpp"
#include "../core/constants.hpp"
#include "../basis/gaussian_basis.hpp"

namespace gem {

/**
 * @brief Abstract base class for quark-quark potentials
 */
class Potential {
public:
    virtual ~Potential() = default;

    /**
     * @brief Evaluate the potential at distance r (in fm)
     */
    virtual Real evaluate(Real r) const = 0;

    /**
     * @brief Build the potential matrix in the Gaussian basis
     */
    virtual Matrix buildMatrix(const GaussianBasis& basis) const = 0;

    /**
     * @brief Get the potential strength/coupling
     */
    virtual Real strength() const = 0;
};

/**
 * @brief Coulomb-like potential: V(r) = -α_s / r × color_factor
 *
 * One-gluon exchange potential.
 */
class CoulombPotential : public Potential {
public:
    /**
     * @brief Construct Coulomb potential
     * @param alpha_s Strong coupling constant
     * @param color_factor Color factor (e.g., -4/3 for meson)
     */
    CoulombPotential(Real alpha_s, Real color_factor);

    Real evaluate(Real r) const override;
    Matrix buildMatrix(const GaussianBasis& basis) const override;
    Real strength() const override { return alpha_s_; }

    void setAlphaS(Real alpha) { alpha_s_ = alpha; }
    void setColorFactor(Real cf) { color_factor_ = cf; }

private:
    Real alpha_s_;
    Real color_factor_;
};

/**
 * @brief Linear confinement potential: V(r) = b × r × color_factor
 *
 * String tension between quarks.
 */
class ConfinementPotential : public Potential {
public:
    /**
     * @brief Construct confinement potential
     * @param b String tension (in GeV^2 or GeV/fm)
     * @param color_factor Color factor
     * @param V0 Constant shift (in MeV)
     */
    ConfinementPotential(Real b, Real color_factor, Real V0 = 0.0);

    Real evaluate(Real r) const override;
    Matrix buildMatrix(const GaussianBasis& basis) const override;
    Real strength() const override { return b_; }

    void setStringTension(Real b) { b_ = b; }
    void setColorFactor(Real cf) { color_factor_ = cf; }
    void setV0(Real v0) { V0_ = v0; }

private:
    Real b_;           // String tension
    Real color_factor_;
    Real V0_;          // Constant shift
};

/**
 * @brief Hyperfine (spin-spin) potential with Gaussian smearing
 *
 * V_hyp(r) = -(8π α_s^hyp / 3m_1 m_2) × (τ³/π^(3/2)) × exp(-τ²r²) × color × spin
 *
 * where τ = μ^b / a, with μ = 2m₁m₂/(m₁+m₂) being the reduced mass.
 * This is the contact term from one-gluon exchange, regularized with
 * a Gaussian smearing instead of a delta function.
 */
class HyperfinePotential : public Potential {
public:
    /**
     * @brief Construct hyperfine potential with mass-dependent tau
     * @param alpha_s_hyp Hyperfine coupling constant
     * @param m1 Mass of quark 1 (MeV)
     * @param m2 Mass of quark 2 (MeV)
     * @param hyp_a Smearing parameter a
     * @param hyp_b_exp Mass exponent b (tau = mu^b / a)
     * @param color_factor Color factor
     * @param spin_factor Spin factor <s_1·s_2>
     */
    HyperfinePotential(Real alpha_s_hyp, Real m1, Real m2,
                       Real hyp_a, Real hyp_b_exp,
                       Real color_factor, Real spin_factor);

    Real evaluate(Real r) const override;
    Matrix buildMatrix(const GaussianBasis& basis) const override;
    Real strength() const override { return alpha_s_; }

    void setMasses(Real m1, Real m2);
    void setSpinFactor(Real sf) { spin_factor_ = sf; }

    Real tau() const { return tau_; }  // Get computed tau value

private:
    Real alpha_s_;
    Real m1_, m2_;       // Quark masses (MeV)
    Real hyp_a_;         // Smearing parameter a
    Real hyp_b_exp_;     // Mass exponent b
    Real tau_;           // Computed smearing: tau = (mu_GeV)^b / a
    Real color_factor_;
    Real spin_factor_;

    void computeTau();   // Compute tau from masses
};

/**
 * @brief Cornell potential combining Coulomb and linear confinement
 *
 * V(r) = -α_s/r + b×r + V0
 */
class CornellPotential : public Potential {
public:
    CornellPotential(Real alpha_s, Real b, Real color_factor, Real V0 = 0.0);

    Real evaluate(Real r) const override;
    Matrix buildMatrix(const GaussianBasis& basis) const override;
    Real strength() const override { return coulomb_.strength(); }

    CoulombPotential& coulomb() { return coulomb_; }
    ConfinementPotential& confinement() { return confinement_; }

private:
    CoulombPotential coulomb_;
    ConfinementPotential confinement_;
};

/**
 * @brief Full quark-quark potential including all terms
 *
 * V_total = V_Coulomb + V_confinement + V_hyperfine
 */
class QuarkPotential {
public:
    QuarkPotential(const PotentialParameters& params, Real m1, Real m2,
                   Real color_factor, Real spin_factor);

    /**
     * @brief Evaluate total potential at distance r
     */
    Real evaluate(Real r) const;

    /**
     * @brief Build total potential matrix
     */
    Matrix buildMatrix(const GaussianBasis& basis) const;

    /**
     * @brief Access individual potential components
     */
    CoulombPotential& coulomb() { return coulomb_; }
    ConfinementPotential& confinement() { return confinement_; }
    HyperfinePotential& hyperfine() { return hyperfine_; }

    const CoulombPotential& coulomb() const { return coulomb_; }
    const ConfinementPotential& confinement() const { return confinement_; }
    const HyperfinePotential& hyperfine() const { return hyperfine_; }

private:
    CoulombPotential coulomb_;
    ConfinementPotential confinement_;
    HyperfinePotential hyperfine_;
};

} // namespace gem

#endif // GEM_POTENTIAL_POTENTIALS_HPP
