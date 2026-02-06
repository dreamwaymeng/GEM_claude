#ifndef GEM_MATRIX_HAMILTONIAN_HPP
#define GEM_MATRIX_HAMILTONIAN_HPP

#include "../core/types.hpp"
#include "../basis/gaussian_basis.hpp"
#include "../potential/potentials.hpp"

namespace gem {

/**
 * @brief Hamiltonian matrix builder for GEM calculations
 *
 * H = T + V = Σ_i p_i²/(2m_i) + Σ_{i<j} V_{ij}
 *
 * In Jacobi coordinates for 2-body:
 * H = T_rel + V(r) where T_rel = p_rel²/(2μ)
 *
 * The matrix elements are computed in the Gaussian basis.
 */
class Hamiltonian {
public:
    /**
     * @brief Get the Hamiltonian matrix
     */
    const Matrix& matrix() const { return H_; }

    /**
     * @brief Get the overlap matrix
     */
    const Matrix& overlapMatrix() const { return S_; }

    /**
     * @brief Get the kinetic energy matrix
     */
    const Matrix& kineticMatrix() const { return T_; }

    /**
     * @brief Get the potential energy matrix
     */
    const Matrix& potentialMatrix() const { return V_; }

    /**
     * @brief Get the basis size
     */
    int size() const { return H_.rows(); }

protected:
    Matrix H_;  // Full Hamiltonian
    Matrix S_;  // Overlap matrix
    Matrix T_;  // Kinetic energy
    Matrix V_;  // Potential energy

    Hamiltonian() = default;
};

/**
 * @brief Two-body Hamiltonian (for mesons)
 *
 * H = p²/(2μ) + V(r)
 * where μ is the reduced mass and V(r) includes all potential terms.
 */
class TwoBodyHamiltonian : public Hamiltonian {
public:
    /**
     * @brief Construct two-body Hamiltonian
     * @param basis Gaussian basis
     * @param reduced_mass Reduced mass μ = m1 m2 / (m1 + m2) in MeV
     * @param potential Quark-quark potential
     */
    TwoBodyHamiltonian(const GaussianBasis& basis, Real reduced_mass,
                       const QuarkPotential& potential);

    /**
     * @brief Build Hamiltonian with specified color and spin factors
     */
    void build(Real color_factor, Real spin_factor);

    /**
     * @brief Rebuild with new potential parameters
     */
    void rebuild(const PotentialParameters& params);

    const GaussianBasis& basis() const { return basis_; }
    Real reducedMass() const { return reduced_mass_; }

private:
    const GaussianBasis& basis_;
    Real reduced_mass_;
    QuarkPotential potential_;
};

/**
 * @brief Three-body Hamiltonian (for baryons)
 *
 * Uses (ρ, λ) Jacobi coordinates.
 * H = T_ρ + T_λ + Σ_{i<j} V_{ij}
 */
class ThreeBodyHamiltonian : public Hamiltonian {
public:
    /**
     * @brief Construct three-body Hamiltonian
     * @param basis_rho Basis for ρ coordinate
     * @param basis_lambda Basis for λ coordinate
     * @param mu_rho Reduced mass for ρ
     * @param mu_lambda Reduced mass for λ
     * @param masses Quark masses {m1, m2, m3}
     * @param params Potential parameters
     */
    ThreeBodyHamiltonian(const GaussianBasis& basis_rho, const GaussianBasis& basis_lambda,
                         Real mu_rho, Real mu_lambda,
                         const std::array<Real, 3>& masses,
                         const PotentialParameters& params);

    /**
     * @brief Build Hamiltonian with specified total spin
     * @param twice_S Total spin times 2 (1 for S=1/2, 3 for S=3/2)
     */
    void build(int twice_S);

private:
    const GaussianBasis& basis_rho_;
    const GaussianBasis& basis_lambda_;
    Real mu_rho_, mu_lambda_;
    std::array<Real, 3> masses_;
    PotentialParameters params_;

    void buildKineticMatrix();
    void buildPotentialMatrix(int twice_S);
};

/**
 * @brief Helper class to compute matrix elements in GEM
 */
class MatrixElements {
public:
    /**
     * @brief Compute kinetic energy matrix element
     * <φ_n|T|φ_m> = (ℏ²/2μ) × 3νnνm/(νn+νm) × <φ_n|φ_m>
     */
    static Real kinetic(const GaussianBasis& basis, int n, int m, Real mu);

    /**
     * @brief Compute Coulomb potential matrix element
     * <φ_n|1/r|φ_m> = Nn Nm × 2π/(νn+νm)
     */
    static Real coulomb(const GaussianBasis& basis, int n, int m);

    /**
     * @brief Compute linear potential matrix element
     * <φ_n|r|φ_m> = Nn Nm × 2π^(3/2)/(νn+νm)^(5/2)
     */
    static Real linear(const GaussianBasis& basis, int n, int m);

    /**
     * @brief Compute Gaussian potential matrix element
     * <φ_n|exp(-αr²)|φ_m> = Nn Nm × (π/(νn+νm+α))^(3/2)
     */
    static Real gaussian(const GaussianBasis& basis, int n, int m, Real alpha);

    /**
     * @brief Compute two-body potential matrix element for three-body system
     * The pair distance r_ij is expressed in terms of (ρ, λ) coordinates
     * r_ij² = a² ρ² + b² λ² + 2ab ρ·λ cos(θ)
     *
     * For s-wave: integrate over angles to get effective radial matrix element
     */
    static Real twoBodyInThreeBody(const GaussianBasis& basis_rho,
                                   const GaussianBasis& basis_lambda,
                                   int n_rho, int n_lambda, int m_rho, int m_lambda,
                                   Real a, Real b,  // Coefficients for ρ and λ
                                   const std::function<Real(Real)>& potential);
};

} // namespace gem

#endif // GEM_MATRIX_HAMILTONIAN_HPP
