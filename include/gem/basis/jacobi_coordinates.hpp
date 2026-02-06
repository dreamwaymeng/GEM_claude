#ifndef GEM_BASIS_JACOBI_COORDINATES_HPP
#define GEM_BASIS_JACOBI_COORDINATES_HPP

#include "../core/types.hpp"
#include <vector>
#include <array>

namespace gem {

/**
 * @brief Jacobi coordinate system for N-body systems
 *
 * For an N-body system, there are N-1 relative Jacobi coordinates.
 *
 * 2-body (meson): r = r_1 - r_2
 *
 * 3-body (baryon):
 *   ρ (rho) = r_1 - r_2                           (relative coordinate)
 *   λ (lambda) = (m_1 r_1 + m_2 r_2)/(m_1+m_2) - r_3   (CM of 12 to 3)
 *
 * 4-body (tetraquark):
 *   Multiple possible arrangements depending on clustering
 *
 * 5-body (pentaquark):
 *   4 Jacobi coordinates
 */
class JacobiCoordinates {
public:
    /**
     * @brief Construct Jacobi coordinates for given particle masses
     * @param masses Vector of particle masses (in MeV)
     */
    explicit JacobiCoordinates(const std::vector<Real>& masses);

    /**
     * @brief Get the number of particles
     */
    int numParticles() const { return masses_.size(); }

    /**
     * @brief Get the number of Jacobi coordinates (N-1)
     */
    int numCoordinates() const { return masses_.size() - 1; }

    /**
     * @brief Get the reduced mass for Jacobi coordinate i
     * μ_i = m_eff,i × m_{remaining} / (m_eff,i + m_{remaining})
     */
    Real reducedMass(int i) const { return reduced_masses_[i]; }

    /**
     * @brief Get the mass of particle i
     */
    Real mass(int i) const { return masses_[i]; }

    /**
     * @brief Get the total mass of the system
     */
    Real totalMass() const { return total_mass_; }

    /**
     * @brief Get transformation coefficient T_{ij,k} for coordinate i
     * relating particle position r_j to Jacobi coordinate ξ_k
     */
    Real transformCoeff(int coord, int particle) const;

    /**
     * @brief Get the Jacobi coordinate matrix for two-particle interaction
     * Returns coefficients a_k, b for: r_i - r_j = Σ_k a_k ξ_k + b ξ_cm
     */
    Vector interactionCoeffs(int i, int j) const;

protected:
    std::vector<Real> masses_;
    std::vector<Real> reduced_masses_;
    Real total_mass_;
    Matrix transform_matrix_;  // Maps particle positions to Jacobi coordinates

    void computeReducedMasses();
    void computeTransformMatrix();
};

/**
 * @brief Specialized 2-body Jacobi coordinates (meson)
 */
class JacobiCoordinates2Body : public JacobiCoordinates {
public:
    JacobiCoordinates2Body(Real m1, Real m2);

    /**
     * @brief Get the single reduced mass μ = m1×m2/(m1+m2)
     */
    Real reducedMass() const { return reduced_masses_[0]; }
};

/**
 * @brief Specialized 3-body Jacobi coordinates (baryon)
 *
 * Standard (ρ, λ) coordinates:
 * ρ = r_1 - r_2
 * λ = (m_1 r_1 + m_2 r_2)/(m_1+m_2) - r_3
 */
class JacobiCoordinates3Body : public JacobiCoordinates {
public:
    JacobiCoordinates3Body(Real m1, Real m2, Real m3);

    /**
     * @brief Get reduced mass for ρ coordinate
     */
    Real mu_rho() const { return reduced_masses_[0]; }

    /**
     * @brief Get reduced mass for λ coordinate
     */
    Real mu_lambda() const { return reduced_masses_[1]; }

    /**
     * @brief Get coefficients for r_i - r_j in terms of (ρ, λ)
     * r_i - r_j = a_ρ × ρ + a_λ × λ
     */
    std::array<Real, 2> pairCoeffs(int i, int j) const;
};

/**
 * @brief Specialized 4-body Jacobi coordinates (tetraquark)
 *
 * Uses diquark-antidiquark clustering:
 * ρ_1 = r_1 - r_2        (quark-quark)
 * ρ_2 = r_3 - r_4        (antiquark-antiquark)
 * R = CM(12) - CM(34)    (diquark-antidiquark)
 */
class JacobiCoordinates4Body : public JacobiCoordinates {
public:
    JacobiCoordinates4Body(Real m1, Real m2, Real m3, Real m4);

    Real mu_rho1() const { return reduced_masses_[0]; }
    Real mu_rho2() const { return reduced_masses_[1]; }
    Real mu_R() const { return reduced_masses_[2]; }

    /**
     * @brief Get coefficients for particle pair in terms of Jacobi coordinates
     */
    std::array<Real, 3> pairCoeffs(int i, int j) const;
};

/**
 * @brief Specialized 5-body Jacobi coordinates (pentaquark)
 */
class JacobiCoordinates5Body : public JacobiCoordinates {
public:
    JacobiCoordinates5Body(Real m1, Real m2, Real m3, Real m4, Real m5);

    /**
     * @brief Get coefficients for particle pair in terms of Jacobi coordinates
     */
    std::array<Real, 4> pairCoeffs(int i, int j) const;
};

/**
 * @brief Raynal-Revai transformation coefficients
 *
 * Used to transform between different Jacobi coordinate sets.
 * Particularly important for antisymmetrization in baryon calculations.
 */
class RaynalRevaiCoeffs {
public:
    /**
     * @brief Compute RR coefficient for 3-body system
     * Transforms from (12)3 to (23)1 arrangement
     */
    static Real coefficient_3body(int n_rho, int l_rho, int n_lambda, int l_lambda,
                                  int np_rho, int lp_rho, int np_lambda, int lp_lambda,
                                  Real mass1, Real mass2, Real mass3);

    /**
     * @brief Compute the transformation angle for 3-body RR
     */
    static Real transformationAngle(Real m1, Real m2, Real m3);

private:
    static Real moshinskyBracket(int n1, int l1, int n2, int l2,
                                  int N, int L, int nr, int lr, int l_total);
};

} // namespace gem

#endif // GEM_BASIS_JACOBI_COORDINATES_HPP
