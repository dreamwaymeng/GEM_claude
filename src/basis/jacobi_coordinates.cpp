#include "gem/basis/jacobi_coordinates.hpp"
#include <cmath>
#include <stdexcept>
#include <numeric>

namespace gem {

// JacobiCoordinates base class implementation

JacobiCoordinates::JacobiCoordinates(const std::vector<Real>& masses)
    : masses_(masses), reduced_masses_(masses.size() - 1), total_mass_(0.0) {
    if (masses.size() < 2) {
        throw std::invalid_argument("Need at least 2 particles for Jacobi coordinates");
    }
    total_mass_ = std::accumulate(masses.begin(), masses.end(), 0.0);
    computeReducedMasses();
    computeTransformMatrix();
}

void JacobiCoordinates::computeReducedMasses() {
    int n = masses_.size();

    // For standard sequential Jacobi coordinates:
    // ξ_1 = r_1 - r_2, μ_1 = m_1 m_2 / (m_1 + m_2)
    // ξ_2 = CM(1,2) - r_3, μ_2 = (m_1+m_2) m_3 / (m_1+m_2+m_3)
    // etc.

    Real cumulative_mass = masses_[0];
    for (int i = 0; i < n - 1; ++i) {
        Real m_next = masses_[i + 1];
        if (i == 0) {
            // First coordinate: between particles 1 and 2
            reduced_masses_[i] = masses_[0] * masses_[1] / (masses_[0] + masses_[1]);
        } else {
            // Subsequent coordinates: CM of previous cluster to next particle
            reduced_masses_[i] = cumulative_mass * m_next / (cumulative_mass + m_next);
        }
        cumulative_mass += m_next;
    }
}

void JacobiCoordinates::computeTransformMatrix() {
    int n = masses_.size();
    transform_matrix_ = Matrix::Zero(n - 1, n);

    // Row i gives coefficients for Jacobi coordinate ξ_i in terms of particle positions
    // ξ_i = Σ_j T_{ij} r_j

    // First Jacobi coordinate: ξ_0 = r_0 - r_1
    transform_matrix_(0, 0) = 1.0;
    transform_matrix_(0, 1) = -1.0;

    // Subsequent coordinates
    Real cumulative_mass = masses_[0] + masses_[1];
    for (int i = 1; i < n - 1; ++i) {
        // ξ_i = CM(0..i) - r_{i+1}
        // CM(0..i) = Σ_{j=0}^{i} m_j r_j / Σ_{j=0}^{i} m_j
        for (int j = 0; j <= i; ++j) {
            transform_matrix_(i, j) = masses_[j] / cumulative_mass;
        }
        transform_matrix_(i, i + 1) = -1.0;
        cumulative_mass += masses_[i + 1];
    }
}

Real JacobiCoordinates::transformCoeff(int coord, int particle) const {
    if (coord < 0 || coord >= numCoordinates() ||
        particle < 0 || particle >= numParticles()) {
        throw std::out_of_range("Invalid coordinate or particle index");
    }
    return transform_matrix_(coord, particle);
}

Vector JacobiCoordinates::interactionCoeffs(int i, int j) const {
    // r_i - r_j in terms of Jacobi coordinates
    // Need to invert the transformation or compute directly
    int n_coords = numCoordinates();
    Vector coeffs = Vector::Zero(n_coords);

    // For each Jacobi coordinate, the contribution to r_i - r_j
    // is given by the inverse transformation
    // This is a simplified version for sequential clustering

    if (i == j) return coeffs;

    // For a 2-body system: r_0 - r_1 = ξ_0
    if (numParticles() == 2) {
        coeffs(0) = (i == 0) ? 1.0 : -1.0;
        return coeffs;
    }

    // For 3-body and beyond, need more careful treatment
    // This is handled in specialized classes
    return coeffs;
}

// JacobiCoordinates2Body implementation

JacobiCoordinates2Body::JacobiCoordinates2Body(Real m1, Real m2)
    : JacobiCoordinates({m1, m2}) {}

// JacobiCoordinates3Body implementation

JacobiCoordinates3Body::JacobiCoordinates3Body(Real m1, Real m2, Real m3)
    : JacobiCoordinates({m1, m2, m3}) {}

std::array<Real, 2> JacobiCoordinates3Body::pairCoeffs(int i, int j) const {
    // Coefficients (a_ρ, a_λ) such that r_i - r_j = a_ρ × ρ + a_λ × λ
    // where ρ = r_1 - r_2, λ = CM(12) - r_3

    std::array<Real, 2> coeffs = {0.0, 0.0};
    Real m1 = masses_[0], m2 = masses_[1], m3 = masses_[2];
    Real m12 = m1 + m2;

    if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
        // r_1 - r_2 = ρ
        coeffs[0] = (i == 0) ? 1.0 : -1.0;
        coeffs[1] = 0.0;
    } else if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
        // r_1 - r_3 = r_1 - CM(12) + CM(12) - r_3
        //           = r_1 - (m1 r_1 + m2 r_2)/m12 + λ
        //           = m2/m12 × (r_1 - r_2) + λ
        //           = m2/m12 × ρ + λ
        coeffs[0] = (i == 0) ? m2 / m12 : -m2 / m12;
        coeffs[1] = (i == 0) ? 1.0 : -1.0;
    } else if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
        // r_2 - r_3 = r_2 - CM(12) + CM(12) - r_3
        //           = r_2 - (m1 r_1 + m2 r_2)/m12 + λ
        //           = -m1/m12 × (r_1 - r_2) + λ
        //           = -m1/m12 × ρ + λ
        coeffs[0] = (i == 1) ? -m1 / m12 : m1 / m12;
        coeffs[1] = (i == 1) ? 1.0 : -1.0;
    }

    return coeffs;
}

// JacobiCoordinates4Body implementation

JacobiCoordinates4Body::JacobiCoordinates4Body(Real m1, Real m2, Real m3, Real m4)
    : JacobiCoordinates({m1, m2, m3, m4}) {
    // Override reduced masses for diquark-antidiquark clustering
    // ρ_1 = r_1 - r_2, μ_1 = m1 m2 / (m1 + m2)
    // ρ_2 = r_3 - r_4, μ_2 = m3 m4 / (m3 + m4)
    // R = CM(12) - CM(34), μ_R = (m1+m2)(m3+m4) / (m1+m2+m3+m4)

    reduced_masses_[0] = m1 * m2 / (m1 + m2);
    reduced_masses_[1] = m3 * m4 / (m3 + m4);
    reduced_masses_[2] = (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4);

    // Update transform matrix for diquark clustering
    transform_matrix_ = Matrix::Zero(3, 4);
    Real m12 = m1 + m2;
    Real m34 = m3 + m4;

    // ρ_1 = r_1 - r_2
    transform_matrix_(0, 0) = 1.0;
    transform_matrix_(0, 1) = -1.0;

    // ρ_2 = r_3 - r_4
    transform_matrix_(1, 2) = 1.0;
    transform_matrix_(1, 3) = -1.0;

    // R = CM(12) - CM(34)
    transform_matrix_(2, 0) = m1 / m12;
    transform_matrix_(2, 1) = m2 / m12;
    transform_matrix_(2, 2) = -m3 / m34;
    transform_matrix_(2, 3) = -m4 / m34;
}

std::array<Real, 3> JacobiCoordinates4Body::pairCoeffs(int i, int j) const {
    std::array<Real, 3> coeffs = {0.0, 0.0, 0.0};
    Real m1 = masses_[0], m2 = masses_[1], m3 = masses_[2], m4 = masses_[3];
    Real m12 = m1 + m2;
    Real m34 = m3 + m4;

    // Determine which pair and compute coefficients
    int pair = (i < j) ? (i * 10 + j) : (j * 10 + i);
    Real sign = (i < j) ? 1.0 : -1.0;

    switch (pair) {
        case 1:   // (0,1): r_1 - r_2 = ρ_1
            coeffs[0] = sign;
            break;
        case 23:  // (2,3): r_3 - r_4 = ρ_2
            coeffs[1] = sign;
            break;
        case 2:   // (0,2): r_1 - r_3 = m2/m12 ρ_1 - m4/m34 ρ_2 + R
            coeffs[0] = sign * m2 / m12;
            coeffs[1] = -sign * m4 / m34;
            coeffs[2] = sign;
            break;
        case 3:   // (0,3): r_1 - r_4 = m2/m12 ρ_1 + m3/m34 ρ_2 + R
            coeffs[0] = sign * m2 / m12;
            coeffs[1] = sign * m3 / m34;
            coeffs[2] = sign;
            break;
        case 12:  // (1,2): r_2 - r_3 = -m1/m12 ρ_1 - m4/m34 ρ_2 + R
            coeffs[0] = -sign * m1 / m12;
            coeffs[1] = -sign * m4 / m34;
            coeffs[2] = sign;
            break;
        case 13:  // (1,3): r_2 - r_4 = -m1/m12 ρ_1 + m3/m34 ρ_2 + R
            coeffs[0] = -sign * m1 / m12;
            coeffs[1] = sign * m3 / m34;
            coeffs[2] = sign;
            break;
    }

    return coeffs;
}

// JacobiCoordinates5Body implementation

JacobiCoordinates5Body::JacobiCoordinates5Body(Real m1, Real m2, Real m3, Real m4, Real m5)
    : JacobiCoordinates({m1, m2, m3, m4, m5}) {}

std::array<Real, 4> JacobiCoordinates5Body::pairCoeffs(int i, int j) const {
    // Simplified implementation - full version would need specific clustering choice
    std::array<Real, 4> coeffs = {0.0, 0.0, 0.0, 0.0};

    if (i == j || i < 0 || j < 0 || i > 4 || j > 4) {
        return coeffs;
    }

    // Use transform matrix to compute coefficients
    // r_i - r_j = Σ_k c_k ξ_k
    // This requires inverting the Jacobi transformation

    // For now, return placeholder - full implementation would compute
    // coefficients based on specific clustering scheme

    return coeffs;
}

// RaynalRevaiCoeffs implementation

Real RaynalRevaiCoeffs::transformationAngle(Real m1, Real m2, Real m3) {
    // The transformation angle θ between (12)3 and (23)1 clustering
    // tan(θ) = √(m1 m3 / (m2 (m1 + m2 + m3)))
    Real total = m1 + m2 + m3;
    return std::atan(std::sqrt(m1 * m3 / (m2 * total)));
}

Real RaynalRevaiCoeffs::coefficient_3body(int n_rho, int l_rho, int n_lambda, int l_lambda,
                                          int np_rho, int lp_rho, int np_lambda, int lp_lambda,
                                          Real mass1, Real mass2, Real mass3) {
    // This is a simplified placeholder
    // Full implementation would use Moshinsky brackets and transformation angle

    // For l=0 case only (s-wave)
    if (l_rho != 0 || l_lambda != 0 || lp_rho != 0 || lp_lambda != 0) {
        return 0.0;  // Higher partial waves not implemented
    }

    Real theta = transformationAngle(mass1, mass2, mass3);
    Real c = std::cos(theta);
    Real s = std::sin(theta);

    // For n=0 case (ground state)
    if (n_rho == 0 && n_lambda == 0 && np_rho == 0 && np_lambda == 0) {
        return 1.0;  // Identity transformation for ground state
    }

    // For general case, need Talmi-Moshinsky brackets
    // This is a complex calculation involving Clebsch-Gordan coefficients
    // and generalized Laguerre polynomials

    return 0.0;  // Placeholder for full implementation
}

Real RaynalRevaiCoeffs::moshinskyBracket(int n1, int l1, int n2, int l2,
                                          int N, int L, int nr, int lr, int l_total) {
    // Moshinsky bracket <n1 l1, n2 l2; l_total | N L, nr lr; l_total>
    // This is a complex calculation - placeholder for full implementation
    return 0.0;
}

} // namespace gem
