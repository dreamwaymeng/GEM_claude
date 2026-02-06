## Gaussian Expansion Method (GEM)

The Gaussian Expansion Method (GEM) is a variational approach for solving few-body quantum mechanical problems. The wave function is expanded using Gaussian basis functions:

$$
\psi_{lm}(\vec{r}) = \sum_{n=1} c_{nl} \phi_{nlm}^{G}(\vec{r})
$$

$$
\phi_{nlm}^{G}(\vec{r}) = \phi_{nl}^{G}(r) Y_{lm}(\hat{r})
$$

$$
\phi_{nl}^{G}(r) = N_{nl} r^{l} e^{-\nu_{n} r^{2}}
$$

**Key features:**
- In this work, only the $l=0$ (s-wave) basis is used.
- The basis with $\nu_n$ arranged in a geometric series has been shown to be highly efficient.
- The main advantage of GEM is that all matrix elements can be calculated analytically.

## Analytical Matrix Elements for $l=0$ Gaussian Basis

For the s-wave ($l=0$) normalized Gaussian basis function:

$$
\phi_n(r) = N_n e^{-\nu_n r^2}, \quad N_n = \left( \frac{2\nu_n}{\pi} \right)^{3/4}
$$

All spatial matrix elements can be computed analytically.

### Overlap Integral

$$
\langle \phi_n | \phi_m \rangle = N_n N_m \left( \frac{\pi}{\nu_n + \nu_m} \right)^{3/2}
$$

### Kinetic Energy Integral

The Laplacian of the Gaussian basis function is:

$$
\nabla^2 e^{-\nu r^2} = (4\nu^2 r^2 - 6\nu) e^{-\nu r^2}
$$

The Laplacian matrix element is:

$$
\langle \phi_n | \nabla^2 | \phi_m \rangle = -\frac{6\nu_n \nu_m}{\nu_n + \nu_m} \langle \phi_n | \phi_m \rangle
$$

The kinetic energy operator is $T = -\frac{\hbar^2}{2\mu} \nabla^2$, so:

$$
\langle \phi_n | T | \phi_m \rangle = \frac{\hbar^2}{2\mu} \cdot \frac{6\nu_n \nu_m}{\nu_n + \nu_m} \langle \phi_n | \phi_m \rangle
$$

In units where $r$ is in fm and $\nu$ is in fm$^{-2}$:

$$
\langle \phi_n | T | \phi_m \rangle = \frac{(\hbar c)^2}{2\mu} \cdot \frac{6\nu_n \nu_m}{\nu_n + \nu_m} \langle \phi_n | \phi_m \rangle
$$

where $(\hbar c)^2 = 0.0389$ GeV$^2$·fm$^2$ and $\mu$ is in GeV.

### Coulomb Integral

$$
\langle \phi_n | \frac{1}{r} | \phi_m \rangle = N_n N_m \frac{2\pi}{\nu_n + \nu_m}
$$

### Linear Integral

Using the Gaussian integral $\int_0^\infty r^3 e^{-\alpha r^2} dr = \frac{1}{2\alpha^2}$:

$$
\langle \phi_n | r | \phi_m \rangle = N_n N_m \frac{2\pi}{(\nu_n + \nu_m)^2}
$$

### Gaussian Potential Integral

For the smeared delta function in the hyperfine potential:

$$
\langle \phi_n | e^{-\alpha r^2} | \phi_m \rangle = N_n N_m \left( \frac{\pi}{\nu_n + \nu_m + \alpha} \right)^{3/2}
$$

### Summary Table

| Matrix Element | Formula |
|----------------|---------|
| $\langle \phi_n \| \phi_m \rangle$ | $N_n N_m \left( \frac{\pi}{\nu_n + \nu_m} \right)^{3/2}$ |
| $\langle \phi_n \| T \| \phi_m \rangle$ | $\frac{(\hbar c)^2}{2\mu} \cdot \frac{6\nu_n \nu_m}{\nu_n + \nu_m} \langle \phi_n \| \phi_m \rangle$ |
| $\langle \phi_n \| 1/r \| \phi_m \rangle$ | $N_n N_m \frac{2\pi}{\nu_n + \nu_m}$ |
| $\langle \phi_n \| r \| \phi_m \rangle$ | $N_n N_m \frac{2\pi}{(\nu_n + \nu_m)^2}$ |
| $\langle \phi_n \| e^{-\alpha r^2} \| \phi_m \rangle$ | $N_n N_m \left( \frac{\pi}{\nu_n + \nu_m + \alpha} \right)^{3/2}$ |

## Quark Potential Model

The minimal quark model includes three interaction terms: Coulomb interaction, hyperfine interaction, and linear confinement.

### Coulomb Potential (One-Gluon Exchange)

$$
V_{\text{Coulomb}}(r) = \frac{\alpha_s^{\text{coul}}}{r} \frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2}
$$

For mesons: $\frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2} = -\frac{4}{3}$, giving an attractive potential.

### Hyperfine Potential (Spin-Spin Interaction)

$$
V_{\text{hyp}}(r) = -\frac{8\pi\alpha_s^{\text{hyp}}}{3m_i m_j} \frac{\tau^3}{\pi^{3/2}} e^{-\tau^2 r^2} \frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2} \, \mathbf{s}_i \cdot \mathbf{s}_j
$$

where $\tau$ depends on the quark masses:

$$
\tau = \frac{(2\mu_{ij})^{b_{exp}}}{a}, \quad \mu_{ij} = \frac{m_i m_j}{m_i + m_j}
$$

### Confinement Potential

$$
V_{\text{conf}}(r) = \left( -\frac{3b}{4} r + V_c \right) \frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2}
$$

### Physical Interpretation

- The Coulomb and confinement interactions depend only on the color-electric interaction, where $\boldsymbol{\lambda}_i$ are the Gell-Mann matrices in color space.
- The hyperfine interaction is a color-magnetic interaction that depends on the quark spins $\mathbf{s}_i$.
- **Important:** The coupling constants $\alpha_s^{\text{coul}}$ and $\alpha_s^{\text{hyp}}$ are different, reflecting different effective couplings for color-electric and color-magnetic interactions.
- Due to color confinement, hadrons composed of quarks must form color singlet states.
- Conventional hadrons include mesons ($q\bar{q}$) and baryons ($qqq$). Multiquark states such as tetraquarks and pentaquarks are also possible.

## Model Parameters

All parameters are in units of GeV (or appropriate powers of GeV).

### Quark Masses

| Quark | Mass (GeV) |
|-------|------------|
| u, d  | 0.315      |
| s     | 0.577      |
| c     | 1.836      |
| b     | 5.227      |

### Potential Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| $\alpha_s^{\text{coul}}$ | 0.380175 | Coulomb coupling constant |
| $\alpha_s^{\text{hyp}}$  | 1.39568  | Hyperfine coupling constant |
| $a$       | 1.6553   | Hyperfine smearing parameter |
| $b_{exp}$ | 0.2204   | Hyperfine mass exponent |
| $b$       | 0.1653   | String tension (GeV²) |
| $V_c$     | 0.624075 | Confinement constant (GeV) |

### Hyperfine Smearing Parameter

The smearing parameter $\tau$ for the hyperfine interaction depends on the reduced mass of the quark pair:

$$
\tau_{ij} = \frac{(2\mu_{ij})^{b_{exp}}}{a}
$$

where $\mu_{ij} = m_i m_j/(m_i + m_j)$ is the reduced mass.

**Examples:**

- For $c\bar{c}$: $\mu = m_c/2 = 0.918$ GeV, $\tau = (2 \times 0.918)^{0.2204}/1.6553 = 1.836^{0.2204}/1.6553 \approx 0.69$ GeV
- For $u\bar{d}$: $\mu = m_u/2 = 0.1575$ GeV, $\tau = (2 \times 0.1575)^{0.2204}/1.6553 = 0.315^{0.2204}/1.6553 \approx 0.47$ GeV

## Calculating Hadron Masses with the Quark Model

To obtain hadron masses, one solves the non-relativistic Schrödinger equation for the (anti)quark system:

1. **Construct the wave function:** Write complete discrete wave functions for hadrons in terms of quark degrees of freedom, subject to the constraints of color singlet, and fixed total spin and isospin.

2. **Calculate discrete matrix elements:** Using these wave functions, compute the matrix elements of $\boldsymbol{\lambda}_i \cdot \boldsymbol{\lambda}_j$ and $\boldsymbol{\lambda}_i \cdot \boldsymbol{\lambda}_j \, \mathbf{s}_i \cdot \mathbf{s}_j$.

3. **Expand spatial wave functions:** Expand the spatial wave function in Gaussian bases defined in different Jacobi coordinates. With GEM, all spatial matrix elements can be calculated analytically.

4. **Combine matrix elements:** The total matrix element is obtained by combining the discrete (color, spin, flavor) and spatial contributions.

5. **Solve the eigenvalue problem:** The masses of multiquark states are obtained by solving the generalized eigenvalue problem.

## Color and Spin Factors

### Meson ($q\bar{q}$)

Color singlet: $\frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2} = -\frac{4}{3}$

Spin factors:
- $S=0$ (pseudoscalar): $\langle \mathbf{s}_1 \cdot \mathbf{s}_2 \rangle = -\frac{3}{4}$
- $S=1$ (vector): $\langle \mathbf{s}_1 \cdot \mathbf{s}_2 \rangle = +\frac{1}{4}$

### Baryon ($qqq$)

Color singlet (antisymmetric): $\frac{\boldsymbol{\lambda}_i}{2} \cdot \frac{\boldsymbol{\lambda}_j}{2} = -\frac{2}{3}$ for each pair

## Jacobi Coordinates for Three-Body Systems

For a 3-body system with particles 1, 2, 3 having masses $m_1$, $m_2$, $m_3$, there are three equivalent Jacobi coordinate sets corresponding to different clusterings.

### Definition of Jacobi Coordinate Sets

**Set $c=1$: clustering (23)1**

$$
\boldsymbol{\rho}_1 = \mathbf{r}_2 - \mathbf{r}_3, \quad \boldsymbol{\lambda}_1 = \frac{m_2 \mathbf{r}_2 + m_3 \mathbf{r}_3}{m_2 + m_3} - \mathbf{r}_1
$$

**Set $c=2$: clustering (31)2**

$$
\boldsymbol{\rho}_2 = \mathbf{r}_3 - \mathbf{r}_1, \quad \boldsymbol{\lambda}_2 = \frac{m_3 \mathbf{r}_3 + m_1 \mathbf{r}_1}{m_3 + m_1} - \mathbf{r}_2
$$

**Set $c=3$: clustering (12)3**

$$
\boldsymbol{\rho}_3 = \mathbf{r}_1 - \mathbf{r}_2, \quad \boldsymbol{\lambda}_3 = \frac{m_1 \mathbf{r}_1 + m_2 \mathbf{r}_2}{m_1 + m_2} - \mathbf{r}_3
$$

### Reduced Masses

For each coordinate set $c$, the reduced masses are:

$$
\mu_{\rho,c} = \frac{m_i m_j}{m_i + m_j}, \quad \mu_{\lambda,c} = \frac{(m_i + m_j) m_k}{m_i + m_j + m_k}
$$

where $(i,j)$ are the clustered pair and $k$ is the spectator.

For example, in set $c=3$ (12)3:

$$
\mu_\rho = \frac{m_1 m_2}{m_1 + m_2}, \quad \mu_\lambda = \frac{(m_1 + m_2) m_3}{M}
$$

where $M = m_1 + m_2 + m_3$ is the total mass.

### Interparticle Distances in Jacobi Coordinates

The key challenge is expressing all interparticle distances $r_{ij} = |\mathbf{r}_i - \mathbf{r}_j|$ in terms of a single Jacobi set.

In set $c=3$ (12)3:

$$
\mathbf{r}_1 - \mathbf{r}_2 = \boldsymbol{\rho}_3
$$

$$
\mathbf{r}_1 - \mathbf{r}_3 = \frac{m_2}{m_1+m_2} \boldsymbol{\rho}_3 + \boldsymbol{\lambda}_3
$$

$$
\mathbf{r}_2 - \mathbf{r}_3 = -\frac{m_1}{m_1+m_2} \boldsymbol{\rho}_3 + \boldsymbol{\lambda}_3
$$

Therefore:

$$
r_{12}^2 = \rho_3^2
$$

$$
r_{13}^2 = \left(\frac{m_2}{M_{12}}\right)^2 \rho_3^2 + \lambda_3^2 + 2\frac{m_2}{M_{12}} \rho_3 \lambda_3 \cos\theta
$$

$$
r_{23}^2 = \left(\frac{m_1}{M_{12}}\right)^2 \rho_3^2 + \lambda_3^2 - 2\frac{m_1}{M_{12}} \rho_3 \lambda_3 \cos\theta
$$

where $M_{12} = m_1 + m_2$ and $\theta$ is the angle between $\boldsymbol{\rho}_3$ and $\boldsymbol{\lambda}_3$.

**The Problem:** The cross term $\rho \lambda \cos\theta$ couples the two coordinates, making matrix element evaluation non-trivial for pairs 13 and 23.

## Raynal-Revai Transformation

The Raynal-Revai (RR) transformation relates Jacobi coordinates defined with different clusterings. This is essential for computing matrix elements of all pair potentials.

### Mass-Scaled Jacobi Coordinates

To simplify the kinetic energy, we introduce mass-scaled coordinates:

$$
\mathbf{x}_c = \sqrt{\mu_{\rho,c}} \, \boldsymbol{\rho}_c, \quad \mathbf{y}_c = \sqrt{\mu_{\lambda,c}} \, \boldsymbol{\lambda}_c
$$

In these coordinates, the kinetic energy becomes:

$$
T = \frac{1}{2M} P_{\text{CM}}^2 + \frac{1}{2}\left( \dot{\mathbf{x}}^2 + \dot{\mathbf{y}}^2 \right)
$$

### Orthogonal Transformation

The transformation from set $c$ to set $c'$ is an orthogonal rotation:

$$
\begin{pmatrix} \mathbf{x}_{c'} \\ \mathbf{y}_{c'} \end{pmatrix} =
\begin{pmatrix} \cos\varphi_{cc'} & \sin\varphi_{cc'} \\ -\sin\varphi_{cc'} & \cos\varphi_{cc'} \end{pmatrix}
\begin{pmatrix} \mathbf{x}_c \\ \mathbf{y}_c \end{pmatrix}
$$

### Transformation Angles

The transformation angle $\varphi_{cc'}$ depends on the masses. For the transformation from set 3 (12)3 to set 1 (23)1:

$$
\cos\varphi_{31} = -\sqrt{\frac{m_2^2 M}{(m_1+m_2)(m_2+m_3)(m_1+m_3)}}
$$

$$
\sin\varphi_{31} = \sqrt{\frac{m_1 m_3 M}{(m_1+m_2)(m_2+m_3)(m_1+m_3)}}
$$

**General formula:** The $d$-coefficients are defined as:

$$
d_k^{(ij)} = \sqrt{\frac{m_k M_{ij}}{\mu_{ij} M}}, \quad M_{ij} = m_i + m_j
$$

The transformation angle from (12)3 to (23)1 is:

$$
\cos\varphi = -\frac{m_2}{\sqrt{M_{12} M_{23}}}, \quad \sin\varphi = \sqrt{\frac{m_1 m_3}{M_{12} M_{23}}}
$$

### Transformation of Jacobi Coordinates

In terms of the original (non-mass-scaled) coordinates:

$$
\boldsymbol{\rho}_1 = -C_\rho \boldsymbol{\rho}_3 - S_\rho \boldsymbol{\lambda}_3
$$

$$
\boldsymbol{\lambda}_1 = S_\lambda \boldsymbol{\rho}_3 - C_\lambda \boldsymbol{\lambda}_3
$$

where the coefficients are:

$$
C_\rho = \sqrt{\frac{\mu_{\rho,3}}{\mu_{\rho,1}}} \cos\varphi, \quad S_\rho = \sqrt{\frac{\mu_{\lambda,3}}{\mu_{\rho,1}}} \sin\varphi
$$

$$
S_\lambda = \sqrt{\frac{\mu_{\rho,3}}{\mu_{\lambda,1}}} \sin\varphi, \quad C_\lambda = \sqrt{\frac{\mu_{\lambda,3}}{\mu_{\lambda,1}}} \cos\varphi
$$

## Raynal-Revai Coefficients for Gaussian Basis

### Transformation of Gaussian Wave Functions

For the product Gaussian basis:

$$
\Phi_{nm}(\boldsymbol{\rho}, \boldsymbol{\lambda}) = \phi_n(\rho) \phi_m(\lambda) = N_n N_m e^{-\nu_n \rho^2 - \nu_m \lambda^2}
$$

Under the RR transformation, this becomes:

$$
\Phi_{nm}^{(3)}(\boldsymbol{\rho}_3, \boldsymbol{\lambda}_3) = \sum_{n'm'} \mathcal{R}_{nm,n'm'}(\varphi) \Phi_{n'm'}^{(1)}(\boldsymbol{\rho}_1, \boldsymbol{\lambda}_1)
$$

### Explicit Formula for l=0 (S-wave)

For s-wave ($l_\rho = l_\lambda = 0$), the RR coefficient has a closed form. When transforming a Gaussian:

$$
e^{-\nu_\rho \rho_3^2 - \nu_\lambda \lambda_3^2}
$$

Using the inverse transformation:

$$
\rho_3 = -C_\rho^{-1}(\rho_1 + S_\rho C_\lambda^{-1} \lambda_1)
$$

The transformed Gaussian becomes:

$$
\exp\left[ -A \rho_1^2 - B \lambda_1^2 - 2C \rho_1 \lambda_1 \cos\theta' \right]
$$

where:

$$
A = \nu_\rho C_\rho^2 + \nu_\lambda S_\lambda^2
$$

$$
B = \nu_\rho S_\rho^2 + \nu_\lambda C_\lambda^2
$$

$$
C = (\nu_\lambda C_\lambda S_\lambda - \nu_\rho C_\rho S_\rho)
$$

### Cross-Term and Angular Dependence

The cross term $C \rho_1 \lambda_1 \cos\theta'$ means the transformed wave function is **not separable** in the new coordinates. This is the fundamental complication.

**Special case:** When $C = 0$, the wave function remains separable. This occurs when:

$$
\nu_\lambda C_\lambda S_\lambda = \nu_\rho C_\rho S_\rho
$$

or equivalently:

$$
\frac{\nu_\lambda}{\nu_\rho} = \frac{C_\rho S_\rho}{C_\lambda S_\lambda}
$$

## Practical Approaches for Three-Body Calculations

### Approach 1: Direct Angular Integration

For s-wave Gaussians, integrate over the angle $\theta$ analytically:

$$
\langle \Phi_{nm} | V(r_{13}) | \Phi_{n'm'} \rangle = \int_0^\infty \int_0^\infty \int_0^\pi V(r_{13}) |\Phi_{nm}|^2 \rho^2 \lambda^2 \sin\theta \, d\rho \, d\lambda \, d\theta
$$

For Gaussian potentials $V(r) = e^{-\alpha r^2}$, the angular integral gives:

$$
\int_0^\pi e^{-\alpha \cdot 2ab\rho\lambda\cos\theta} \sin\theta \, d\theta = \frac{\sinh(2\alpha ab\rho\lambda)}{\alpha ab\rho\lambda}
$$

This leads to modified matrix elements involving error functions.

### Approach 2: RR Transformation Method

1. **For $V_{12}(r_{12})$:** Compute directly in set (12)3 where $r_{12} = \rho_3$.

2. **For $V_{23}(r_{23})$:**
   - Transform basis to set (23)1 where $r_{23} = \rho_1$
   - Compute matrix elements in set 1
   - Transform back using RR coefficients

3. **For $V_{13}(r_{13})$:**
   - Transform basis to set (13)2 where $r_{13} = \rho_2$
   - Compute matrix elements in set 2
   - Transform back

### Matrix Element Transformation Formula

$$
\langle \Phi_{nm}^{(3)} | V | \Phi_{n'm'}^{(3)} \rangle = \sum_{\alpha\beta\gamma\delta} \mathcal{R}_{nm,\alpha\beta} \langle \Phi_{\alpha\beta}^{(1)} | V | \Phi_{\gamma\delta}^{(1)} \rangle \mathcal{R}_{n'm',\gamma\delta}^*
$$

For the GEM with $N$ basis functions, this involves $N^4$ terms, but can be computed efficiently using matrix operations.

## Equal Mass Simplification

For equal masses ($m_1 = m_2 = m_3 = m$), the formulas simplify significantly:

$$
\mu_\rho = \frac{m}{2}, \quad \mu_\lambda = \frac{2m}{3}
$$

The transformation angle becomes:

$$
\cos\varphi = -\frac{1}{2}, \quad \sin\varphi = \frac{\sqrt{3}}{2}, \quad \varphi = \frac{2\pi}{3} = 120°
$$

The interparticle distances are:

$$
r_{12} = \rho, \quad r_{13} = \sqrt{\frac{\rho^2}{4} + \lambda^2 + \rho\lambda\cos\theta}, \quad r_{23} = \sqrt{\frac{\rho^2}{4} + \lambda^2 - \rho\lambda\cos\theta}
$$

The three Jacobi sets are related by $120°$ rotations, reflecting the $S_3$ permutation symmetry of identical particles.

## Summary: Algorithm for Baryon Calculations

1. **Choose a reference Jacobi set** (e.g., (12)3)

2. **Build the kinetic energy matrix** (separable in any Jacobi set):
   $$T_{ij} = T_\rho(i_\rho, j_\rho) S_\lambda(i_\lambda, j_\lambda) + S_\rho(i_\rho, j_\rho) T_\lambda(i_\lambda, j_\lambda)$$

3. **Build the potential matrix for pair 12** directly:
   $$V_{12} = V_\rho(i_\rho, j_\rho) S_\lambda(i_\lambda, j_\lambda)$$

4. **Build $V_{23}$ and $V_{13}$** using either:
   - Angular integration approach, or
   - RR transformation to diagonal coordinate set

5. **Add rest masses and solve** the generalized eigenvalue problem:
   $$(T + V_{12} + V_{13} + V_{23}) \mathbf{c} = E \, S \, \mathbf{c}$$

## Generalization to N-Body Systems

For tetraquarks (4-body) and pentaquarks (5-body), we need a systematic framework to handle the increased complexity. The key insight is that the Raynal-Revai transformation generalizes to **orthogonal rotations in (N-1)-dimensional space**.

### Jacobi Trees: Systematic Clustering

For an N-body system, the Jacobi coordinates are defined by a **binary tree structure** (Jacobi tree) that specifies the clustering order.

**Definition:** A Jacobi tree for N particles is a binary tree where:
- Each leaf represents a particle
- Each internal node represents the center of mass of its children
- The N-1 edges correspond to the N-1 Jacobi coordinates

**Example for 4-body:**

```
Tree A: ((12)(34))          Tree B: (((12)3)4)

       R                           λ₂
      / \                         / \
    ρ₁   ρ₂                      λ₁   4
    /\    /\                     /\
   1  2  3  4                   ρ   3
                               / \
                              1   2
```

Tree A (diquark-antidiquark):
- $\boldsymbol{\rho}_1 = \mathbf{r}_1 - \mathbf{r}_2$
- $\boldsymbol{\rho}_2 = \mathbf{r}_3 - \mathbf{r}_4$
- $\mathbf{R} = \text{CM}(12) - \text{CM}(34)$

Tree B (sequential):
- $\boldsymbol{\rho} = \mathbf{r}_1 - \mathbf{r}_2$
- $\boldsymbol{\lambda}_1 = \text{CM}(12) - \mathbf{r}_3$
- $\boldsymbol{\lambda}_2 = \text{CM}(123) - \mathbf{r}_4$

### General Properties

For an N-body system:

| Property | Formula |
|----------|---------|
| Number of Jacobi coordinates | $N - 1$ |
| Number of particle pairs | $\frac{N(N-1)}{2}$ |
| Number of distinct Jacobi trees | $(2N-3)!! = 1 \times 3 \times 5 \times \cdots \times (2N-3)$ |
| Dimension of kinematic space | $3(N-1)$ for 3D |

**Examples:**

| System | N | Jacobi coords | Pairs | Trees |
|--------|---|---------------|-------|-------|
| Meson | 2 | 1 | 1 | 1 |
| Baryon | 3 | 2 | 3 | 3 |
| Tetraquark | 4 | 3 | 6 | 15 |
| Pentaquark | 5 | 4 | 10 | 105 |

### Mass-Scaled Coordinates and Kinematic Rotations

**Definition:** The mass-scaled Jacobi coordinates are:

$$
\boldsymbol{\xi}_i = \sqrt{\mu_i} \, \boldsymbol{\eta}_i
$$

where $\boldsymbol{\eta}_i$ is the $i$-th Jacobi coordinate and $\mu_i$ is its reduced mass.

In these coordinates, the kinetic energy is diagonal:

$$
T = \frac{1}{2M} P_{\text{CM}}^2 + \frac{1}{2} \sum_{i=1}^{N-1} \dot{\boldsymbol{\xi}}_i^2
$$

**Key Theorem:** The transformation between any two Jacobi trees is an **orthogonal rotation** in the $(N-1)$-dimensional space of mass-scaled coordinates:

$$
\boldsymbol{\xi}' = \mathbf{O} \boldsymbol{\xi}
$$

where $\mathbf{O}$ is an $(N-1) \times (N-1)$ orthogonal matrix ($\mathbf{O}^T \mathbf{O} = \mathbf{I}$).

### Computing the Transformation Matrix

The transformation matrix $\mathbf{O}$ can be computed from the **T-coefficients**:

$$
\mathbf{r}_i - \mathbf{R}_{\text{CM}} = \sum_{k=1}^{N-1} T_{ik} \boldsymbol{\eta}_k
$$

where $\mathbf{R}_{\text{CM}}$ is the center of mass.

For two Jacobi trees with T-coefficient matrices $\mathbf{T}$ and $\mathbf{T}'$:

$$
\mathbf{O} = \mathbf{D}'^{-1/2} (\mathbf{T}')^{-1} \mathbf{T} \, \mathbf{D}^{1/2}
$$

where $\mathbf{D} = \text{diag}(\mu_1, \mu_2, \ldots, \mu_{N-1})$.

## Tetraquark (4-Body) Jacobi Coordinates

### Diquark-Antidiquark Clustering ((12)(34))

For a tetraquark $q_1 q_2 \bar{q}_3 \bar{q}_4$, the natural clustering is:

$$
\boldsymbol{\rho}_1 = \mathbf{r}_1 - \mathbf{r}_2 \quad \text{(quark-quark)}
$$

$$
\boldsymbol{\rho}_2 = \mathbf{r}_3 - \mathbf{r}_4 \quad \text{(antiquark-antiquark)}
$$

$$
\mathbf{R} = \frac{m_1 \mathbf{r}_1 + m_2 \mathbf{r}_2}{m_1 + m_2} - \frac{m_3 \mathbf{r}_3 + m_4 \mathbf{r}_4}{m_3 + m_4} \quad \text{(diquark-antidiquark)}
$$

### Reduced Masses

$$
\mu_{\rho_1} = \frac{m_1 m_2}{m_1 + m_2}, \quad \mu_{\rho_2} = \frac{m_3 m_4}{m_3 + m_4}
$$

$$
\mu_R = \frac{(m_1 + m_2)(m_3 + m_4)}{M}, \quad M = m_1 + m_2 + m_3 + m_4
$$

### Interparticle Distances

In the ((12)(34)) Jacobi set:

$$
\mathbf{r}_1 - \mathbf{r}_2 = \boldsymbol{\rho}_1 \quad \text{(diagonal)}
$$

$$
\mathbf{r}_3 - \mathbf{r}_4 = \boldsymbol{\rho}_2 \quad \text{(diagonal)}
$$

$$
\mathbf{r}_1 - \mathbf{r}_3 = \frac{m_2}{M_{12}} \boldsymbol{\rho}_1 - \frac{m_4}{M_{34}} \boldsymbol{\rho}_2 + \mathbf{R}
$$

$$
\mathbf{r}_1 - \mathbf{r}_4 = \frac{m_2}{M_{12}} \boldsymbol{\rho}_1 + \frac{m_3}{M_{34}} \boldsymbol{\rho}_2 + \mathbf{R}
$$

$$
\mathbf{r}_2 - \mathbf{r}_3 = -\frac{m_1}{M_{12}} \boldsymbol{\rho}_1 - \frac{m_4}{M_{34}} \boldsymbol{\rho}_2 + \mathbf{R}
$$

$$
\mathbf{r}_2 - \mathbf{r}_4 = -\frac{m_1}{M_{12}} \boldsymbol{\rho}_1 + \frac{m_3}{M_{34}} \boldsymbol{\rho}_2 + \mathbf{R}
$$

where $M_{12} = m_1 + m_2$ and $M_{34} = m_3 + m_4$.

### Classification of Pair Potentials

| Pair | Jacobi dependence | Direct in ((12)(34))? |
|------|-------------------|----------------------|
| 1-2 | $\rho_1$ only | Yes |
| 3-4 | $\rho_2$ only | Yes |
| 1-3 | $\rho_1, \rho_2, R$ | No - needs transformation |
| 1-4 | $\rho_1, \rho_2, R$ | No - needs transformation |
| 2-3 | $\rho_1, \rho_2, R$ | No - needs transformation |
| 2-4 | $\rho_1, \rho_2, R$ | No - needs transformation |

### Alternative Clustering: ((13)(24))

For pairs 1-3 and 2-4, use the Jacobi set:

$$
\boldsymbol{\rho}'_1 = \mathbf{r}_1 - \mathbf{r}_3, \quad \boldsymbol{\rho}'_2 = \mathbf{r}_2 - \mathbf{r}_4
$$

$$
\mathbf{R}' = \text{CM}(13) - \text{CM}(24)
$$

### The 3×3 Transformation Matrix

The transformation from ((12)(34)) to ((13)(24)) is a 3D rotation:

$$
\begin{pmatrix} \boldsymbol{\xi}'_1 \\ \boldsymbol{\xi}'_2 \\ \boldsymbol{\Xi}' \end{pmatrix} =
\mathbf{O}_{4}
\begin{pmatrix} \boldsymbol{\xi}_1 \\ \boldsymbol{\xi}_2 \\ \boldsymbol{\Xi} \end{pmatrix}
$$

where $\mathbf{O}_4$ is a 3×3 orthogonal matrix determined by the masses.

**For equal masses** ($m_1 = m_2 = m_3 = m_4 = m$):

$$
\mathbf{O}_4 = \frac{1}{2} \begin{pmatrix} 1 & 1 & \sqrt{2} \\ 1 & 1 & -\sqrt{2} \\ -\sqrt{2} & \sqrt{2} & 0 \end{pmatrix}
$$

## Pentaquark (5-Body) Jacobi Coordinates

### Baryon-Meson Clustering (((12)3)(45))

For a pentaquark $q_1 q_2 q_3 \bar{q}_4 q_5$:

$$
\boldsymbol{\rho} = \mathbf{r}_1 - \mathbf{r}_2 \quad \text{(quark pair in baryon)}
$$

$$
\boldsymbol{\lambda} = \text{CM}(12) - \mathbf{r}_3 \quad \text{(third quark to pair)}
$$

$$
\boldsymbol{\sigma} = \mathbf{r}_4 - \mathbf{r}_5 \quad \text{(meson internal)}
$$

$$
\mathbf{R} = \text{CM}(123) - \text{CM}(45) \quad \text{(baryon-meson)}
$$

### The 4×4 Transformation Matrix

Transformations between 5-body Jacobi trees are 4D rotations:

$$
\boldsymbol{\xi}' = \mathbf{O}_5 \boldsymbol{\xi}
$$

where $\mathbf{O}_5$ is a 4×4 orthogonal matrix.

### Number of Pair Potentials

With 5 particles, there are $\binom{5}{2} = 10$ pair potentials:
- Some pairs are diagonal in the chosen Jacobi set
- Others require transformation to different Jacobi trees

## Systematic Algorithm for N-Body Calculations

### Step 1: Choose Reference Jacobi Tree

Select a Jacobi tree $\mathcal{T}_0$ that makes the most important pairs diagonal.

**Heuristic:** Choose clustering that diagonalizes the strongest interactions (e.g., diquarks for tetraquarks).

### Step 2: Enumerate All Pairs

List all $\frac{N(N-1)}{2}$ particle pairs and classify:
- **Diagonal pairs:** $r_{ij}$ equals one Jacobi coordinate
- **Non-diagonal pairs:** $r_{ij}$ depends on multiple coordinates

### Step 3: Find Optimal Jacobi Trees for Each Pair

For each non-diagonal pair $(i,j)$, find a Jacobi tree $\mathcal{T}_{ij}$ where $r_{ij}$ is diagonal.

### Step 4: Compute Transformation Matrices

For each required transformation $\mathcal{T}_0 \to \mathcal{T}_{ij}$, compute the orthogonal matrix $\mathbf{O}_{ij}$.

### Step 5: Build the Hamiltonian

**Kinetic energy:** (same in all Jacobi trees due to orthogonal transformation)

$$
T = \sum_{k=1}^{N-1} \frac{(\hbar c)^2}{2\mu_k} \cdot \frac{6\nu_{n_k}\nu_{m_k}}{\nu_{n_k}+\nu_{m_k}} \langle\phi_{n_k}|\phi_{m_k}\rangle
$$

**Potential energy for diagonal pair:**

$$
V_{ij}^{\text{diag}} = V(\rho_k) \otimes \mathbf{S}_{\text{other}}
$$

where $\mathbf{S}_{\text{other}}$ is the overlap matrix for all other coordinates.

**Potential energy for non-diagonal pair:**

$$
V_{ij} = \mathbf{O}_{ij}^T V_{ij}^{\text{diag}} \mathbf{O}_{ij}
$$

### Step 6: Tensor Product Basis

The full basis is a tensor product:

$$
|\Phi_{\mathbf{n}}\rangle = |\phi_{n_1}\rangle \otimes |\phi_{n_2}\rangle \otimes \cdots \otimes |\phi_{n_{N-1}}\rangle
$$

Total basis size: $N_{\text{basis}}^{N-1}$ (grows exponentially with particle number).

**Practical limits:**
- Tetraquark: $20^3 = 8000$ basis states
- Pentaquark: $15^4 = 50625$ basis states

## Gaussian Transformation Under Orthogonal Rotation

### The Fundamental Formula

A product Gaussian in the original coordinates:

$$
\Phi(\boldsymbol{\xi}) = \exp\left( -\sum_{i=1}^{N-1} \nu_i \xi_i^2 \right) = \exp\left( -\boldsymbol{\xi}^T \mathbf{N} \boldsymbol{\xi} \right)
$$

where $\mathbf{N} = \text{diag}(\nu_1, \ldots, \nu_{N-1})$.

Under the orthogonal transformation $\boldsymbol{\xi}' = \mathbf{O} \boldsymbol{\xi}$:

$$
\Phi(\boldsymbol{\xi}) = \exp\left( -\boldsymbol{\xi}'^T \mathbf{O} \mathbf{N} \mathbf{O}^T \boldsymbol{\xi}' \right) = \exp\left( -\boldsymbol{\xi}'^T \mathbf{N}' \boldsymbol{\xi}' \right)
$$

where $\mathbf{N}' = \mathbf{O} \mathbf{N} \mathbf{O}^T$.

### Non-Separability

**Key observation:** Unless $\mathbf{N}$ is proportional to the identity (all $\nu_i$ equal), the matrix $\mathbf{N}'$ is generally **not diagonal**.

This means the transformed Gaussian is **not separable** in the new coordinates — it contains cross terms like $\xi'_i \xi'_j$.

### Expansion in New Basis

The non-separable Gaussian must be expanded in the separable basis of the new coordinate system:

$$
\exp\left( -\boldsymbol{\xi}'^T \mathbf{N}' \boldsymbol{\xi}' \right) = \sum_{\mathbf{n}'} c_{\mathbf{n}'} \prod_{i=1}^{N-1} \phi_{n'_i}(\xi'_i)
$$

The coefficients $c_{\mathbf{n}'}$ are the **generalized Raynal-Revai coefficients**.

### Efficient Computation Strategy

Instead of computing RR coefficients explicitly, use:

1. **Matrix representation:** Work with matrices in the tensor product space
2. **Sparse structure:** Exploit the fact that $\mathbf{O}$ couples only certain basis states
3. **Numerical transformation:** Transform the full Hamiltonian matrix numerically

## Summary: Scaling and Complexity

| System | Jacobi dim | Pairs | Basis size ($N_b=15$) | Matrix elements |
|--------|-----------|-------|----------------------|-----------------|
| Meson | 1 | 1 | 15 | 225 |
| Baryon | 2 | 3 | 225 | ~50,000 |
| Tetraquark | 3 | 6 | 3,375 | ~11 million |
| Pentaquark | 4 | 10 | 50,625 | ~2.5 billion |

**Computational strategies for large systems:**
1. Use symmetry (identical particles, isospin) to reduce basis size
2. Stochastic variational method instead of full diagonalization
3. Select important configurations based on physical intuition
4. Use correlated Gaussians with optimized nonlinear parameters
