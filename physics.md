## Gaussian expansion method (GEM)
* GEM is a method to solve the few-body problem by expand the wave funciton with the following gaussian wave functions
$$
\psi_{lm}(\vec{r})	=\sum_{n=1}c_{nl}\phi_{nlm}^{G}(\vec{r})\\
\phi_{nlm}^{G}(\vec{r})	=\phi_{nl}^{G}(r)Y_{lm}(\hat{r})\\
\phi_{nl}^{G}	=N_{nl}r^{l}e^{-\nu_{n}r^{2}}
$$
* In this work, only l=0 basis is used.
* It is shown that the basis with $\nu_n$ appear as the geometric series is very effeicient.
* the advantage of GEM is the matrix element could be solved analytically.
## Quark potential model
* The mimimal quark model including three part: clomb interaction, hyperfine intraction and linear confinment interaction
$$
V_{clomb}=\frac{\alpha_{s}}{r}\frac{\lambda_{i}}{2}\frac{\lambda_{j}}{2}\\
V_{hyp}(r)=-\frac{8\pi\alpha_{s}}{3m_{i}m_{j}}\frac{\tau^{3}}{\pi^{3/2}}e^{-\tau^{2}r^{2}}\frac{\lambda_{i}}{2}\frac{\lambda_{j}}{2}\mathbf{s}_{i}\cdot\mathbf{s}_{j}\\
V_{cof}(r)=\left(-\frac{3b}{4}r+V_{c}\right)\frac{\lambda_{i}}{2}\frac{\lambda_{j}}{2}
$$ 
* the clomb interaction and confinment only depend or the color electric interaction, where $\lambda_i$ is the Gell-mann matrix in the color space
*  the hyperfine interaction is color magnetic interaction, where the spin $s_i$ will be dpedented. 
* due to the color confinment, the hadrons composed of quark-antiquark should be color singlet. 
* The conventional hadrons composed of the conventional meson composed of $q\bar{q}$ and baryons $qqq$. However, there are also multqiaurk state, for example, tetraquark state, pentaquark state.

## using quark model to calculate the hadron mass
* to get hadron mass, one could slove a nonrelativistic schroding equation of (anti)quarks
* one first write down the complete discrete wave functions of hadrons in the degree of freedome of quarks with the conctraint the totoal color singlent, and fixed total spin and isospin. 
* using these discrete wave function, one can calculate the matrix element of $\lambda_i\cdot \lambda_j$ and $\lambda_i\cdot \lambda_j s_i\cdot s_j$.
* For the spatal wave function, one could expanded it into gaussian bassis in different jacobi. with these wave function, one can caculate the matrix element in the spatial part anlaytical in the GEM.
* The total matrix element could be get by combined the matrix element in the discrete part and spatial part.
* the mass of the multiquark state could be obtained by solve a general eigenvalue probelm.