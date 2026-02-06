#include <iostream>
#include <iomanip>
#include "gem/hadrons/baryon.hpp"
#include "gem/core/constants.hpp"

using namespace gem;

int main() {
    std::cout << "\n=== Charmed Baryon Masses ===\n\n";
    
    // Lambda_c: udc, S=1/2
    Baryon lambda_c = BaryonFactory::lambda_c();
    Real lambda_c_mass = lambda_c.calculateMass() * constants::GEV_TO_MEV;
    std::cout << "Lambda_c (udc, S=1/2): " << std::fixed << std::setprecision(1) 
              << lambda_c_mass << " MeV (exp: 2286 MeV)";
    if (lambda_c.isValid()) std::cout << "  [OK]";
    std::cout << "\n";
    
    // Sigma_c: uuc, S=1/2
    Baryon sigma_c = BaryonFactory::sigma_c();
    Real sigma_c_mass = sigma_c.calculateMass() * constants::GEV_TO_MEV;
    std::cout << "Sigma_c (uuc, S=1/2):  " << std::fixed << std::setprecision(1) 
              << sigma_c_mass << " MeV (exp: 2453 MeV)";
    if (sigma_c.isValid()) std::cout << "  [OK]";
    std::cout << "\n";
    
    // Sigma_c*: uuc, S=3/2
    Baryon sigma_c_star(QuarkFlavor::Up, QuarkFlavor::Up, QuarkFlavor::Charm, 3);
    Real sigma_c_star_mass = sigma_c_star.calculateMass() * constants::GEV_TO_MEV;
    std::cout << "Sigma_c* (uuc, S=3/2): " << std::fixed << std::setprecision(1) 
              << sigma_c_star_mass << " MeV (exp: 2518 MeV)";
    if (sigma_c_star.isValid()) std::cout << "  [OK]";
    std::cout << "\n";
    
    return 0;
}
