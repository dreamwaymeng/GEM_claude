#ifndef GEM_CORE_CONFIG_HPP
#define GEM_CORE_CONFIG_HPP

#include "types.hpp"
#include "constants.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <stdexcept>

namespace gem {

class Config {
public:
    Config() {
        setDefaults();
    }

    void setDefaults() {
        basis_.n_basis = constants::default_basis::N_BASIS;
        basis_.r_min = constants::default_basis::R_MIN;
        basis_.r_max = constants::default_basis::R_MAX;

        potential_.alpha_s_coul = constants::default_potential::ALPHA_S_COULOMB;
        potential_.alpha_s_hyp = constants::default_potential::ALPHA_S_HYPERFINE;
        potential_.b = constants::default_potential::B_STRING;
        potential_.V0 = constants::default_potential::V0;
        potential_.hyp_a = constants::default_potential::HYPERFINE_A;
        potential_.hyp_b_exp = constants::default_potential::HYPERFINE_B_EXP;

        quark_masses_[QuarkFlavor::Up] = constants::quark_mass::UP;
        quark_masses_[QuarkFlavor::Down] = constants::quark_mass::DOWN;
        quark_masses_[QuarkFlavor::Strange] = constants::quark_mass::STRANGE;
        quark_masses_[QuarkFlavor::Charm] = constants::quark_mass::CHARM;
        quark_masses_[QuarkFlavor::Bottom] = constants::quark_mass::BOTTOM;
    }

    // Load configuration from a simple key=value file
    bool loadFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }

        std::string line;
        while (std::getline(file, line)) {
            // Skip comments and empty lines
            if (line.empty() || line[0] == '#' || line[0] == '/') {
                continue;
            }

            auto pos = line.find('=');
            if (pos == std::string::npos) continue;

            std::string key = trim(line.substr(0, pos));
            std::string value = trim(line.substr(pos + 1));

            parseKeyValue(key, value);
        }

        return true;
    }

    // Getters
    const BasisParameters& getBasisParameters() const { return basis_; }
    const PotentialParameters& getPotentialParameters() const { return potential_; }

    Real getQuarkMass(QuarkFlavor flavor) const {
        auto it = quark_masses_.find(flavor);
        if (it != quark_masses_.end()) {
            return it->second;
        }
        return constants::getQuarkMass(flavor);
    }

    // Setters
    void setBasisParameters(const BasisParameters& params) { basis_ = params; }
    void setPotentialParameters(const PotentialParameters& params) { potential_ = params; }
    void setQuarkMass(QuarkFlavor flavor, Real mass) { quark_masses_[flavor] = mass; }

    void setAlphaSCoulomb(Real alpha) { potential_.alpha_s_coul = alpha; }
    void setAlphaSHyperfine(Real alpha) { potential_.alpha_s_hyp = alpha; }
    void setStringTension(Real b) { potential_.b = b; }
    void setV0(Real v0) { potential_.V0 = v0; }

    void setBasisSize(int n) { basis_.n_basis = n; }
    void setRRange(Real rmin, Real rmax) {
        basis_.r_min = rmin;
        basis_.r_max = rmax;
    }

private:
    BasisParameters basis_;
    PotentialParameters potential_;
    std::map<QuarkFlavor, Real> quark_masses_;

    static std::string trim(const std::string& str) {
        size_t first = str.find_first_not_of(" \t\n\r");
        if (first == std::string::npos) return "";
        size_t last = str.find_last_not_of(" \t\n\r");
        return str.substr(first, last - first + 1);
    }

    void parseKeyValue(const std::string& key, const std::string& value) {
        Real val = 0.0;
        try {
            val = std::stod(value);
        } catch (...) {
            return;  // Skip invalid values
        }

        if (key == "n_basis") basis_.n_basis = static_cast<int>(val);
        else if (key == "r_min") basis_.r_min = val;
        else if (key == "r_max") basis_.r_max = val;
        else if (key == "alpha_s_coul" || key == "alpha_s") potential_.alpha_s_coul = val;
        else if (key == "alpha_s_hyp") potential_.alpha_s_hyp = val;
        else if (key == "b" || key == "string_tension") potential_.b = val;
        else if (key == "hyp_a") potential_.hyp_a = val;
        else if (key == "hyp_b_exp") potential_.hyp_b_exp = val;
        else if (key == "V0") potential_.V0 = val;
        else if (key == "m_up" || key == "m_u") quark_masses_[QuarkFlavor::Up] = val;
        else if (key == "m_down" || key == "m_d") quark_masses_[QuarkFlavor::Down] = val;
        else if (key == "m_strange" || key == "m_s") quark_masses_[QuarkFlavor::Strange] = val;
        else if (key == "m_charm" || key == "m_c") quark_masses_[QuarkFlavor::Charm] = val;
        else if (key == "m_bottom" || key == "m_b") quark_masses_[QuarkFlavor::Bottom] = val;
    }
};

} // namespace gem

#endif // GEM_CORE_CONFIG_HPP
