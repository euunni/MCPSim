#include "MCPConfig.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace MCPSim {

Config& Config::getInstance() {
    static Config instance;
    return instance;
}

void Config::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + filename);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        std::istringstream iss(line);
        std::string key;
        double value;
        
        if (std::getline(iss, key, '=')) {
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            
            if (iss >> value) {
                config_values_[key] = value;
            }
        }
    }

    config_values_["c_c"] = (1.6 * config_values_["diff_pot"]) / ((config_values_["x1"] - config_values_["x0"]) * config_values_["m"]);
    config_values_["c_s"] = (1.6 * config_values_["diff_pot"]) / ((config_values_["x1"] - config_values_["x0"]) * config_values_["m"]);
    config_values_["R"] = config_values_["dia"] / 2;
    config_values_["I_strip"] = config_values_["diff_pot"] / config_values_["Resistance"];
    config_values_["c_c2"] = (1.6 * config_values_["diff_pot"]) / ((config_values_["x3"] - config_values_["x2"]) * config_values_["m"]);
    config_values_["c_s2"] = (1.6 * config_values_["diff_pot"]) / ((config_values_["x3"] - config_values_["x2"]) * config_values_["m"]);
    config_values_["R2"] = config_values_["dia"] / 2;
    config_values_["I_strip2"] = config_values_["diff_pot"] / config_values_["Resistance"];
}

double Config::get(const std::string& param) {
    auto it = config_values_.find(param);
    if (it != config_values_.end()) {
        return it->second;
    }
    throw std::runtime_error("Parameter not found in config: " + param);
}

} // namespace MCPSim 