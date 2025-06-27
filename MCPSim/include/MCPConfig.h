#pragma once

#include <string>
#include <unordered_map>
#include <memory>

namespace MCPSim {

class Config {
public:
    static Config& getInstance();

    double get(const std::string& param);
    void loadFromFile(const std::string& filename);

private:
    Config() = default;
    ~Config() = default;
    Config(Config&) = delete;
    Config& operator=(Config&) = delete;

    std::unordered_map<std::string, double> config_values_;
};

} // namespace MCPSim 