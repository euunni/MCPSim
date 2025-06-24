#include <iostream>
#include <filesystem>
#include "MCPSimulation.h"
#include "MCPConfig.h"
#include "TError.h"

int main(int argc, char* argv[]) {
    gErrorIgnoreLevel = kError;
    
    try {
        // Get absolute path to config file
        std::filesystem::path config_path = std::filesystem::absolute("/home/jangh/MCPSim/v3_track/MCPSim/config/config.txt");
        std::cout << "Loading config from: " << config_path << std::endl;
        
        auto& config = MCPSim::Config::getInstance();
        config.loadFromFile(config_path.string());
        
        // Print some config values to verify
        std::cout << "Config values:" << std::endl;
        std::cout << "diff_pot = " << config.get("diff_pot") << std::endl;
        std::cout << "m = " << config.get("m") << std::endl;
        std::cout << "x0 = " << config.get("x0") << std::endl;
        std::cout << "x1 = " << config.get("x1") << std::endl;
        std::cout << "x2 = " << config.get("x2") << std::endl;
        std::cout << "c_c = " << config.get("c_c") << std::endl;
        std::cout << "c_e = " << config.get("c_e") << std::endl;
        std::cout << "alpha = " << config.get("alpha") << std::endl;
        std::cout << "gama_ts0 = " << config.get("gama_ts0") << std::endl;

        std::string outFile = "/home/jangh/MCPSim/v3_track/output/400µm600V_3.5_0.5_saturation_0.8_each_5ps_" + 
                             std::string(argv[1]) + "_" + std::string(argv[2]) + ".root";

        MCPSim::Simulation sim;
        std::cout << "\nStarting simulation with 4.21 eV photoelectron..." << std::endl;
        auto results = sim.Run(4.21);  // 4.21 eV photoelectron
        sim.Save(results, 4.21, outFile);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 
