#include <iostream>
#include <filesystem>
#include "MCPSimulation.h"
#include "MCPConfig.h"
#include "TError.h"

int main(int argc, char* argv[]) {
    gErrorIgnoreLevel = kError;
    
    try {
        // Get absolute path to config file
        std::filesystem::path config_path = std::filesystem::absolute("/home/jangh/MCPSim/v3_geometry_3/MCPSim/MCPSim/config/config.txt");
        std::cout << "Loading config from: " << config_path << std::endl;
        
        auto& config = MCPSim::Config::getInstance();
        config.loadFromFile(config_path.string());

        std::string outFile = "/home/jangh/MCPSim/v3_geometry_3/MCPSim/output/pbs/250723/400µm600V_3.5_0.5_saturation_0.8_each_5ps_" + 
                             std::string(argv[1]) + "_" + std::string(argv[2]) + ".root";

        MCPSim::Simulation sim;
        std::cout << "\nStarting simulation with 4.21 eV photoelectron..." << std::endl;
        auto results = sim.Run(4.21);  // 4.21 eV photoelectron
        sim.Save(results, 4.21, outFile);

        std::cout << "Simulation completed. Results saved to: " << outFile << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 
