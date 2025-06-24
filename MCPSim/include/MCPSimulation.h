#ifndef MCP_SIMULATION_H
#define MCP_SIMULATION_H

#include <vector>
#include <memory>
#include <string>
#include "MCPPhysics.h"
#include "MCPEventData.h"
#include "MCPRootManager.h"
#include <unordered_map>

namespace MCPSim {

class Simulation {
public:
    Simulation();
    
    std::vector<Matrix3x3> Run(double initial_energy);
    void Save(const std::vector<Matrix3x3>& results, double initial_energy, const std::string& rootfile);

private:
    mcp::Event ConvertEvent(const std::vector<Matrix3x3>& results, double initial_energy);
    
    // Function to record the trajectory of electrons outside the pore
    void TrackElectronOutsidePore(const Matrix3x3& start, const Matrix3x3& end, int trackID, double cts);
    
    // Create new electron (add track)
    int CreateElectron(int parentID, float time,
                     float posX, float posY, float posZ,
                     float velX, float velY, float velZ,
                     float energy, int procType);
    
    // Update electron position (add step)
    void AddElectronStep(int trackID, float time,
                       float posX, float posY, float posZ,
                       float velX, float velY, float velZ,
                       float energy);
    
    // Electron termination (anode hit or end)
    void FinalizeElectron(int trackID, int status, float time,
                        float posX, float posY, float posZ,
                        float velX, float velY, float velZ,
                        float energy);

    std::unique_ptr<Physics> physics;
    double x0, x1, x2;
    double q, limite, I_strip, c_c, m;
    int nextTrackID_ = 0;                             // unique track id counter
    mcp::Track tracks_;                               // Actual electron information
    mcp::Step steps_;                                 // Trajectory point information
    std::unordered_map<int, int> trackIDMap_;         // External trackID -> internal trackID mapping
};

struct SimElectron {
    Matrix3x3 M;       // 3x3 state matrix (pos, vel, misc)
    int trackID;       // unique identifier
    int parentID;      // parent trackID (-1 for primary)
};

} // namespace MCPSim 

#endif // MCP_SIMULATION_H 