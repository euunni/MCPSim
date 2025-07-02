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
    int CreateElectron(int parentID, float time,
                     float posX, float posY, float posZ,
                     float velX, float velY, float velZ,
                     float energy, int procType);
    void AddElectronStep(int trackID, float time,
                       float posX, float posY, float posZ,
                       float velX, float velY, float velZ,
                       float energy);

private:
    mcp::Event ConvertEvent(const std::vector<Matrix3x3>& results, double initial_energy);
    
    // Function to record the trajectory of electrons outside the pore
    void TrackElectronOutsidePore(const Matrix3x3& start, const Matrix3x3& end, int trackID, double cts);
    
    // Electron termination (anode hit or end)
    void FinalizeElectron(int trackID, int status, float time,
                        float posX, float posY, float posZ,
                        float velX, float velY, float velZ,
                        float energy);

    // --- Ensemble-based processing helpers (under migration) ---
    void ProcessMCP1Ensemble(double& cts1, double alpha1, double x0, double x1,
                             double R, double dia, double pas, double m, double E0,
                             int& nextTrackBaseID);

    std::unique_ptr<Physics> physics;
    double x0, x1, x2;
    double q, limite, I_strip, c_c, m;
    int nextTrackID_ = 0;                             // unique track id counter
    mcp::Track tracks_;                               // Actual electron information
    mcp::Step steps_;                                 // Trajectory point information
    std::unordered_map<int, int> trackIDMap_;         // External trackID -> internal trackID mapping

    // --- New data structures for v3_track-style event-driven algorithm ---
    // MCP-1
    std::vector<Matrix3x3> E1_emi_;      // electrons that have a pending wall collision (sorted by t_hit in (2,1))
    std::vector<Matrix3x3> E1_nonemi_;   // electrons inside the pore but currently not destined to hit wall
    std::vector<Matrix3x3> G1_;          // electrons traveling in GAP1 (x1→x2)

    // MCP-2
    std::vector<Matrix3x3> E2_emi_;
    std::vector<Matrix3x3> E2_nonemi_;
    std::vector<Matrix3x3> G2_;          // electrons traveling in GAP2 (x3→x4)
};

struct SimElectron {
    Matrix3x3 M;       // 3x3 state matrix (pos, vel, misc)
    int trackID;       // unique identifier
    int parentID;      // parent trackID (-1 for primary)
    int channel;       // pore index
};

} // namespace MCPSim 

#endif // MCP_SIMULATION_H 