#ifndef MCP_SIMULATION_H
#define MCP_SIMULATION_H

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <iostream>

#include "MCPPhysics.h"
#include "MCPEventData.h"
#include "MCPRootManager.h"

namespace MCPSim {

using Matrix3x3 = Eigen::Matrix<double,3,3>;

// ---------- Output control (Track / Node / Step) ----------
enum class OutputLevel { kTrack = 0, kNode = 1, kStep = 2 };
// Return current output level read from configuration (default = kNode)
OutputLevel GetOutputLevel();

class Simulation {
public:
    Simulation();
    
    // main entry
    std::vector<Matrix3x3> Run(double initial_energy);
    void Save(const std::vector<Matrix3x3>& results,
              double initial_energy,
              const std::string& rootfile);
    mcp::Event ConvertEvent(const std::vector<Matrix3x3>& results,
                            double initial_energy);

    // track helpers
    int  CreateElectron(int parentID, float time,
                        float x, float y, float z,
                        float vx, float vy, float vz,
                     float energy, int procType);
    void AddElectronStep(int trackID, float time,
                         float x, float y, float z,
                         float vx, float vy, float vz,
                       float energy);

private:
    // helpers
    void TrackElectronOutsidePore(const Matrix3x3& start,
                                  const Matrix3x3& end,
                                  int trackID, double cts);
    void RecordNode(int trackID, const Matrix3x3& M);   // <— NEW helper: store single point depending on level
    void FinalizeElectron(int trackID, int status, float time,
                          float x, float y, float z,
                          float vx, float vy, float vz,
                        float energy);

    // ---------- members ----------
    std::unique_ptr<Physics> physics;
    mcp::Track              tracks_;
    mcp::Step               steps_;
    std::unordered_map<int,int> trackIDMap_;
    int   nextTrackID_{1};

    // electrons reaching anode
    std::vector<Matrix3x3> anode_hits_;

    // statistics: secondaries (parentID!=-1) that exit MCP-1 into GAP-1
    int nSecOutMCP1_ = 0;

    // --- options ---
    bool recordGapSteps_ = true;   // if false, TrackElectronOutsidePore() becomes no-op

public:
    void SetRecordGapSteps(bool on){ recordGapSteps_ = on; }
};

} // namespace MCPSim 
#endif
