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
    
    // 포어 밖에서 전자의 궤적을 기록하는 함수
    void TrackElectronOutsidePore(const Matrix3x3& start, const Matrix3x3& end, int trackID, double cts);
    
    // 새 전자 생성 (트랙 추가)
    int CreateElectron(int parentID, float time,
                     float posX, float posY, float posZ,
                     float velX, float velY, float velZ,
                     float energy, int procType);
    
    // 전자 위치 업데이트 (스텝 추가)
    void AddElectronStep(int trackID, float time,
                       float posX, float posY, float posZ,
                       float velX, float velY, float velZ,
                       float energy, bool isInteraction = false);
    
    // 전자 종료 (애노드 도달 또는 종료)
    void FinalizeElectron(int trackID, int status, float time,
                        float posX, float posY, float posZ,
                        float velX, float velY, float velZ,
                        float energy);

    std::unique_ptr<Physics> physics;
    double x0, x1, x2;
    double q, limite, I_strip, c_c, m;
    int nextTrackID_ = 0;                             // unique track id counter
    mcp::Track tracks_;                               // 실제 전자 정보
    mcp::Step steps_;                                 // 궤적 포인트 정보
    std::unordered_map<int, int> trackIDMap_;         // 외부 trackID -> 내부 trackID 매핑
};

struct SimElectron {
    Matrix3x3 M;       // 3x3 state matrix (pos, vel, misc)
    int trackID;       // unique identifier
    int parentID;      // parent trackID (-1 for primary)
};

} // namespace MCPSim 

#endif // MCP_SIMULATION_H 