#ifndef MCP_EVENT_DATA_H
#define MCP_EVENT_DATA_H

#include <vector>
#include <cstddef>
#include "TObject.h"

namespace mcp {

class EventInfo : public TObject {
public:
    float initialEnergy;

    EventInfo();
    void Reset();
    
    ClassDef(EventInfo, 2);
};

class ConfigParameters : public TObject {
public:
    // Physical constants
    float diff_pot, m, x0, x1, x2, c_e;
    float pas, dia, alpha, limite, E0, Resistance, q;
    
    // Secondary emission parameters  
    float P_inf, P_1e, W, Ee, p, P_1r, r, Er;
    float e1, e2, r1, r2, t1, t2, t3, t4;
    float gama_ts0, E_ts, s;
    
    ConfigParameters();
    void Reset();
    void LoadFromConfig();  
    
    ClassDef(ConfigParameters, 1);  
};

// 물리 프로세스 타입 상수 정의
enum ProcessType {
    PROCESS_BACKSCATTER = 0,  // 후방 산란
    PROCESS_REDIFFUSED = 1,   // 재확산
    PROCESS_SECONDARY = 2     // 이차 전자 방출
};

// 트랙 클래스: 실제 전자 정보를 저장
class Track : public TObject {
public:
    int nTracks;                   // 총 트랙 수
    std::vector<int> trackID;      // 트랙 ID
    std::vector<int> parentID;     // 부모 트랙 ID
    std::vector<float> birthTime;  // 생성 시간
    std::vector<float> birthPosX, birthPosY, birthPosZ;  // 생성 위치
    std::vector<float> birthVelX, birthVelY, birthVelZ;  // 초기 속도
    std::vector<float> birthEnergy;  // 초기 에너지
    std::vector<int> processType;  // 생성 프로세스 타입
    std::vector<int> isAnode;      // 애노드 도달 여부 (0: 도달하지 않음, 1: 도달함)
    std::vector<float> finalTime;  // 최종 시간 (애노드 도달 또는 종료 시간)
    std::vector<float> finalPosX, finalPosY, finalPosZ;  // 최종 위치
    std::vector<float> finalVelX, finalVelY, finalVelZ;  // 최종 속도
    std::vector<float> finalEnergy;  // 최종 에너지
    std::vector<float> transitTime;  // 이동 시간 (finalTime - birthTime)
    
    Track();
    virtual ~Track();  // 가상 소멸자 추가
    void Reset();
    
    // 새 트랙 추가 (전자 생성)
    int AddTrack(int parentID, float time,
                float posX, float posY, float posZ,
                float velX, float velY, float velZ,
                float energy, int procType);
    
    // 트랙 종료 (애노드 도달 또는 종료)
    void FinalizeTrack(int trackID, int isAnodeHit, float time,
                      float posX, float posY, float posZ,
                      float velX, float velY, float velZ,
                      float energy);
                      
    // 트랙 ID로 인덱스 찾기
    int FindTrackIndex(int trackID) const;
    
    ClassDef(Track, 1);
};

// 스텝 클래스: 전자 궤적 정보를 저장
class Step : public TObject {
public:
    int nSteps;                    // 총 스텝 수
    std::vector<int> trackID;      // 어떤 트랙의 스텝인지
    std::vector<float> time;       // 시간
    std::vector<float> posX, posY, posZ;  // 위치
    std::vector<float> velX, velY, velZ;  // 속도
    std::vector<float> energy;     // 에너지
    std::vector<bool> isInteraction;  // 물리적 상호작용 발생 여부
    
    Step();
    virtual ~Step();  // 가상 소멸자 추가
    void Reset();
    
    // 스텝 추가
    void AddStep(int trackID, float time,
                float posX, float posY, float posZ,
                float velX, float velY, float velZ,
                float energy, bool interaction = false);
                
    ClassDef(Step, 1);
};

class Event : public TObject {
public:
    EventInfo eventInfo;
    ConfigParameters config;
    Track tracks;        // 실제 전자 정보
    Step steps;          // 궤적 포인트 정보
    
    Event();
    void Reset();
    
    // 애노드에 도달한 전자 수 계산
    int GetAnodeElectronCount() const;
    
    // 애노드에 도달한 전자 정보 얻기 (이전 버전 호환성)
    void GetAnodeElectrons(std::vector<int>& anodeTrackID,
                          std::vector<float>& finalPosX, 
                          std::vector<float>& finalPosY, 
                          std::vector<float>& finalPosZ,
                          std::vector<float>& finalEnergy,
                          std::vector<float>& transitTime,
                          std::vector<int>& processType,
                          std::vector<float>& Vx, 
                          std::vector<float>& Vy, 
                          std::vector<float>& Vz) const;
    
    ClassDef(Event, 2);  // 버전 증가
};

} // namespace mcp

#endif // MCP_EVENT_DATA_H 