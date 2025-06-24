#include "MCPEventData.h"
#include "MCPConfig.h"

ClassImp(mcp::EventInfo);
ClassImp(mcp::ConfigParameters);
ClassImp(mcp::Track);
ClassImp(mcp::Step);
ClassImp(mcp::Event);

namespace mcp {

EventInfo::EventInfo() {
    Reset();
}

void EventInfo::Reset() {
    initialEnergy = 0.f;
}

ConfigParameters::ConfigParameters() {
    Reset();
}

void ConfigParameters::Reset() {
    diff_pot = m = x0 = x1 = x2 = c_e = 0.f;
    pas = dia = alpha = limite = E0 = Resistance = q = 0.f;
    P_inf = P_1e = W = Ee = p = P_1r = r = Er = 0.f;
    e1 = e2 = r1 = r2 = t1 = t2 = t3 = t4 = 0.f;
    gama_ts0 = E_ts = s = 0.f;
}

void ConfigParameters::LoadFromConfig() {
    auto& config = MCPSim::Config::getInstance();
    
    // Physical constants
    diff_pot = config.get("diff_pot");
    m = config.get("m");
    x0 = config.get("x0");
    x1 = config.get("x1");
    x2 = config.get("x2");
    c_e = config.get("c_e");
    pas = config.get("pas");
    dia = config.get("dia");
    alpha = config.get("alpha");
    limite = config.get("limite");
    E0 = config.get("E0");
    Resistance = config.get("Resistance");
    q = config.get("q");
    
    // Secondary emission parameters
    P_inf = config.get("P_inf");
    P_1e = config.get("P_1e");
    W = config.get("W");
    Ee = config.get("Ee");
    p = config.get("p");
    P_1r = config.get("P_1r");
    r = config.get("r");
    Er = config.get("Er");
    e1 = config.get("e1");
    e2 = config.get("e2");
    r1 = config.get("r1");
    r2 = config.get("r2");
    t1 = config.get("t1");
    t2 = config.get("t2");
    t3 = config.get("t3");
    t4 = config.get("t4");
    gama_ts0 = config.get("gama_ts0");
    E_ts = config.get("E_ts");
    s = config.get("s");
}

Track::Track() {
    Reset();
}

Track::~Track() {
}

void Track::Reset() {
    nTracks = 0;
    trackID.clear();
    parentID.clear();
    birthTime.clear();
    birthPosX.clear();
    birthPosY.clear();
    birthPosZ.clear();
    birthVelX.clear();
    birthVelY.clear();
    birthVelZ.clear();
    birthEnergy.clear();
    processType.clear();
    isAnode.clear();
    finalTime.clear();
    finalPosX.clear();
    finalPosY.clear();
    finalPosZ.clear();
    finalVelX.clear();
    finalVelY.clear();
    finalVelZ.clear();
    finalEnergy.clear();
    transitTime.clear();
}

int Track::AddTrack(int parID, float time,
                  float posX, float posY, float posZ,
                  float velX, float velY, float velZ,
                  float energy, int procType) {
    int trID = nTracks + 1; // Start from 1
    nTracks++;
    
    trackID.push_back(trID);
    parentID.push_back(parID);
    birthTime.push_back(time);
    birthPosX.push_back(posX);
    birthPosY.push_back(posY);
    birthPosZ.push_back(posZ);
    birthVelX.push_back(velX);
    birthVelY.push_back(velY);
    birthVelZ.push_back(velZ);
    birthEnergy.push_back(energy);
    processType.push_back(procType);
    isAnode.push_back(0); // 초기값은 0 (애노드에 도달하지 않음)
    
    finalTime.push_back(time);
    finalPosX.push_back(posX);
    finalPosY.push_back(posY);
    finalPosZ.push_back(posZ);
    finalVelX.push_back(velX);
    finalVelY.push_back(velY);
    finalVelZ.push_back(velZ);
    finalEnergy.push_back(energy);
    transitTime.push_back(0.0f);
    
    return trID;
}

void Track::FinalizeTrack(int trID, int isAnodeHit, float time,
                        float posX, float posY, float posZ,
                        float velX, float velY, float velZ,
                        float energy) {
    int idx = FindTrackIndex(trID);
    if (idx >= 0) {
        isAnode[idx] = isAnodeHit;
        finalTime[idx] = time;
        finalPosX[idx] = posX;
        finalPosY[idx] = posY;
        finalPosZ[idx] = posZ;
        finalVelX[idx] = velX;
        finalVelY[idx] = velY;
        finalVelZ[idx] = velZ;
        finalEnergy[idx] = energy;
        transitTime[idx] = time - birthTime[idx]; 
    }
}

int Track::FindTrackIndex(int trID) const {
    for (int i = 0; i < nTracks; i++) {
        if (trackID[i] == trID) {
            return i;
        }
    }
    return -1;
}

Step::Step() {
    Reset();
}

Step::~Step() {
}

void Step::Reset() {
    nSteps = 0;
    trackID.clear();
    time.clear();
    posX.clear();
    posY.clear();
    posZ.clear();
    velX.clear();
    velY.clear();
    velZ.clear();
    energy.clear();
    isInteraction.clear();
}

void Step::AddStep(int trID, float t,
                 float pX, float pY, float pZ,
                 float vX, float vY, float vZ,
                 float e, bool interaction) {
    nSteps++;
    trackID.push_back(trID);
    time.push_back(t);
    posX.push_back(pX);
    posY.push_back(pY);
    posZ.push_back(pZ);
    velX.push_back(vX);
    velY.push_back(vY);
    velZ.push_back(vZ);
    energy.push_back(e);
    isInteraction.push_back(interaction);
}

Event::Event() {
    Reset();
}

void Event::Reset() {
    eventInfo.Reset();
    config.Reset();
    tracks.Reset();
    steps.Reset();
}

int Event::GetAnodeElectronCount() const {
    int count = 0;
    for (int i = 0; i < tracks.nTracks; i++) {
        if (tracks.isAnode[i] == 1) {
            count++;
        }
    }
    return count;
}

void Event::GetAnodeElectrons(std::vector<int>& anodeTrackID,
                            std::vector<float>& finalPosX, 
                            std::vector<float>& finalPosY, 
                            std::vector<float>& finalPosZ,
                            std::vector<float>& finalEnergy,
                            std::vector<float>& transitTime,
                            std::vector<int>& processType,
                            std::vector<float>& Vx, 
                            std::vector<float>& Vy, 
                            std::vector<float>& Vz) const {
    anodeTrackID.clear();
    finalPosX.clear();
    finalPosY.clear();
    finalPosZ.clear();
    finalEnergy.clear();
    transitTime.clear();
    processType.clear();
    Vx.clear();
    Vy.clear();
    Vz.clear();
    
    for (int i = 0; i < tracks.nTracks; i++) {
        if (tracks.isAnode[i] == 1) {
            anodeTrackID.push_back(tracks.trackID[i]);
            finalPosX.push_back(tracks.finalPosX[i]);
            finalPosY.push_back(tracks.finalPosY[i]);
            finalPosZ.push_back(tracks.finalPosZ[i]);
            finalEnergy.push_back(tracks.finalEnergy[i]);
            transitTime.push_back(tracks.transitTime[i]);
            processType.push_back(tracks.processType[i]);
            Vx.push_back(tracks.finalVelX[i]);
            Vy.push_back(tracks.finalVelY[i]);
            Vz.push_back(tracks.finalVelZ[i]);
        }
    }
}

} // namespace mcp 