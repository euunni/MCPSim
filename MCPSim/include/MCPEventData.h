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

// Physical process type
enum ProcessType {
    PROCESS_BACKSCATTER = 0,  // Backscattering
    PROCESS_REDIFFUSED = 1,   // Rediffused electron
    PROCESS_SECONDARY = 2     // Secondary emission
};

// Track class: Store actual electron information
class Track : public TObject {
public:
    int nTracks;                   // Total number of tracks
    std::vector<int> trackID;      // Track ID
    std::vector<int> parentID;     // Parent track ID
    std::vector<float> birthTime;  // Birth time
    std::vector<float> birthPosX, birthPosY, birthPosZ;  // Birth position
    std::vector<float> birthVelX, birthVelY, birthVelZ;  // Initial velocity
    std::vector<float> birthEnergy;  // Initial energy
    std::vector<int> processType;  // Creation process type
    std::vector<int> isAnode;      // Anode hit (0: not hit, 1: hit)
    std::vector<float> finalTime;  // Final time (anode hit or end time)
    std::vector<float> finalPosX, finalPosY, finalPosZ;  // Final position
    std::vector<float> finalVelX, finalVelY, finalVelZ;  // Final velocity
    std::vector<float> finalEnergy;  // Final energy
    
    Track();
    virtual ~Track();  // Virtual destructor
    void Reset();
    
    // Add new track (electron creation)
    int AddTrack(int parentID, float time,
                float posX, float posY, float posZ,
                float velX, float velY, float velZ,
                float energy, int procType);
    
    // Finalize track (anode hit or end)
    void FinalizeTrack(int trackID, int isAnodeHit, float time,
                      float posX, float posY, float posZ,
                      float velX, float velY, float velZ,
                      float energy);
                      
    // Find track index by track ID
    int FindTrackIndex(int trackID) const;
    
    ClassDef(Track, 1);
};

// Step class: Store electron trajectory information
class Step : public TObject {
public:
    int nSteps;                    // Total number of steps
    std::vector<int> trackID;      // Which track's step
    std::vector<float> time;       // Time
    std::vector<float> posX, posY, posZ;  // Position
    std::vector<float> velX, velY, velZ;  // Velocity
    std::vector<float> energy;     // Energy
    
    Step();
    virtual ~Step();
    void Reset();
    
    // Add step
    void AddStep(int trackID, float time,
                float posX, float posY, float posZ,
                float velX, float velY, float velZ,
                float energy);
                
    ClassDef(Step, 1);
};

class Event : public TObject {
public:
    EventInfo eventInfo;
    ConfigParameters config;
    Track tracks;        // Actual electron information
    Step steps;          // Trajectory point information
    
    Event();
    void Reset();
    
    // Calculate number of electrons that hit the anode
    int GetAnodeElectronCount() const;
    
    // Get information about electrons that hit the anode (for backward compatibility)
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
    
    ClassDef(Event, 2);
};

} // namespace mcp

#endif // MCP_EVENT_DATA_H 