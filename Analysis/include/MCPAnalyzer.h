#ifndef MCP_ANALYZER_H
#define MCP_ANALYZER_H

#include <string>
#include <vector>
#include "../../MCPSim/include/MCPEventData.h"

class MCPAnalyzer {
public:
    MCPAnalyzer();
    ~MCPAnalyzer();
    
    // Data loading - unified function
    bool Load(const std::string& filename);
    
    // Event access
    void SelectEvent(int eventIndex);
    const mcp::Event* GetEvent() const;
    const mcp::ConfigParameters* GetConfig() const;  // Access configuration parameters
    
    // Basic analysis functions - data extraction
    std::vector<float> GetEnergy(bool anodeOnly = false) const;  // Returns energy data
    std::vector<float> GetTime(bool anodeOnly = false) const;    // Returns time data
    std::vector<float> GetTransitTime(bool anodeOnly = false) const; // Returns transit time data
    
    // Velocity data (x,y,z components)
    std::vector<float> GetVx(bool anodeOnly = false) const;
    std::vector<float> GetVy(bool anodeOnly = false) const;
    std::vector<float> GetVz(bool anodeOnly = false) const;
    
    // Position data (x,y,z components)
    std::vector<float> GetPosX(bool anodeOnly = false) const;
    std::vector<float> GetPosY(bool anodeOnly = false) const;
    std::vector<float> GetPosZ(bool anodeOnly = false) const;
    
    // Electron count calculation (unified function)
    int GetEleCount(bool inAnode = false) const;
    
private:
    std::vector<mcp::Event*> events_;
    int currentEventIndex_;
    
    // Internal loading functions
    bool LoadRootFile(const std::string& filename);
};

#endif // MCP_ANALYZER_H 