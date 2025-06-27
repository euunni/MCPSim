#ifndef MCP_ROOT_MANAGER_H
#define MCP_ROOT_MANAGER_H

#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "MCPEventData.h"

class MCPRootManager {
public:
    MCPRootManager();
    ~MCPRootManager();

    bool OpenFile(const std::string& fname, const std::string& mode = "RECREATE");
    void Close();

    bool WriteEvt(const mcp::Event& evt);
    bool ReadEvt(mcp::Event& evt, int evtID);
    
    int GetNumEvts() const;

private:
    void CreateTrees();
    void ClearBranches();
    void FillBranchesFromEvent(const mcp::Event& evt);
    void FillEventFromBranches(mcp::Event& evt);

    TFile* fFile;        
    TTree* fEvtTree;     
    
    mcp::Event* fEventObj;
};

#endif // MCP_ROOT_MANAGER_H 