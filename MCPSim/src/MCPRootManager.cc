#include "MCPRootManager.h"
#include <iostream>
#include <string>
#include <cstring>

MCPRootManager::MCPRootManager() 
    : fFile(nullptr), fEvtTree(nullptr), fEventObj(nullptr) {
    
    fEventObj = new mcp::Event();
}

MCPRootManager::~MCPRootManager() {
    Close();
    delete fEventObj;
}

bool MCPRootManager::OpenFile(const std::string& fname, const std::string& mode) {
    if (fFile) {
        Close();
    }
    
    fFile = new TFile(fname.c_str(), mode.c_str());
    if (!fFile || fFile->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file: " << fname << std::endl;
        delete fFile;
        fFile = nullptr;
        return false;
    }
    
    if (mode == "READ" || mode == "UPDATE") {
        fEvtTree = (TTree*)fFile->Get("Events");
        
        if (!fEvtTree) {
            std::cerr << "Error: Cannot find Events tree" << std::endl;
            return false;
        }
        
        // Track/Step 구조로 브랜치 연결
        fEvtTree->SetBranchAddress("EventInfo", &fEventObj->eventInfo);
        fEvtTree->SetBranchAddress("Config", &fEventObj->config);
        fEvtTree->SetBranchAddress("Tracks", &fEventObj->tracks);
        fEvtTree->SetBranchAddress("Steps", &fEventObj->steps);
    } else {
        CreateTrees();
    }
    
    return true;
}

void MCPRootManager::Close() {
    if (fFile) {
        fFile->Write();
        fFile->Close();
        delete fFile;
        fFile = nullptr;
    }
    fEvtTree = nullptr;
}

void MCPRootManager::CreateTrees() {
    fEvtTree = new TTree("Events", "MCP Simulation Events");
    
    fEvtTree->Branch("EventInfo", &fEventObj->eventInfo);
    fEvtTree->Branch("Config", &fEventObj->config);
    fEvtTree->Branch("Tracks", &fEventObj->tracks);
    fEvtTree->Branch("Steps", &fEventObj->steps);
}

void MCPRootManager::ClearBranches() {
    if (fEventObj) {
        fEventObj->Reset();
    }
}

void MCPRootManager::FillBranchesFromEvent(const mcp::Event& evt) {
    *fEventObj = evt;
}

void MCPRootManager::FillEventFromBranches(mcp::Event& evt) {
    evt.Reset();
    evt = *fEventObj;
}

bool MCPRootManager::WriteEvt(const mcp::Event& evt) {
    if (!fEvtTree) return false;
    
    ClearBranches();
    FillBranchesFromEvent(evt);
    fEvtTree->Fill();
    return true;
}

bool MCPRootManager::ReadEvt(mcp::Event& evt, int evtID) {
    if (!fEvtTree || evtID >= fEvtTree->GetEntries()) return false;
    
    fEvtTree->GetEntry(evtID);
    FillEventFromBranches(evt);
    
    return true;
}

int MCPRootManager::GetNumEvts() const {
    return fEvtTree ? fEvtTree->GetEntries() : 0;
} 