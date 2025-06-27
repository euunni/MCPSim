#include "../include/MCPAnalyzer.h"
#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"

MCPAnalyzer::MCPAnalyzer() : currentEventIndex_(0) {
}

MCPAnalyzer::~MCPAnalyzer() {
    // Clean up loaded events
    for (auto& event : events_) {
        delete event;
    }
    events_.clear();
}

bool MCPAnalyzer::Load(const std::string& filename) {
    // Only support ROOT files now
    if (filename.find(".root") == std::string::npos) {
        std::cerr << "Only ROOT files are supported" << std::endl;
        return false;
    }
    return LoadRootFile(filename);
}

bool MCPAnalyzer::LoadRootFile(const std::string& filename) {
    // Clean up any previously loaded events
    for (auto& event : events_) {
        delete event;
    }
    events_.clear();
    
    // Open ROOT file
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
    }
    
    // Get event tree
    TTree* tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "Cannot find tree 'Events' in file" << std::endl;
        file->Close();
        return false;
    }
    
    // Set branch addresses
    mcp::EventInfo* eventInfo = nullptr;
    mcp::ConfigParameters* config = nullptr;
    mcp::Track* tracks = nullptr;
    mcp::Step* steps = nullptr;
    
    tree->SetBranchAddress("EventInfo", &eventInfo);
    tree->SetBranchAddress("Config", &config);
    tree->SetBranchAddress("Tracks", &tracks);
    tree->SetBranchAddress("Steps", &steps);
    
    // Load events
    int nEntries = tree->GetEntries();
    events_.reserve(nEntries);
    
    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        mcp::Event* newEvent = new mcp::Event();
        
        if (eventInfo) newEvent->eventInfo = *eventInfo;
        if (config) newEvent->config = *config;
        if (tracks) newEvent->tracks = *tracks;
        if (steps) newEvent->steps = *steps;
        
        events_.push_back(newEvent);
    }
    
    file->Close();
    
    // Set current event to first event
    currentEventIndex_ = 0;
    
    std::cout << "Loaded " << events_.size() << " events from " << filename << std::endl;
    return true;
}

void MCPAnalyzer::SelectEvent(int eventIndex) {
    if (eventIndex >= 0 && eventIndex < events_.size()) {
        currentEventIndex_ = eventIndex;
        std::cout << "Selected event " << eventIndex << std::endl;
    } else {
        std::cerr << "Event index out of range: " << eventIndex << std::endl;
    }
}

const mcp::Event* MCPAnalyzer::GetEvent() const {
    if (events_.empty() || currentEventIndex_ >= events_.size()) {
        return nullptr;
    }
    return events_[currentEventIndex_];
}

const mcp::ConfigParameters* MCPAnalyzer::GetConfig() const {
    const mcp::Event* event = GetEvent();
    if (!event) return nullptr;
    return &(event->config);
}

std::vector<float> MCPAnalyzer::GetEnergy(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Get final energy of electrons (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            result.push_back(event->tracks.finalEnergy[i]);
        }
    }
    
    return result;
}

std::vector<float> MCPAnalyzer::GetTime(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Get final time of electrons (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            result.push_back(event->tracks.finalTime[i]);
        }
    }
    
    return result;
}

std::vector<float> MCPAnalyzer::GetTransitTime(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Calculate transit time as finalTime - birthTime (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            float transitTime = event->tracks.finalTime[i] - event->tracks.birthTime[i];
            result.push_back(transitTime);
        }
    }
    
    return result;
}

std::vector<float> MCPAnalyzer::GetVx(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Get final x-velocity of electrons (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            result.push_back(event->tracks.finalVelX[i]);
        }
    }
    
    return result;
}

std::vector<float> MCPAnalyzer::GetVy(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Get final y-velocity of electrons (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            result.push_back(event->tracks.finalVelY[i]);
        }
    }
    
    return result;
}

std::vector<float> MCPAnalyzer::GetVz(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Get final z-velocity of electrons (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            result.push_back(event->tracks.finalVelZ[i]);
        }
    }
    
    return result;
}

std::vector<float> MCPAnalyzer::GetPosX(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Get final x-position of electrons (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            result.push_back(event->tracks.finalPosX[i]);
        }
    }
    
    return result;
}

std::vector<float> MCPAnalyzer::GetPosY(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Get final y-position of electrons (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            result.push_back(event->tracks.finalPosY[i]);
        }
    }
    
    return result;
}

std::vector<float> MCPAnalyzer::GetPosZ(bool anodeOnly) const {
    std::vector<float> result;
    const mcp::Event* event = GetEvent();
    if (!event) return result;
    
    // Get final z-position of electrons (filtered by anode if requested)
    for (int i = 0; i < event->tracks.nTracks; i++) {
        if (!anodeOnly || event->tracks.isAnode[i] == 1) {
            result.push_back(event->tracks.finalPosZ[i]);
        }
    }
    
    return result;
}

int MCPAnalyzer::GetEleCount(bool inAnode) const {
    const mcp::Event* event = GetEvent();
    if (!event) return 0;
    
    if (inAnode) {
        // Count electrons that hit the anode
        int count = 0;
        for (int i = 0; i < event->tracks.nTracks; i++) {
            if (event->tracks.isAnode[i] == 1) {
                count++;
            }
        }
        return count;
    } else {
        // Return total number of electrons
        return event->tracks.nTracks;
    }
} 