#include <iostream>
#include <string>
#include <vector>
#include "include/functions.h"
#include "../MCPSim/include/MCPEventData.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " <root_file> [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --density     : Create electron density histogram" << std::endl;
    std::cout << "  --animation   : Create electron cascade animation" << std::endl;
    std::cout << "  --all         : Visualize all" << std::endl;
    std::cout << "  --event=N     : Analyze N-th event (default: 0)" << std::endl;
    std::cout << "  --output=NAME : Specify output file name (default: output)" << std::endl;
}

int main(int argc, char** argv) {

    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string inputFileName = argv[1];
    
    bool doDensity = false;
    bool doAnimation = false;
    int eventIndex = 0;
    std::string outputBaseName = "output";
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--density") {
            doDensity = true;
        } else if (arg == "--animation") {
            doAnimation = true;
        } else if (arg == "--all") {
            doDensity = doAnimation = true;
        } else if (arg.find("--event=") == 0) {
            eventIndex = std::stoi(arg.substr(8));
        } else if (arg.find("--output=") == 0) {
            outputBaseName = arg.substr(9);
        } else {
            std::cout << "Unknown option: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // Default: visualize all
    if (!(doDensity || doAnimation)) {
        doDensity = doAnimation = true;
    }
    
    // Load data
    std::vector<mcp::Event*> eventPtrs;
    std::string fullPath = inputFileName;
    if (fullPath.find(".root") == std::string::npos) {
        fullPath = "../output/" + fullPath + ".root";
    }
    
    TFile* file = TFile::Open(fullPath.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << fullPath << std::endl;
        return 1;
    }
    
    TTree* tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "Cannot find tree 'Events' in file: " << fullPath << std::endl;
        file->Close();
        return 1;
    }

    mcp::EventInfo* eventInfo = nullptr;
    mcp::ConfigParameters* config = nullptr;
    mcp::Track* tracks = nullptr;
    mcp::Step* steps = nullptr;
    
    tree->SetBranchAddress("EventInfo", &eventInfo);
    tree->SetBranchAddress("Config", &config);
    tree->SetBranchAddress("Tracks", &tracks);
    tree->SetBranchAddress("Steps", &steps);

    int nEntries = tree->GetEntries();
    eventPtrs.reserve(nEntries);
    
    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        mcp::Event* newEvent = new mcp::Event();
        
        if (eventInfo) newEvent->eventInfo = *eventInfo;
        if (config) newEvent->config = *config;
        if (tracks) newEvent->tracks = *tracks;
        if (steps) newEvent->steps = *steps;
        
        eventPtrs.push_back(newEvent);
    }
    
    file->Close();    

    if (eventPtrs.empty()) {
        std::cerr << "No events loaded from file." << std::endl;
        return 1;
    }

    if (eventIndex >= eventPtrs.size()) {
        std::cerr << "Event index " << eventIndex << " is out of range. Max index: " << (eventPtrs.size() - 1) << std::endl;
        return 1;
    }

    const mcp::Event& selectedEvent = *(eventPtrs[eventIndex]);
    
    std::cout << "Total number of electrons: " << selectedEvent.tracks.nTracks << std::endl;
    std::cout << "Number of electrons reaching the anode: " << selectedEvent.GetAnodeElectronCount() << std::endl;
    std::cout << "Number of trajectory steps: " << selectedEvent.steps.nSteps << std::endl;
    
    // Create electron cascade tree
    ElectronCascadeTree cascadeTree = buildElectronCascadeTree(selectedEvent);
    std::cout << "Electron cascade tree created. Total number of nodes: " << cascadeTree.size() << std::endl;
    
    std::string outputPrefix = outputBaseName + "_event" + std::to_string(eventIndex);
    
    // Visualize
    if (doDensity) {
        std::cout << "Creating electron density histogram..." << std::endl;
        createDensityHistogram(selectedEvent, outputPrefix + "_density.pdf");
    }

    if (doAnimation) {
        std::cout << "Creating electron amplification animation..." << std::endl;
        createCascadeAnimation(selectedEvent, outputPrefix + "_animation");
    }
    
    for (auto ptr : eventPtrs) {
        delete ptr;
    }
    
    std::cout << "All visualization tasks completed." << std::endl;
    return 0;
}
