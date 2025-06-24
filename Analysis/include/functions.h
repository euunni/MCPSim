#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <vector>
#include <map>
#include <functional>
#include "../../MCPSim/include/MCPEventData.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TArrow.h"
#include "TText.h"
#include "TMath.h"
#include "TPolyLine3D.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TGaxis.h"
#include "TSystem.h"
#include "TPolyMarker3D.h"
#include "TKey.h"
#include "TCollection.h"
#include "TView3D.h"

// Data loading function
void loadRootFile(const std::string& filename, std::vector<mcp::Event*>& eventPtrs);

// Create electron amplification tree structure
struct ElectronNode {
    int trackID;
    int parentID;
    float posX, posY, posZ;
    float time;
    float energy;
    std::vector<int> childIDs;
};

typedef std::map<int, ElectronNode> ElectronCascadeTree;

ElectronCascadeTree buildElectronCascadeTree(const mcp::Event& event);

// Visualization functions
// 2D histogram: electron creation density per location in the pore
void createDensityHistogram(const mcp::Event& event, const std::string& outputFileName);

// Create animation of the amplification process in a single pore
void createCascadeAnimation(const mcp::Event& event, const std::string& outputFileName);

#endif // FUNCTIONS_H
