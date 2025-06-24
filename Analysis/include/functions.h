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

// 데이터 로드 함수
void loadRootFile(const std::string& filename, std::vector<mcp::Event*>& eventPtrs);

// 전자 증폭 트리 구조 생성 함수
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

// 시각화 함수들
// 2D 히스토그램: pore 내 위치별 전자 생성 밀도
void createDensityHistogram(const mcp::Event& event, const std::string& outputFileName);

// 단일 pore 내 증폭 과정 애니메이션 생성
void createCascadeAnimation(const mcp::Event& event, const std::string& outputFileName);

#endif // FUNCTIONS_H
