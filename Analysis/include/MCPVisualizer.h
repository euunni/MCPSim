#ifndef MCP_VISUALIZER_H
#define MCP_VISUALIZER_H

#include <string>
#include "MCPAnalyzer.h"
#include "TCanvas.h"

class MCPVisualizer {
public:
    MCPVisualizer(const MCPAnalyzer* analyzer);
    ~MCPVisualizer();
    
    // MCP structure visualization
    void DrawPore(TCanvas* canvas);
    
    // Electron movement visualization - all return canvases for Analysis.cc to save
    TCanvas* DrawTrajectoryXY();   // Returns XY trajectory canvas
    TCanvas* DrawTrajectoryXZ();   // Returns XZ trajectory canvas
    TCanvas* DrawDensityXY();      // Returns XY density canvas
    TCanvas* DrawDensityXZ();      // Returns XZ density canvas
    TCanvas* AnimateCascadeFrame(int frameIndex, int totalFrames);  // Returns single animation frame
    int GetAnimationFrameCount();  // Returns total number of frames needed
    TCanvas* DrawMCP2DWithPoresAndSteps();
    // New: draw pores and only steps belonging to given trackID set (e.g., lineage)
    TCanvas* DrawMCP2DForTracks(const std::vector<int>& trackIDs);
    // Overlay multiple cascades; each inner vector is trackID list for one cascade
    TCanvas* DrawMCP2DOverlay(const std::vector<std::vector<int>>& trackSets);
    // Overlay multiple cascades (3D view); each inner vector is trackID list for one cascade
    TCanvas* DrawMCP3DOverlay(const std::vector<std::vector<int>>& trackSets);
    // Zoomed-in 3D view showing a few pores (±nRange rows) with cascades overlay
    TCanvas* DrawMCP3DZoom(const std::vector<std::vector<int>>& trackSets,
                           int nRange = 2, int nzRange = 2);
    
private:
    const MCPAnalyzer* analyzer_;
    
    // Internal utility functions
    void GetDataRange(float& minX, float& maxX, float& minY, float& maxY, float& minZ, float& maxZ) const;
    void DrawPoreBoundary(TVirtualPad* pad);
    Int_t GetColorByEnergy(float energy);
    
    // Helper functions for density and animation
    void CreateDensityHistograms(float minX, float maxX, float minY, float maxY, float minZ, float maxZ,
                               class TH2F*& hDensityXY, class TH2F*& hDensityXZ);
    void GetTimeRange(float& minTime, float& maxTime, float& minEnergy, float& maxEnergy);
};

#endif // MCP_VISUALIZER_H 