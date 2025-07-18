#include "../include/MCPVisualizer.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TArrow.h"
#include "TText.h"
#include "TMath.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TGaxis.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TView3D.h"
#include "TGLViewer.h"   // for camera control in OpenGL viewer
#include "TSystem.h"
#include "TDirectory.h"
#include "TBox.h"
#include <unordered_set> // Added for DrawMCP2DForTracks
#include "TGeoManager.h"

MCPVisualizer::MCPVisualizer(const MCPAnalyzer* analyzer) : analyzer_(analyzer) {
    // Set ROOT style
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    // Force OpenGL viewer for 3-D pads so that interactive rotation/zoom is available
    gStyle->SetCanvasPreferGL(kTRUE);
}

MCPVisualizer::~MCPVisualizer() {
}

void MCPVisualizer::GetDataRange(float& minX, float& maxX, float& minY, float& maxY, float& minZ, float& maxZ) const {
    minX = minY = minZ = std::numeric_limits<float>::max();
    maxX = maxY = maxZ = -std::numeric_limits<float>::max();
    
    const mcp::Event* event = analyzer_->GetEvent();
    if (!event) return;
    
    // Get range from track data
    for (int i = 0; i < event->tracks.nTracks; i++) {
        // Birth position
        minX = std::min(minX, event->tracks.birthPosX[i]);
        maxX = std::max(maxX, event->tracks.birthPosX[i]);
        minY = std::min(minY, event->tracks.birthPosY[i]);
        maxY = std::max(maxY, event->tracks.birthPosY[i]);
        minZ = std::min(minZ, event->tracks.birthPosZ[i]);
        maxZ = std::max(maxZ, event->tracks.birthPosZ[i]);
        
        // Final position
        minX = std::min(minX, event->tracks.finalPosX[i]);
        maxX = std::max(maxX, event->tracks.finalPosX[i]);
        minY = std::min(minY, event->tracks.finalPosY[i]);
        maxY = std::max(maxY, event->tracks.finalPosY[i]);
        minZ = std::min(minZ, event->tracks.finalPosZ[i]);
        maxZ = std::max(maxZ, event->tracks.finalPosZ[i]);
    }
    
    // Get range from step data if available
    for (int i = 0; i < event->steps.nSteps; i++) {
        minX = std::min(minX, event->steps.posX[i]);
        maxX = std::max(maxX, event->steps.posX[i]);
        minY = std::min(minY, event->steps.posY[i]);
        maxY = std::max(maxY, event->steps.posY[i]);
        minZ = std::min(minZ, event->steps.posZ[i]);
        maxZ = std::max(maxZ, event->steps.posZ[i]);
    }
    
    // If the range is too narrow, expand it
    if (maxX - minX < 1.0) {
        float mid = (minX + maxX) / 2;
        minX = mid - 0.5;
        maxX = mid + 0.5;
    }
    if (maxY - minY < 1.0) {
        float mid = (minY + maxY) / 2;
        minY = mid - 0.5;
        maxY = mid + 0.5;
    }
    if (maxZ - minZ < 1.0) {
        float mid = (minZ + maxZ) / 2;
        minZ = mid - 0.5;
        maxZ = mid + 0.5;
    }
    
    // Add padding exactly as in original
    float padX = (maxX - minX) * 0.1;
    float padY = (maxY - minY) * 0.1;
    float padZ = (maxZ - minZ) * 0.1;
    
    minX -= padX;
    maxX += padX;
    minY -= padY;
    maxY += padY;
    minZ -= padZ;
    maxZ += padZ;
}

Int_t MCPVisualizer::GetColorByEnergy(float energy) {
    // Get energy range
    std::vector<float> energies = analyzer_->GetEnergy();
    if (energies.empty()) return kBlack;
    
    float minEnergy = *std::min_element(energies.begin(), energies.end());
    float maxEnergy = *std::max_element(energies.begin(), energies.end());
    
    // Normalize energy to 0~1
    float normalizedEnergy = (energy - minEnergy) / (maxEnergy - minEnergy);
    if (normalizedEnergy < 0) normalizedEnergy = 0;
    if (normalizedEnergy > 1) normalizedEnergy = 1;
    
    gStyle->SetPalette(kViridis);
    
    // Get color from palette
    Int_t colorIndex = static_cast<Int_t>(normalizedEnergy * 255);
    Int_t color = TColor::GetColorPalette(colorIndex);
    
    return color;
}

void MCPVisualizer::DrawPoreBoundary(TVirtualPad* pad) {
    pad->cd();
    
    // Get data range
    float minX, maxX, minY, maxY, minZ, maxZ;
    GetDataRange(minX, maxX, minY, maxY, minZ, maxZ);
    
    // Get configuration parameters from ROOT file
    const mcp::ConfigParameters* config = analyzer_->GetConfig();
    if (!config) {
        std::cerr << "No configuration data available" << std::endl;
        return;
    }
    
    // Pore parameters (from ROOT file Config branch)
    const float poreStartY = 6.0;  // Y center coordinate at the start point (keep default)
    const float poreDiameter = config->dia;  // dia
    const float poreAngle = 0.13;   // alpha
    const float poreRadius = poreDiameter / 2.0;  // Pore radius
    const float poreStartX = config->x0;  // x0
    const float poreEndX = config->x1;    // x1
    
    // Draw top and bottom lines of the pore
    TLine* topLine = new TLine();
    topLine->SetLineColor(kBlack);
    topLine->SetLineWidth(1);
    
    TLine* bottomLine = new TLine();
    bottomLine->SetLineColor(kBlack);
    bottomLine->SetLineWidth(1);
    
    // Calculate start and end points (apply slope)
    float startTopY = poreStartY + poreRadius;
    float startBottomY = poreStartY - poreRadius;
    
    float endTopY = poreStartY + poreRadius + (poreEndX - poreStartX) * tan(poreAngle);
    float endBottomY = poreStartY - poreRadius + (poreEndX - poreStartX) * tan(poreAngle);
    
    topLine->DrawLine(poreStartX, startTopY, poreEndX, endTopY);
    bottomLine->DrawLine(poreStartX, startBottomY, poreEndX, endBottomY);

    // Draw start and end lines of the pore
    TLine* startLine = new TLine(poreStartX, startBottomY, poreStartX, startTopY);
    startLine->SetLineColor(kBlack);
    startLine->SetLineWidth(1);
    startLine->SetLineStyle(2); 
    startLine->Draw();
    
    TLine* endLine = new TLine(poreEndX, endBottomY, poreEndX, endTopY);
    endLine->SetLineColor(kBlack);
    endLine->SetLineWidth(1);
    endLine->SetLineStyle(2);
    endLine->Draw();
}

void MCPVisualizer::DrawPore(TCanvas* canvas) {
    if (!canvas) return;
    canvas->cd();
    
    // Clear canvas
    canvas->Clear();
    
    // Draw pore boundary
    DrawPoreBoundary(canvas);
}

TCanvas* MCPVisualizer::DrawTrajectoryXY() {
    // Get event data
    const mcp::Event* event = analyzer_->GetEvent();
    if (!event) {
        std::cerr << "No event data available" << std::endl;
        return nullptr;
    }
    
    std::cout << "Creating XY trajectory visualization..." << std::endl;
    
    // Get data range
    float minX, maxX, minY, maxY, minZ, maxZ;
    GetDataRange(minX, maxX, minY, maxY, minZ, maxZ);
    
    // Create canvas with unique name
    static int trajCounter = 0;
    trajCounter++;
    
    TCanvas* c = new TCanvas(Form("c_traj_xy_%d", trajCounter), "Electron Trajectories (XY)", 800, 600);
    c->SetRightMargin(0.15);  // Make room for color bar
    
    // XY projection frame
    TH2F* frameXY = new TH2F(Form("frameXY_%d", trajCounter), "Electron Trajectories (XY);X [#mum];Y [#mum]",
                           100, minX, maxX, 100, minY, maxY);
    frameXY->SetStats(0);
    frameXY->Draw();
    
    // Draw pore boundary
    DrawPoreBoundary(c);
    
    // Limit the number of points to draw for performance
    const int maxPointsToDraw = 10000;
    int stepsToDraw = event->steps.nSteps;
    int stepInterval = 1;
    
    // If too many steps, use sampling
    if (stepsToDraw > maxPointsToDraw) {
        stepInterval = stepsToDraw / maxPointsToDraw + 1;
    }
    
    // Get energy range for color mapping
    std::vector<float> energies = analyzer_->GetEnergy();
    float minEnergy = *std::min_element(energies.begin(), energies.end());
    float maxEnergy = *std::max_element(energies.begin(), energies.end());
    
    // Draw trajectories (XY) with markers
    for (int i = 0; i < event->steps.nSteps; i += stepInterval) {
        float x = event->steps.posX[i];
        float y = event->steps.posY[i];
        float energy = event->steps.energy[i];
        
        TMarker* marker = new TMarker(x, y, 20);
        marker->SetMarkerSize(0.5);
        marker->SetMarkerColor(GetColorByEnergy(energy));
        marker->Draw();
    }
    
    // Create color bar manually using TGaxis
    c->Update();
    
    // Get canvas coordinates for color bar
    Double_t x1 = 0.85;  // Left edge of color bar (in NDC coordinates)
    Double_t x2 = 0.90;  // Right edge of color bar
    Double_t y1 = 0.1;   // Bottom edge
    Double_t y2 = 0.9;   // Top edge
    
    // Convert NDC to user coordinates
    Double_t ux1, uy1, ux2, uy2;
    c->GetRange(ux1, uy1, ux2, uy2);
    Double_t px1 = ux1 + x1 * (ux2 - ux1);
    Double_t px2 = ux1 + x2 * (ux2 - ux1);
    Double_t py1 = uy1 + y1 * (uy2 - uy1);
    Double_t py2 = uy1 + y2 * (uy2 - uy1);
    
    // Draw colored rectangles for color bar
    Int_t nColorBins = 50;
    Double_t binHeight = (py2 - py1) / nColorBins;
    
    for (Int_t i = 0; i < nColorBins; i++) {
        Float_t energy = minEnergy + (maxEnergy - minEnergy) * i / nColorBins;
        Int_t color = GetColorByEnergy(energy);
        
        TBox* box = new TBox(px1, py1 + i * binHeight, px2, py1 + (i + 1) * binHeight);
        box->SetFillColor(color);
        box->SetLineColor(color);
        box->Draw();
    }
    
    // Add axis labels
    TGaxis* colorAxis = new TGaxis(px2, py1, px2, py2, minEnergy, maxEnergy, 510, "+L");
    colorAxis->SetTitle("Energy [eV]");
    colorAxis->SetTitleOffset(1.2);
    colorAxis->SetLabelSize(0.03);
    colorAxis->SetTitleSize(0.03);
    colorAxis->Draw();
    
    // Add information text
    TPaveText* info = new TPaveText(0.15, 0.85, 0.55, 0.90, "NDC");
    info->SetFillColor(0);
    info->SetTextAlign(12);
    info->SetTextSize(0.035);
    info->AddText(Form("Total electrons: %d, Anode hits: %d", 
                     event->tracks.nTracks, event->GetAnodeElectronCount()));
    info->Draw();
    
    std::cout << "XY trajectory visualization created." << std::endl;
    
    return c;
}

TCanvas* MCPVisualizer::DrawTrajectoryXZ() {
    // Get event data
    const mcp::Event* event = analyzer_->GetEvent();
    if (!event) {
        std::cerr << "No event data available" << std::endl;
        return nullptr;
    }
    
    std::cout << "Creating XZ trajectory visualization..." << std::endl;
    
    // Get data range
    float minX, maxX, minY, maxY, minZ, maxZ;
    GetDataRange(minX, maxX, minY, maxY, minZ, maxZ);
    
    // Create canvas with unique name
    static int trajCounter = 0;
    trajCounter++;
    
    TCanvas* c = new TCanvas(Form("c_traj_xz_%d", trajCounter), "Electron Trajectories (XZ)", 800, 600);
    c->SetRightMargin(0.15);  // Make room for color bar
    
    // XZ projection frame
    TH2F* frameXZ = new TH2F(Form("frameXZ_%d", trajCounter), "Electron Trajectories (XZ);X [#mum];Z [#mum]",
                           100, minX, maxX, 100, minZ, maxZ);
    frameXZ->SetStats(0);
    frameXZ->Draw();
    
    // Limit the number of points to draw for performance
    const int maxPointsToDraw = 10000;
    int stepsToDraw = event->steps.nSteps;
    int stepInterval = 1;
    
    // If too many steps, use sampling
    if (stepsToDraw > maxPointsToDraw) {
        stepInterval = stepsToDraw / maxPointsToDraw + 1;
    }
    
    // Get energy range for color mapping
    std::vector<float> energies = analyzer_->GetEnergy();
    float minEnergy = *std::min_element(energies.begin(), energies.end());
    float maxEnergy = *std::max_element(energies.begin(), energies.end());
    
    // Draw trajectories (XZ) with markers
    for (int i = 0; i < event->steps.nSteps; i += stepInterval) {
        float x = event->steps.posX[i];
        float z = event->steps.posZ[i];
        float energy = event->steps.energy[i];
        
        TMarker* marker = new TMarker(x, z, 20);
        marker->SetMarkerSize(0.5);
        marker->SetMarkerColor(GetColorByEnergy(energy));
        marker->Draw();
    }
    
    // Create color bar manually using TGaxis
    c->Update();
    
    // Get canvas coordinates for color bar
    Double_t x1 = 0.85;  // Left edge of color bar (in NDC coordinates)
    Double_t x2 = 0.90;  // Right edge of color bar
    Double_t y1 = 0.1;   // Bottom edge
    Double_t y2 = 0.9;   // Top edge
    
    // Convert NDC to user coordinates
    Double_t ux1, uy1, ux2, uy2;
    c->GetRange(ux1, uy1, ux2, uy2);
    Double_t px1 = ux1 + x1 * (ux2 - ux1);
    Double_t px2 = ux1 + x2 * (ux2 - ux1);
    Double_t py1 = uy1 + y1 * (uy2 - uy1);
    Double_t py2 = uy1 + y2 * (uy2 - uy1);
    
    // Draw colored rectangles for color bar
    Int_t nColorBins = 50;
    Double_t binHeight = (py2 - py1) / nColorBins;
    
    for (Int_t i = 0; i < nColorBins; i++) {
        Float_t energy = minEnergy + (maxEnergy - minEnergy) * i / nColorBins;
        Int_t color = GetColorByEnergy(energy);
        
        TBox* box = new TBox(px1, py1 + i * binHeight, px2, py1 + (i + 1) * binHeight);
        box->SetFillColor(color);
        box->SetLineColor(color);
        box->Draw();
    }
    
    // Add axis labels
    TGaxis* colorAxis = new TGaxis(px2, py1, px2, py2, minEnergy, maxEnergy, 510, "+L");
    colorAxis->SetTitle("Energy [eV]");
    colorAxis->SetTitleOffset(1.2);
    colorAxis->SetLabelSize(0.03);
    colorAxis->SetTitleSize(0.03);
    colorAxis->Draw();
    
    // Add information text
    TPaveText* info = new TPaveText(0.15, 0.85, 0.55, 0.90, "NDC");
    info->SetFillColor(0);
    info->SetTextAlign(12);
    info->SetTextSize(0.035);
    info->AddText(Form("Total electrons: %d, Anode hits: %d", 
                     event->tracks.nTracks, event->GetAnodeElectronCount()));
    info->Draw();
    
    std::cout << "XZ trajectory visualization created." << std::endl;
    
    return c;
}

void MCPVisualizer::GetTimeRange(float& minTime, float& maxTime, float& minEnergy, float& maxEnergy) {
    const mcp::Event* event = analyzer_->GetEvent();
    if (!event) return;
    
    minTime = minEnergy = std::numeric_limits<float>::max();
    maxTime = maxEnergy = -std::numeric_limits<float>::max();
    
    // Use track information to calculate time range
    for (int i = 0; i < event->tracks.nTracks; i++) {
        minTime = std::min(minTime, event->tracks.birthTime[i]);
        maxTime = std::max(maxTime, event->tracks.finalTime[i]);
        minEnergy = std::min(minEnergy, event->tracks.birthEnergy[i]);
        maxEnergy = std::max(maxEnergy, event->tracks.finalEnergy[i]);
    }
    
    // Use step information to calculate time range
    for (int i = 0; i < event->steps.nSteps; i++) {
        minTime = std::min(minTime, event->steps.time[i]);
        maxTime = std::max(maxTime, event->steps.time[i]);
        minEnergy = std::min(minEnergy, event->steps.energy[i]);
        maxEnergy = std::max(maxEnergy, event->steps.energy[i]);
    }
}

void MCPVisualizer::CreateDensityHistograms(float minX, float maxX, float minY, float maxY, float minZ, float maxZ,
                                           TH2F*& hDensityXY, TH2F*& hDensityXZ) {
    const mcp::Event* event = analyzer_->GetEvent();
    if (!event) return;
    
    // Use timestamp to create unique names
    static int histCounter = 0;
    histCounter++;
    
    hDensityXY = new TH2F(Form("hDensityXY_%d", histCounter), "XY Plane Electron Density",
                         200, minX, maxX, 200, -200, 300);  
    hDensityXY->SetStats(0);  
    hDensityXY->GetXaxis()->SetTitle("X (#mum)");
    hDensityXY->GetYaxis()->SetTitle("Y (#mum)");
    hDensityXY->GetXaxis()->SetTitleOffset(1.2);  
    hDensityXY->GetYaxis()->SetTitleOffset(1.5); 
    hDensityXY->GetXaxis()->SetLabelSize(0.04);  
    hDensityXY->GetYaxis()->SetLabelSize(0.04);  
    hDensityXY->GetXaxis()->SetTitleSize(0.04);  
    hDensityXY->GetYaxis()->SetTitleSize(0.04);  
    
    hDensityXZ = new TH2F(Form("hDensityXZ_%d", histCounter), "XZ Plane Electron Density",
                         200, minX, maxX, 200, minZ, maxZ);
    hDensityXZ->SetStats(0);    
    hDensityXZ->GetXaxis()->SetTitle("X (#mum)");
    hDensityXZ->GetYaxis()->SetTitle("Z (#mum)");
    hDensityXZ->GetXaxis()->SetTitleOffset(1.2);  
    hDensityXZ->GetYaxis()->SetTitleOffset(1.5); 
    hDensityXZ->GetXaxis()->SetLabelSize(0.04);  
    hDensityXZ->GetYaxis()->SetLabelSize(0.04);  
    hDensityXZ->GetXaxis()->SetTitleSize(0.04);  
    hDensityXZ->GetYaxis()->SetTitleSize(0.04);  
    
    // Fill histograms with track and step data
    for (int i = 0; i < event->tracks.nTracks; i++) {
        // Start point
        float x = event->tracks.birthPosX[i];
        float y = event->tracks.birthPosY[i];
        float z = event->tracks.birthPosZ[i];
        
        if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ) {
            hDensityXY->Fill(x, y);
            hDensityXZ->Fill(x, z);
        }
        
        // End point
        x = event->tracks.finalPosX[i];
        y = event->tracks.finalPosY[i];
        z = event->tracks.finalPosZ[i];
        
        if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ) {
            hDensityXY->Fill(x, y);
            hDensityXZ->Fill(x, z);
        }
    }
    
    for (int i = 0; i < event->steps.nSteps; i++) {
        float x = event->steps.posX[i];
        float y = event->steps.posY[i];
        float z = event->steps.posZ[i];
        
        if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ) {
            hDensityXY->Fill(x, y);
            hDensityXZ->Fill(x, z);
        }
    }
}

TCanvas* MCPVisualizer::DrawDensityXY() {
    const mcp::Event* event = analyzer_->GetEvent();
    if (!event) {
        std::cerr << "No event data available" << std::endl;
        return nullptr;
    }
    
    // Get data range
    float minX, maxX, minY, maxY, minZ, maxZ;
    GetDataRange(minX, maxX, minY, maxY, minZ, maxZ);
    
    // Create histograms
    TH2F* hDensityXY = nullptr;
    TH2F* hDensityXZ = nullptr;
    CreateDensityHistograms(minX, maxX, minY, maxY, minZ, maxZ, hDensityXY, hDensityXZ);
    
    // Create canvas
    TCanvas* canvas = new TCanvas("canvasDensityXY", "XY Plane Electron Density", 1000, 800);
    canvas->SetLeftMargin(0.15);  
    canvas->SetRightMargin(0.15);  
    canvas->SetBottomMargin(0.15); 
    canvas->SetTopMargin(0.1);     

    gStyle->SetPalette(kViridis);
    
    hDensityXY->Draw("COLZ");
    
    // Draw pore boundary
    DrawPoreBoundary(canvas);
    
    // Get time range
    float minTime, maxTime, minEnergy, maxEnergy;
    GetTimeRange(minTime, maxTime, minEnergy, maxEnergy);
    
    // Add time range information
    TPaveText* timeInfo = new TPaveText(0.15, 0.85, 0.55, 0.90, "NDC");
    timeInfo->SetFillColor(0);
    timeInfo->SetTextAlign(12);
    timeInfo->SetTextSize(0.035);
    timeInfo->AddText(Form("Time range: %.2f - %.2f ps", minTime, maxTime));
    timeInfo->Draw();
    
    // Add electron count information
    TPaveText* electronInfo = new TPaveText(0.15, 0.80, 0.55, 0.85, "NDC");
    electronInfo->SetFillColor(0);
    electronInfo->SetTextAlign(12);
    electronInfo->SetTextSize(0.035);
    electronInfo->AddText(Form("Total electrons: %d, Anode hits: %d", 
                             event->tracks.nTracks, event->GetAnodeElectronCount()));
    electronInfo->Draw();
    
    // Clean up unused histogram
    delete hDensityXZ;
    
    return canvas;
}

TCanvas* MCPVisualizer::DrawDensityXZ() {
    const mcp::Event* event = analyzer_->GetEvent();
    if (!event) {
        std::cerr << "No event data available" << std::endl;
        return nullptr;
    }
    
    // Get data range
    float minX, maxX, minY, maxY, minZ, maxZ;
    GetDataRange(minX, maxX, minY, maxY, minZ, maxZ);
    
    // Create histograms
    TH2F* hDensityXY = nullptr;
    TH2F* hDensityXZ = nullptr;
    CreateDensityHistograms(minX, maxX, minY, maxY, minZ, maxZ, hDensityXY, hDensityXZ);
    
    // Create canvas
    TCanvas* canvas = new TCanvas("canvasDensityXZ", "XZ Plane Electron Density", 1000, 800);
    canvas->SetLeftMargin(0.15);  
    canvas->SetRightMargin(0.15);  
    canvas->SetBottomMargin(0.15); 
    canvas->SetTopMargin(0.1);     

    gStyle->SetPalette(kViridis);
    
    hDensityXZ->Draw("COLZ");
    
    // Get time range
    float minTime, maxTime, minEnergy, maxEnergy;
    GetTimeRange(minTime, maxTime, minEnergy, maxEnergy);
    
    // Add time range information
    TPaveText* timeInfo = new TPaveText(0.15, 0.85, 0.55, 0.90, "NDC");
    timeInfo->SetFillColor(0);
    timeInfo->SetTextAlign(12);
    timeInfo->SetTextSize(0.035);
    timeInfo->AddText(Form("Time range: %.2f - %.2f ps", minTime, maxTime));
    timeInfo->Draw();
    
    // Add electron count information
    TPaveText* electronInfo = new TPaveText(0.15, 0.80, 0.55, 0.85, "NDC");
    electronInfo->SetFillColor(0);
    electronInfo->SetTextAlign(12);
    electronInfo->SetTextSize(0.035);
    electronInfo->AddText(Form("Total electrons: %d, Anode hits: %d", 
                             event->tracks.nTracks, event->GetAnodeElectronCount()));
    electronInfo->Draw();
    
    // Clean up unused histogram
    delete hDensityXY;
    
    return canvas;
}

int MCPVisualizer::GetAnimationFrameCount() {
    return 15;  // Same as original
}

TCanvas* MCPVisualizer::AnimateCascadeFrame(int frameIndex, int totalFrames) {
    const mcp::Event* event = analyzer_->GetEvent();
    if (!event) {
        std::cerr << "No event data available" << std::endl;
        return nullptr;
    }
    
    // Get data range
    float minX, maxX, minY, maxY, minZ, maxZ;
    GetDataRange(minX, maxX, minY, maxY, minZ, maxZ);
    
    // Get time and energy range
    float minTime, maxTime, minEnergy, maxEnergy;
    GetTimeRange(minTime, maxTime, minEnergy, maxEnergy);
    
    // Calculate current time for this frame
    float timeStep = (maxTime - minTime) / totalFrames;
    float currentTime = minTime + frameIndex * timeStep;
    float timeWindow = 0.1f;
    
    gStyle->SetPalette(kViridis);
    
    TCanvas* canvas = new TCanvas(Form("canvasAnim_%d", frameIndex), 
                                Form("MCP Pore Electron Cascade (Time: %.2f ps)", currentTime), 
                                1300, 1100);  
    canvas->SetLeftMargin(0.15);  
    canvas->SetRightMargin(0.05);  
    canvas->SetBottomMargin(0.15); 
    canvas->SetTopMargin(0.1);     
    canvas->Divide(2, 2, 0.01, 0.01);  
    
    // XY plane - use unique names for each frame
    canvas->cd(1);
    gPad->SetLeftMargin(0.15);  
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.1);     
    
    TH2F* hXY = new TH2F(Form("hXY_%d", frameIndex), Form("XY View (Time: %.2f ps)", currentTime),
                        100, minX, maxX, 100, -600, 600);
    hXY->SetStats(0);  
    hXY->GetXaxis()->SetTitle("X (#mum)");
    hXY->GetYaxis()->SetTitle("Y (#mum)");
    hXY->GetXaxis()->SetTitleOffset(1.2);  
    hXY->GetYaxis()->SetTitleOffset(1.5);  
    hXY->GetXaxis()->SetLabelSize(0.04);  
    hXY->GetYaxis()->SetLabelSize(0.04);  
    hXY->GetXaxis()->SetTitleSize(0.04);  
    hXY->GetYaxis()->SetTitleSize(0.04);  
    hXY->Draw();
    
    // Draw pore boundary (XY plane)
    DrawPoreBoundary(canvas->cd(1));
    
    // XZ plane - use unique names for each frame
    canvas->cd(2);
    gPad->SetLeftMargin(0.15);  
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.1);     
    
    TH2F* hXZ = new TH2F(Form("hXZ_%d", frameIndex), "XZ View",
                        100, minX, maxX, 100, -600, 600);
    hXZ->SetStats(0);  
    hXZ->GetXaxis()->SetTitle("X (#mum)");
    hXZ->GetYaxis()->SetTitle("Z (#mum)");
    hXZ->GetXaxis()->SetTitleOffset(1.2);  
    hXZ->GetYaxis()->SetTitleOffset(1.5);  
    hXZ->GetXaxis()->SetLabelSize(0.04);  
    hXZ->GetYaxis()->SetLabelSize(0.04);  
    hXZ->GetXaxis()->SetTitleSize(0.04);  
    hXZ->GetYaxis()->SetTitleSize(0.04);  
    hXZ->Draw();
    
    // ZY plane (axis position exchange) - use unique names for each frame
    canvas->cd(3);
    gPad->SetLeftMargin(0.15);  
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.1);     
    
    TH2F* hYZ = new TH2F(Form("hYZ_%d", frameIndex), "ZY View",  
                        100, -600, 600, 100, -600, 600);  
    hYZ->SetStats(0);  
    hYZ->GetXaxis()->SetTitle("Z (#mum)");  
    hYZ->GetYaxis()->SetTitle("Y (#mum)");
    hYZ->GetXaxis()->SetTitleOffset(1.2);  
    hYZ->GetYaxis()->SetTitleOffset(1.5);  
    hYZ->GetXaxis()->SetLabelSize(0.04);  
    hYZ->GetYaxis()->SetLabelSize(0.04);  
    hYZ->GetXaxis()->SetTitleSize(0.04);  
    hYZ->GetYaxis()->SetTitleSize(0.04);  
    hYZ->Draw();
    
    // 3D view - use unique names for each frame
    canvas->cd(4);
    gPad->SetLeftMargin(0.15);  
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.18);  
    gPad->SetTopMargin(0.1);     
    gPad->SetTheta(30);
    gPad->SetPhi(30);
    
    TH3F* h3D = new TH3F(Form("h3D_%d", frameIndex), "3D View",
                        10, minX, maxX, 10, -600, 600, 10, -600, 600);
    h3D->SetStats(0);  
    h3D->GetXaxis()->SetTitle("X (#mum)");
    h3D->GetYaxis()->SetTitle("Z (#mum)");
    h3D->GetZaxis()->SetTitle("Y (#mum)");
    h3D->GetXaxis()->SetTitleOffset(2.0);  
    h3D->GetYaxis()->SetTitleOffset(2.2);  
    h3D->GetZaxis()->SetTitleOffset(2.8);  
    h3D->SetTitleSize(0.04);  
    
    h3D->Draw("BOX");
    
    TView* view = gPad->GetView();
    if (view) {
        view->RotateView(35, 45);  
    }
    
    int pointsDrawn = 0;
    
    // Show tracks in the current time window
    for (int i = 0; i < event->tracks.nTracks; i++) {
        float birthTime = event->tracks.birthTime[i];
        float finalTime = event->tracks.finalTime[i];
        
        // Show only tracks in the current time window
        if ((birthTime <= currentTime + timeWindow && finalTime >= currentTime - timeWindow) ||
            (birthTime >= currentTime - timeWindow && birthTime <= currentTime + timeWindow) ||
            (finalTime >= currentTime - timeWindow && finalTime <= currentTime + timeWindow)) {
            
            // Select the point closer to the current time among the start and end points
            float time, x, y, z, energy;
            if (std::abs(birthTime - currentTime) < std::abs(finalTime - currentTime)) {
                // The start point is closer
                time = event->tracks.birthTime[i];
                x = event->tracks.birthPosX[i];
                y = event->tracks.birthPosY[i];
                z = event->tracks.birthPosZ[i];
                energy = event->tracks.birthEnergy[i];
            } else {
                // The end point is closer
                time = event->tracks.finalTime[i];
                x = event->tracks.finalPosX[i];
                y = event->tracks.finalPosY[i];
                z = event->tracks.finalPosZ[i];
                energy = event->tracks.finalEnergy[i];
            }
            
            Int_t color = GetColorByEnergy(energy);
            
            canvas->cd(1);
            TMarker* mXY = new TMarker(x, y, 20);
            mXY->SetMarkerColor(color);
            mXY->SetMarkerSize(0.8);
            mXY->Draw();
            
            canvas->cd(2);
            TMarker* mXZ = new TMarker(x, z, 20);
            mXZ->SetMarkerColor(color);
            mXZ->SetMarkerSize(0.8);
            mXZ->Draw();
            
            canvas->cd(3);
            TMarker* mYZ = new TMarker(z, y, 20);
            mYZ->SetMarkerColor(color);
            mYZ->SetMarkerSize(0.8);
            mYZ->Draw();
            
            canvas->cd(4);
            TPolyMarker3D* pm3d = new TPolyMarker3D(1);
            pm3d->SetPoint(0, x, z, y);
            pm3d->SetMarkerColor(color);
            pm3d->SetMarkerStyle(20);
            pm3d->SetMarkerSize(0.8);
            pm3d->Draw();
            
            pointsDrawn++;
        }
    }
    
    // Show steps in the current time window
    for (int i = 0; i < event->steps.nSteps; i++) {
        float time = event->steps.time[i];
        
        // Show only steps in the current time window
        if (time >= currentTime - timeWindow && time <= currentTime + timeWindow) {
            float x = event->steps.posX[i];
            float y = event->steps.posY[i];
            float z = event->steps.posZ[i];
            float energy = event->steps.energy[i];
            
            Int_t color = GetColorByEnergy(energy);
            
            canvas->cd(1);
            TMarker* mXY = new TMarker(x, y, 20);
            mXY->SetMarkerColor(color);
            mXY->SetMarkerSize(0.8);
            mXY->Draw();
            
            canvas->cd(2);
            TMarker* mXZ = new TMarker(x, z, 20);
            mXZ->SetMarkerColor(color);
            mXZ->SetMarkerSize(0.8);
            mXZ->Draw();
            
            canvas->cd(3);
            TMarker* mYZ = new TMarker(z, y, 20);
            mYZ->SetMarkerColor(color);
            mYZ->SetMarkerSize(0.8);
            mYZ->Draw();
            
            canvas->cd(4);
            TPolyMarker3D* pm3d = new TPolyMarker3D(1);
            pm3d->SetPoint(0, x, z, y);
            pm3d->SetMarkerColor(color);
            pm3d->SetMarkerStyle(20);
            pm3d->SetMarkerSize(0.8);
            pm3d->Draw();
            
            pointsDrawn++;
        }
    }
    
    canvas->cd(1);
    TPaveText* info = new TPaveText(0.15, 0.85, 0.55, 0.90, "NDC");
    info->SetFillColor(0);
    info->SetTextAlign(12);
    info->SetTextSize(0.035);
    info->AddText(Form("Time: %.2f ps (Window: %.2f ps)", currentTime, timeWindow));
    info->Draw();
    
    TPaveText* stats = new TPaveText(0.15, 0.80, 0.55, 0.85, "NDC");
    stats->SetFillColor(0);
    stats->SetTextAlign(12);
    stats->SetTextSize(0.035);
    stats->AddText(Form("Points drawn: %d", pointsDrawn));
    stats->Draw();
    
    return canvas;
} 

TCanvas* MCPVisualizer::DrawMCP2DWithPoresAndSteps() {
    // 1. 캔버스/프레임 생성
    TCanvas* c = new TCanvas("c_mcp2d", "MCP1/2 Pores and Steps", 900, 600);
    TH2F* frame = new TH2F("frame", "MCP1/2 Pores and Steps;X [#mum];Y [#mum]", 100, 0, 3500, 100, -100, 300);
    frame->SetStats(0);
    frame->Draw();

    // 2. config 파라미터 읽기
    const mcp::ConfigParameters* config = analyzer_->GetConfig();
    if (!config) return c;
    float x0 = config->x0, x1 = config->x1, x2 = config->x2, x3 = config->x3, x4 = config->x4;
    float alpha1 = config->alpha1, alpha2 = config->alpha2;
    float dia = config->dia;
    float pas = config->pas;
    float r = dia / 2.0;
    float pitch = dia + pas;

    // 3. MCP1, MCP2 포어 경계 그리기
    float pore_center_y0 = 6.0; // 중심 포어 y=6.0에 맞춤
    auto clip_line = [](float x0, float y0, float x1, float y1, float y_min, float y_max, float& cx0, float& cy0, float& cx1, float& cy1) {
        cx0 = x0; cy0 = y0;
        cx1 = x1; cy1 = y1;
        // y0가 아래로 벗어나면
        if (cy0 < y_min) {
            cx0 = x0 + (x1 - x0) * (y_min - y0) / (y1 - y0);
            cy0 = y_min;
        }
        // y0가 위로 벗어나면
        if (cy0 > y_max) {
            cx0 = x0 + (x1 - x0) * (y_max - y0) / (y1 - y0);
            cy0 = y_max;
        }
        // y1가 아래로 벗어나면
        if (cy1 < y_min) {
            cx1 = x0 + (x1 - x0) * (y_min - y0) / (y1 - y0);
            cy1 = y_min;
        }
        // y1가 위로 벗어나면
        if (cy1 > y_max) {
            cx1 = x0 + (x1 - x0) * (y_max - y0) / (y1 - y0);
            cy1 = y_max;
        }
    };
    auto draw_pore_array = [&](float x_start, float x_end, float alpha, int color) {
        float y_min = -100.0, y_max = 300.0;
        int n_min = static_cast<int>(std::floor((y_min - pore_center_y0) / pitch)) - 2;
        int n_max = static_cast<int>(std::ceil((y_max - pore_center_y0) / pitch)) + 2;
        for (int n = n_min; n <= n_max; ++n) {
            float y0 = pore_center_y0 + n * pitch;
            float y0_top = y0 + r, y0_bot = y0 - r;
            float y1_top = y0 + r + (x_end - x_start) * tan(alpha);
            float y1_bot = y0 - r + (x_end - x_start) * tan(alpha);
            // top wall
            float cx0, cy0, cx1, cy1;
            clip_line(x_start, y0_top, x_end, y1_top, y_min, y_max, cx0, cy0, cx1, cy1);
            if ((cy0 >= y_min && cy0 <= y_max) || (cy1 >= y_min && cy1 <= y_max)) {
                TLine* l1 = new TLine(cx0, cy0, cx1, cy1);
                l1->SetLineColor(color);
                l1->SetLineWidth(1);
                l1->Draw();
            }
            // bottom wall
            clip_line(x_start, y0_bot, x_end, y1_bot, y_min, y_max, cx0, cy0, cx1, cy1);
            if ((cy0 >= y_min && cy0 <= y_max) || (cy1 >= y_min && cy1 <= y_max)) {
                TLine* l2 = new TLine(cx0, cy0, cx1, cy1);
                l2->SetLineColor(color);
                l2->SetLineWidth(1);
                l2->Draw();
            }
        }
    };
    draw_pore_array(x0, x1, alpha1, kGray+1);
    draw_pore_array(x2, x3, alpha2, kGray+1);

    // 4. 모든 step 누적 포인트 찍기
    const mcp::Event* event = analyzer_->GetEvent();
    if (event) {
        for (int i = 0; i < event->steps.nSteps; ++i) {
            float x = event->steps.posX[i];
            float y = event->steps.posY[i];
            TMarker* m = new TMarker(x, y, 7);
            m->SetMarkerColor(kBlack);
            m->SetMarkerStyle(20);
            m->SetMarkerSize(0.3);
            m->Draw();
        }
    }
    return c;
} 

// ----------------------------------------------------------------------
// Draw pores and steps for a subset of tracks (specified by trackID list)
// ----------------------------------------------------------------------
TCanvas* MCPVisualizer::DrawMCP2DForTracks(const std::vector<int>& trackIDs){
    // Build a lookup set for fast membership test
    std::unordered_set<int> sel(trackIDs.begin(), trackIDs.end());
 
    // 1. canvas/frame (reuse same ranges)
    TCanvas* c = new TCanvas("c_mcp2d_sel", "MCP Pores and Steps (selected)", 900,600);
    TH2F* frame = new TH2F("f_sel","MCP Pores and Steps;X [#mum];Y [#mum]",100,0,3500,100,-100,300);
    frame->SetStats(0);
    frame->Draw();
 
    // 2. draw pore boundaries (reuse helper lambda)
    const mcp::ConfigParameters* config = analyzer_->GetConfig();
    if(!config) return c;
    float x0=config->x0,x1=config->x1,x2=config->x2,x3=config->x3;
    float alpha1=config->alpha1,alpha2=config->alpha2;
    float dia=config->dia,pas=config->pas,r=dia/2.0,pitch=dia+pas;
    float pore_center_y0=6.0;
    auto clip_line=[&](float x0,float y0,float x1,float y1,float y_min,float y_max,
                        float& cx0,float& cy0,float& cx1,float& cy1){
         cx0=x0;cy0=y0;cx1=x1;cy1=y1;
         if(cy0<y_min){cx0=x0+(x1-x0)*(y_min-y0)/(y1-y0);cy0=y_min;}
         if(cy0>y_max){cx0=x0+(x1-x0)*(y_max-y0)/(y1-y0);cy0=y_max;}
         if(cy1<y_min){cx1=x0+(x1-x0)*(y_min-y0)/(y1-y0);cy1=y_min;}
         if(cy1>y_max){cx1=x0+(x1-x0)*(y_max-y0)/(y1-y0);cy1=y_max;}
     };
     auto draw_pore=[&](float xs,float xe,float alpha){
         float y_min=-100,y_max=300;
         int n_min=int(std::floor((y_min-pore_center_y0)/pitch))-2;
         int n_max=int(std::ceil((y_max-pore_center_y0)/pitch))+2;
         for(int n=n_min;n<=n_max;++n){
             float y0=pore_center_y0+n*pitch;
             float y0t=y0+r, y0b=y0-r;
             float y1t=y0+r+(xe-xs)*std::tan(alpha);
             float y1b=y0-r+(xe-xs)*std::tan(alpha);
             float cx0,cy0,cx1,cy1;
             clip_line(xs,y0t,xe,y1t,y_min,y_max,cx0,cy0,cx1,cy1);
             if((cy0>=y_min&&cy0<=y_max)||(cy1>=y_min&&cy1<=y_max)){
                 TLine* l=new TLine(cx0,cy0,cx1,cy1); l->SetLineColor(kGray+1); l->Draw();}
             clip_line(xs,y0b,xe,y1b,y_min,y_max,cx0,cy0,cx1,cy1);
             if((cy0>=y_min&&cy0<=y_max)||(cy1>=y_min&&cy1<=y_max)){
                 TLine* l=new TLine(cx0,cy0,cx1,cy1); l->SetLineColor(kGray+1); l->Draw();}
         }
     };
     draw_pore(x0,x1,alpha1);
     draw_pore(x2,x3,alpha2);
 
     // 3. plot only steps whose trackID in sel
     const mcp::Event* evt = analyzer_->GetEvent();
     if(evt){
         for(int i=0;i<evt->steps.nSteps;++i){
             if(sel.find(evt->steps.trackID[i])==sel.end()) continue;
             float x=evt->steps.posX[i];
             float y=evt->steps.posY[i];
             TMarker* m=new TMarker(x,y,7);
             m->SetMarkerColor(kBlack);
             m->SetMarkerStyle(20);
             m->SetMarkerSize(0.3);
             m->Draw();
         }
     }
     return c;
 } 

// ----------------------------------------------------------------------
// Overlay cascades: draw pores once, then each track-set with different color
// ----------------------------------------------------------------------
TCanvas* MCPVisualizer::DrawMCP2DOverlay(const std::vector<std::vector<int>>& trackSets){
    // base frame
    TCanvas* c = new TCanvas("c_mcp2d_overlay","Cascade Overlay",900,600);
    TH2F* frame = new TH2F("f_overlay","MCP Cascades;X [#mum];Y [#mum]",100,0,3500,100,-100,300);
    frame->SetStats(0);
    frame->Draw();

    // pore boundaries (reuse from earlier lambda)
    const mcp::ConfigParameters* cfg = analyzer_->GetConfig();
    if(!cfg) return c;
    float x0=cfg->x0,x1=cfg->x1,x2=cfg->x2,x3=cfg->x3;
    float alpha1=cfg->alpha1,alpha2=cfg->alpha2;
    float dia=cfg->dia,pas=cfg->pas,r=dia/2.0,pitch=dia+pas;
    float y0c=6.0;
    auto clip=[&](float x0,float y0,float x1,float y1,float ymin,float ymax,float& cx0,float& cy0,float& cx1,float& cy1){
        cx0=x0;cy0=y0;cx1=x1;cy1=y1;
        if(cy0<ymin){cx0=x0+(x1-x0)*(ymin-y0)/(y1-y0);cy0=ymin;}
        if(cy0>ymax){cx0=x0+(x1-x0)*(ymax-y0)/(y1-y0);cy0=ymax;}
        if(cy1<ymin){cx1=x0+(x1-x0)*(ymin-y0)/(y1-y0);cy1=ymin;}
        if(cy1>ymax){cx1=x0+(x1-x0)*(ymax-y0)/(y1-y0);cy1=ymax;}
    };
    auto drawP=[&](float xs,float xe,float alpha){
        float ymin=-100, ymax=300;
        int nmin=int(std::floor((ymin-y0c)/pitch))-2;
        int nmax=int(std::ceil((ymax-y0c)/pitch))+2;
        for(int n=nmin;n<=nmax;++n){
            float y00=y0c+n*pitch;
            float yt0=y00+r, yb0=y00-r;
            float yt1=y00+r+(xe-xs)*std::tan(alpha);
            float yb1=y00-r+(xe-xs)*std::tan(alpha);
            float cx0,cy0,cx1,cy1;
            clip(xs,yt0,xe,yt1,ymin,ymax,cx0,cy0,cx1,cy1);
            if((cy0>=ymin&&cy0<=ymax)||(cy1>=ymin&&cy1<=ymax)){
                TLine* l = new TLine(cx0,cy0,cx1,cy1);
                l->Draw();
            }
            clip(xs,yb0,xe,yb1,ymin,ymax,cx0,cy0,cx1,cy1);
            if((cy0>=ymin&&cy0<=ymax)||(cy1>=ymin&&cy1<=ymax)){
                TLine* l2 = new TLine(cx0,cy0,cx1,cy1);
                l2->Draw();
            }
        }
    };
    drawP(x0,x1,alpha1);
    drawP(x2,x3,alpha2);

    // color list
    const int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+1, kViolet};
    int nColors = sizeof(colors)/sizeof(int);

    const mcp::Event* evt = analyzer_->GetEvent();
    if(!evt) return c;

    for(size_t s=0;s<trackSets.size();++s){
        int col = colors[s % nColors];
        std::unordered_set<int> sel(trackSets[s].begin(), trackSets[s].end());
        for(int i=0;i<evt->steps.nSteps;++i){
            if(sel.find(evt->steps.trackID[i])==sel.end()) continue;
            TMarker* m = new TMarker(evt->steps.posX[i], evt->steps.posY[i], 20);
            m->SetMarkerColor(col);
            m->SetMarkerSize(0.4);
            m->Draw();
        }
    }
    return c;
} 

// ----------------------------------------------------------------------
// 3D overlay view of cascades
// ----------------------------------------------------------------------
TCanvas* MCPVisualizer::DrawMCP3DOverlay(const std::vector<std::vector<int>>& trackSets){
    // Create canvas
    TCanvas* c = new TCanvas("c_mcp3d_overlay","MCP Cascades (3D)", 900, 700);
    c->cd();

    // Axis ranges – keep X same as 2D frame, Y and Z identical ranges for square aspect
    double xMin = 0.0, xMax = 3500.0;
    double yMin = -600.0, yMax = 600.0;
    double zMin = -600.0, zMax = 600.0;

    // Dummy histogram just to draw the 3-D box & axes
    TH3F* hFrame = new TH3F("h3d_frame","MCP Cascades (3D);X [#mum];Z [#mum];Y [#mum]",
                           10,xMin,xMax,
                           10,zMin,zMax,
                           10,yMin,yMax);
    hFrame->SetStats(0);
    hFrame->Draw("BOX");

    // Set some view angles for initial inspection (user can rotate interactively)
    gPad->SetTheta(25);
    gPad->SetPhi(35);

    // Colour palette for different cascades (same as 2D)
    const int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+1, kViolet};
    int nColors = sizeof(colors)/sizeof(int);

    // Access event data
    const mcp::Event* evt = analyzer_->GetEvent();
    if(!evt) return c;

    // Build look-up sets for quick membership testing per cascade
    std::vector<std::unordered_set<int>> cascadeSel;
    cascadeSel.reserve(trackSets.size());
    for(const auto& v : trackSets){
        cascadeSel.emplace_back(v.begin(), v.end());
    }

    // Iterate over all steps once, draw marker according to membership
    for(int i=0;i<evt->steps.nSteps;++i){
        int trackId = evt->steps.trackID[i];
        // Determine which cascade this step belongs to
        int cascadeIdx = -1;
        for(size_t cIdx=0;cIdx<cascadeSel.size();++cIdx){
            if(cascadeSel[cIdx].find(trackId)!=cascadeSel[cIdx].end()){
                cascadeIdx = static_cast<int>(cIdx);
                break;
            }
        }
        if(cascadeIdx<0) continue; // step not in any selected cascade

        // Grab coordinates
        double x = evt->steps.posX[i];
        double y = evt->steps.posY[i];
        double z = evt->steps.posZ[i];

        // Apply basic zoom cuts (optional): keep within axis limits
        if(x < xMin || x > xMax) continue;
        if(y < yMin || y > yMax) continue;
        if(z < zMin || z > zMax) continue;

        // Draw marker – re-use simple 1-point poly-marker to avoid memory of large arrays
        TPolyMarker3D* pm = new TPolyMarker3D(1);
        pm->SetPoint(0,x,z,y);
        pm->SetMarkerColor(colors[cascadeIdx % nColors]);
        pm->SetMarkerStyle(20);
        pm->SetMarkerSize(0.6);
        pm->Draw();
    }

    // Optionally draw pore walls as grey hatched regions along full Z span for context
    const mcp::ConfigParameters* cfg = analyzer_->GetConfig();
    if(cfg){
        double x0=cfg->x0, x1=cfg->x1, x2=cfg->x2, x3=cfg->x3;
        double alpha1=cfg->alpha1, alpha2=cfg->alpha2;
        double dia = cfg->dia;
        double r = dia/2.0;
        double pitch = dia + cfg->pas;
        double yMid=6.0;

        // Helper lambda: draw filled planar wall (vertical rectangle extruded along Z)
        auto drawWall = [&](double xs, double ys, double xe, double ye){
            // Build rectangle (xs,ys) -> (xe,ye) and extrude between zMin,zMax
            const int nPts = 5; // close the polygon
            TPolyLine3D* surf = new TPolyLine3D(nPts);
            surf->SetPoint(0, xs, zMin, ys);
            surf->SetPoint(1, xs, zMax, ys);
            surf->SetPoint(2, xe, zMax, ye);
            surf->SetPoint(3, xe, zMin, ye);
            surf->SetPoint(4, xs, zMin, ys);

            // Plain grey outline (no fill) to avoid any transparency path
            surf->SetLineColor(kGray+1);
            surf->SetLineWidth(1);
            surf->Draw();
        };

        // Draw walls for several pore rows in range
        int nMin = static_cast<int>(std::floor((yMin - yMid)/pitch)) - 2;
        int nMax = static_cast<int>(std::ceil((yMax - yMid)/pitch)) + 2;
        for(int n=nMin; n<=nMax; ++n){
            double y0 = yMid + n*pitch;
            // MCP1
            double ysTop = y0 + r;
            double ysBot = y0 - r;
            double yeTop = y0 + r + (x1 - x0)*std::tan(alpha1);
            double yeBot = y0 - r + (x1 - x0)*std::tan(alpha1);
            drawWall(x0, ysTop, x1, yeTop);
            drawWall(x0, ysBot, x1, yeBot);
            // MCP2
            ysTop = y0 + r;
            ysBot = y0 - r;
            yeTop = y0 + r + (x3 - x2)*std::tan(alpha2);
            yeBot = y0 - r + (x3 - x2)*std::tan(alpha2);
            drawWall(x2, ysTop, x3, yeTop);
            drawWall(x2, ysBot, x3, yeBot);
        }
    }

    return c;
} 

// ----------------------------------------------------------------------
// Zoomed-in 3D view focusing on a handful of pores around centre
// ----------------------------------------------------------------------
TCanvas* MCPVisualizer::DrawMCP3DZoom(const std::vector<std::vector<int>>& trackSets,
                                      int nRange, int nzRange){
    // Configuration parameters
    const mcp::ConfigParameters* cfg = analyzer_->GetConfig();
    if(!cfg) return nullptr;
    const double x0 = cfg->x0, x1 = cfg->x1, x2 = cfg->x2, x3 = cfg->x3;
    const double alpha1 = cfg->alpha1, alpha2 = cfg->alpha2;
    const double dia = cfg->dia;
    const double pitch = dia + cfg->pas;
    const double R = dia/2.0;
    const double y0c = 6.0;            // central pore row reference (n=0)

    // Axis limits (tight)
    double xMin = x0 - 20.0;           // entrance of MCP1 a bit before
    double xMax = x3 + 20.0;           // exit of MCP2 a bit after
    double zMin = y0c + (-nRange-1)*pitch - R*1.5;
    double zMax = y0c + ( nRange+1)*pitch + R*1.5;
    double yMin = (-nzRange-1)*pitch - R*1.5;
    double yMax = ( nzRange+1)*pitch + R*1.5;

    // Canvas (no bounding TH3F – geometry itself provides context)
    TCanvas* c = new TCanvas("c_mcp3d_zoom","MCP Zoom", 900, 700);
    // We avoid relying on camera/view rotation (problematic on some ROOT builds)
    // Instead, we will rotate the entire geometry by +90° around X so that
    // the default view (phi=30°, theta=30°, roll=0°) already shows the desired
    // chevron orientation.  If you need a different orientation, simply adjust
    // the rotation matrix below (RotateX / Y / Z).

    // Build colour mapping
    const int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+1, kViolet};
    const int nColors = sizeof(colors)/sizeof(int);

    // Prepare cascade membership sets
    std::vector<std::unordered_set<int>> selSets;
    selSets.reserve(trackSets.size());
    for(const auto& vec : trackSets){ selSets.emplace_back(vec.begin(), vec.end()); }

    // ------------------------------------------------------------------
    //  Build geometry with a single reusable TGeoManager (avoid multiple
    //  managers which often leads to dangling pointers inside GL viewer)
    // ------------------------------------------------------------------

    // Ensure a single global GeoManager instance is used ------------------------------------------------
    TGeoManager* geo = nullptr;
    if (gGeoManager) {
        // Re-use existing manager after cleaning previous geometry
        gGeoManager->GetListOfVolumes()->Delete();      // delete old volumes
        gGeoManager->GetListOfShapes()->Delete();       // delete old shapes
        gGeoManager->GetListOfMatrices()->Delete();     // delete old matrices
        gGeoManager->SetTopVolume(nullptr);             // detach previous top
        geo = gGeoManager;
    } else {
        // First time: create the global manager
        geo = new TGeoManager("mcp_geo","MCP geometry");
    }

    // simple vacuum material/medium
    TGeoMaterial* matVac = new TGeoMaterial("vacuum", 0,0,0);
    TGeoMedium*  medVac = new TGeoMedium("vac", 1, matVac);

    // world box (half-lengths a bit larger than view frustum)
    double wdx = (xMax - xMin)/2 + 50;
    double wdy = (yMax - yMin)/2 + 50;
    double wdz = (zMax - zMin)/2 + 50;
    // Build visible world box --------------------------------------------------------------------------
    TGeoVolume* world = geo->MakeBox("world", medVac, wdx, wdy, wdz);
    world->SetLineColor(kGray+1); // visible wireframe

    // ------------------------------------------------------------------
    // Wrap the world volume in a rotated assembly so that the default
    // camera angles already show the desired chevron orientation.
    // ------------------------------------------------------------------

    // Rotation: +90 deg around global X (swap Y/Z)
    auto *rotScene = new TGeoRotation();
    rotScene->RotateX(90);       // adjust as needed

    auto *trScene  = new TGeoCombiTrans(0, 0, 0, rotScene);

    // Assembly volume acts as new top container
    TGeoVolume *top = geo->MakeVolumeAssembly("top");
    top->AddNode(world, 1, trScene);
    geo->SetTopVolume(top);

    // Helper lambda to add one MCP section (tube array)
    auto addMcpTubes = [&](double xs, double xe, double alpha){
        double len = xe - xs;                  // physical length along X
        double halfLen = len / 2.0;

        // base tube volume (axis along global X after rotation)
        TString tubeName = Form("tube_%.0f", xs);
        TGeoVolume* tube = geo->MakeTube(tubeName, medVac, 0, R, halfLen);
        tube->SetLineColor(kGray+2);
        tube->SetFillColor(kGray+2);          // semi-transparent faces
        // Removed transparency to avoid libAfterImage segfaults
        // tube->SetTransparency(90);            // 0=opaque,100=invisible (30% visible)

        // build rotation: local Z (tube axis) -> global X; then tilt by alpha
        auto* rot = new TGeoRotation();
        rot->RotateY(90);                     // align local Z to global X
        rot->RotateZ(alpha*180.0/TMath::Pi());       // apply pore tilt

        // place tubes in requested Y/Z range
        for(int n=-nRange; n<=nRange; ++n){
            for(int nz=-nzRange; nz<=nzRange; ++nz){
                double xMid = xs + halfLen;
                double zMid = y0c + n*pitch + std::tan(alpha)*(xMid - xs);
                double yMid = nz * pitch;

                auto* comb = new TGeoCombiTrans(xMid, yMid, zMid, rot);
                comb->RegisterYourself();           // needed if reused
                world->AddNode(tube, (n+nRange)*(2*nzRange+1)+ (nz+nzRange), comb);
            }
        }
    };

    addMcpTubes(x0, x1, alpha1);
    addMcpTubes(x2, x3, alpha2);

    // Keep world box but draw only its wireframe
    // world->SetLineColor(kGray+1);
    // world->SetFillStyle(0);

    // Finalise and draw geometry -----------------------------------------------------------------------
    geo->CloseGeometry();

    // Draw geometry first to initialise GL viewer (use rotated top volume)
    top->Draw("gl");
    gPad->Update();   // ensure TGLViewer is created

    // Pad will be updated later once geometry is drawn (see below).

    // ------------------------------------------------------------------
    // Draw helper XYZ axes (X-Z-Y order)
    // ------------------------------------------------------------------
    auto drawAxis = [&](double x0,double y0,double z0,
                         double x1,double y1,double z1,
                         Color_t col){
        TPolyLine3D* ax = new TPolyLine3D(2);
        ax->SetPoint(0,x0,y0,z0);
        ax->SetPoint(1,x1,y1,z1);
        ax->SetLineColor(col);
        ax->SetLineWidth(2);
        ax->Draw();
    };

    // Choose an axis origin near the world corner
    double axX0 = xMin;
    double axY0 = yMin;
    double axZ0 = zMin;
    double axLen = 0.30*(xMax - xMin);

    // X-axis (red) - horizontal right
    drawAxis(axX0, axY0, axZ0,  axX0+axLen, axY0, axZ0,  kRed);
    // Y-axis (green) - vertical down (reversed)
    drawAxis(axX0, axY0, axZ0,  axX0, axY0-axLen, axZ0, kGreen+2);
    // Z-axis (blue) - vertical up
    drawAxis(axX0, axY0, axZ0,  axX0, axY0, axZ0+axLen, kBlue);

    // After geometry is visible, overlay electron steps inside zoom box
    const mcp::Event* evt = analyzer_->GetEvent();
    if(evt){
        for(int i=0;i<evt->steps.nSteps;++i){
            double x = evt->steps.posX[i];
            double z = evt->steps.posY[i];
            double y = evt->steps.posZ[i];
            if(x<xMin||x>xMax||y<yMin||y>yMax||z<zMin||z>zMax) continue;

            int trackId = evt->steps.trackID[i];
            int cIdx=-1;
            for(size_t k=0;k<selSets.size();++k)
                if(selSets[k].count(trackId)){ cIdx=k; break; }
            if(cIdx<0) continue;
            TPolyMarker3D* pm = new TPolyMarker3D(1);
            pm->SetPoint(0,x,y,z);
            pm->SetMarkerStyle(20);
            pm->SetMarkerSize(0.6);
            pm->SetMarkerColor(colors[cIdx % nColors]);
            pm->Draw();
        }
    }

    std::cout << "MCP1 중심 Y(after rot) = "
          << - (y0c +  std::tan(alpha1)*(x1-x0)/2)
          << "\nMCP2 중심 Y(after rot) = "
          << - (y0c +  std::tan(alpha2)*(x3-x2)/2) << std::endl;

    // Draw thin reference box atop everything
    gPad->Modified(); gPad->Update();

    // Draw thin wireframe box for orientation (12 edges)
    auto drawEdge = [&](double ax,double ay,double az,double bx,double by,double bz){
        TPolyLine3D* line=new TPolyLine3D(2);
        line->SetPoint(0,ax,ay,az);
        line->SetPoint(1,bx,by,bz);
        line->SetLineColor(kGray+1);
        line->SetLineWidth(1);
        line->Draw();
    };
    // bottom rectangle (z=zMin)
    drawEdge(xMin,yMin,zMin, xMax,yMin,zMin);
    drawEdge(xMax,yMin,zMin, xMax,yMax,zMin);
    drawEdge(xMax,yMax,zMin, xMin,yMax,zMin);
    drawEdge(xMin,yMax,zMin, xMin,yMin,zMin);
    // top rectangle (z=zMax)
    drawEdge(xMin,yMin,zMax, xMax,yMin,zMax);
    drawEdge(xMax,yMin,zMax, xMax,yMax,zMax);
    drawEdge(xMax,yMax,zMax, xMin,yMax,zMax);
    drawEdge(xMin,yMax,zMax, xMin,yMin,zMax);
    // vertical edges
    drawEdge(xMin,yMin,zMin, xMin,yMin,zMax);
    drawEdge(xMax,yMin,zMin, xMax,yMin,zMax);
    drawEdge(xMax,yMax,zMin, xMax,yMax,zMax);
    drawEdge(xMin,yMax,zMin, xMin,yMax,zMax);

    return c;
} 
