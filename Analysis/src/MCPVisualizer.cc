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
#include "TView3D.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TBox.h"

MCPVisualizer::MCPVisualizer(const MCPAnalyzer* analyzer) : analyzer_(analyzer) {
    // Set ROOT style
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
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
                        100, minX, maxX, 100, -100, 300);
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
                        100, minX, maxX, 100, minZ, maxZ);
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
                        100, minZ, maxZ, 100, -100, 300);  
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
                        10, minX, maxX, 10, -100, 300, 10, minZ, maxZ);
    h3D->SetStats(0);  
    h3D->GetXaxis()->SetTitle("X (#mum)");
    h3D->GetYaxis()->SetTitle("Y (#mum)");
    h3D->GetZaxis()->SetTitle("Z (#mum)");
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
            pm3d->SetPoint(0, x, y, z);
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
            pm3d->SetPoint(0, x, y, z);
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

    std::cout << "alpha1: " << alpha1 << ", alpha2: " << alpha2 << std::endl;

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