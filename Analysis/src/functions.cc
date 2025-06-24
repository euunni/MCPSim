#include "../include/functions.h"
#include "../../MCPSim/include/MCPEventData.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>

void loadRootFile(const std::string& filename, std::vector<mcp::Event*>& eventPtrs) {
    TFile* file = TFile::Open(("../output/" + filename + ".root").c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }
    
    TTree* tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "Cannot find tree 'Events' in file" << std::endl;
        file->Close();
        return;
    }
    
    // 브랜치 변수 설정 (포인터로 설정)
    mcp::EventInfo* eventInfo = nullptr;
    mcp::ConfigParameters* config = nullptr;
    mcp::Track* tracks = nullptr;
    mcp::Step* steps = nullptr;
    
    // 브랜치 주소 설정
    tree->SetBranchAddress("EventInfo", &eventInfo);
    tree->SetBranchAddress("Config", &config);
    tree->SetBranchAddress("Tracks", &tracks);
    tree->SetBranchAddress("Steps", &steps);

    int nEntries = tree->GetEntries();
    eventPtrs.reserve(nEntries);

    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // 새 이벤트 객체 생성
        mcp::Event* newEvent = new mcp::Event();
        
        // 포인터 할당 (복사 없음)
        if (eventInfo) newEvent->eventInfo = *eventInfo;
        if (config) newEvent->config = *config;
        if (tracks) newEvent->tracks = *tracks;
        if (steps) newEvent->steps = *steps;
        
        eventPtrs.push_back(newEvent);
    }

    file->Close();
    
    // Check data range
    if (!eventPtrs.empty()) {
        const mcp::Event* firstEvent = eventPtrs[0];
        
        float minX = std::numeric_limits<float>::max();
        float maxX = -std::numeric_limits<float>::max();
        float minY = std::numeric_limits<float>::max();
        float maxY = -std::numeric_limits<float>::max();
        float minZ = std::numeric_limits<float>::max();
        float maxZ = -std::numeric_limits<float>::max();
        
        // 트랙 정보 사용
        for (int i = 0; i < firstEvent->tracks.nTracks; i++) {
            // 시작 위치
            minX = std::min(minX, firstEvent->tracks.birthPosX[i]);
            maxX = std::max(maxX, firstEvent->tracks.birthPosX[i]);
            minY = std::min(minY, firstEvent->tracks.birthPosY[i]);
            maxY = std::max(maxY, firstEvent->tracks.birthPosY[i]);
            minZ = std::min(minZ, firstEvent->tracks.birthPosZ[i]);
            maxZ = std::max(maxZ, firstEvent->tracks.birthPosZ[i]);
            
            // 최종 위치
            minX = std::min(minX, firstEvent->tracks.finalPosX[i]);
            maxX = std::max(maxX, firstEvent->tracks.finalPosX[i]);
            minY = std::min(minY, firstEvent->tracks.finalPosY[i]);
            maxY = std::max(maxY, firstEvent->tracks.finalPosY[i]);
            minZ = std::min(minZ, firstEvent->tracks.finalPosZ[i]);
            maxZ = std::max(maxZ, firstEvent->tracks.finalPosZ[i]);
        }
        
        // 스텝 정보 사용
        for (int i = 0; i < firstEvent->steps.nSteps; i++) {
            minX = std::min(minX, firstEvent->steps.posX[i]);
            maxX = std::max(maxX, firstEvent->steps.posX[i]);
            minY = std::min(minY, firstEvent->steps.posY[i]);
            maxY = std::max(maxY, firstEvent->steps.posY[i]);
            minZ = std::min(minZ, firstEvent->steps.posZ[i]);
            maxZ = std::max(maxZ, firstEvent->steps.posZ[i]);
        }
        
        std::cout << "Data range: X[" << minX << ", " << maxX << "], "
                  << "Y[" << minY << ", " << maxY << "], "
                  << "Z[" << minZ << ", " << maxZ << "]" << std::endl;
    }
}

// Set color based on energy
Int_t getEnergyColor(float energy, float minEnergy, float maxEnergy) {
    // Normalize energy to 0~1
    float normalizedEnergy = (energy - minEnergy) / (maxEnergy - minEnergy);
    if (normalizedEnergy < 0) normalizedEnergy = 0;
    if (normalizedEnergy > 1) normalizedEnergy = 1;
    
    gStyle->SetPalette(kViridis);
     
    Int_t colorIndex = static_cast<Int_t>(normalizedEnergy * 255);
    Int_t color = TColor::GetColorPalette(colorIndex);
    
    return color;
}

void drawPoreBoundary(TVirtualPad* pad, float centerX, float centerY, float minX, float maxX, float minZ, float maxZ) {
    pad->cd();
    
    const float poreStartY = 6.0;  // Y center coordinate at the start point
    const float poreDiameter = 10.0;  // Pore diameter
    const float poreAngle = 0.13;  // Pore angle (radian)
    const float poreRadius = poreDiameter / 2.0;  // Pore radius
    const float poreLength = 400.0;  // Pore x-component length (μm)
    
    // Pore start and end point coordinates
    const float poreStartX = minX + (maxX - minX) * 0.1;  // Pore start X coordinate (10% from the left)
    const float poreEndX = poreStartX + poreLength;       // Pore end X coordinate (start point + length)
    
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

// Create electron cascade tree
ElectronCascadeTree buildElectronCascadeTree(const mcp::Event& event) {
    ElectronCascadeTree tree;
    
    // 트랙 정보 사용
    for (int i = 0; i < event.tracks.nTracks; i++) {
        int trackID = event.tracks.trackID[i];
        int parentID = event.tracks.parentID[i];
        
        // Check if it already exists
        if (tree.find(trackID) == tree.end()) {
            ElectronNode node;
            node.trackID = trackID;
            node.parentID = parentID;
            node.posX = event.tracks.birthPosX[i];
            node.posY = event.tracks.birthPosY[i];
            node.posZ = event.tracks.birthPosZ[i];
            node.time = event.tracks.birthTime[i];
            node.energy = event.tracks.birthEnergy[i];
            tree[trackID] = node;
        }
        
        // Add as a child to the parent node
        if (parentID > 0 && tree.find(parentID) != tree.end()) {
            if (std::find(tree[parentID].childIDs.begin(), 
                         tree[parentID].childIDs.end(), 
                         trackID) == tree[parentID].childIDs.end()) {
                tree[parentID].childIDs.push_back(trackID);
            }
        }
    }
    
    std::cout << "전자 증폭 트리 생성 완료. 총 전자 수: " << event.tracks.nTracks << std::endl;
    
    return tree;
}

// Calculate data range
void calculateDataRange(const mcp::Event& event, 
                       float& minX, float& maxX, 
                       float& minY, float& maxY, 
                       float& minZ, float& maxZ) {
    minX = minY = minZ = std::numeric_limits<float>::max();
    maxX = maxY = maxZ = -std::numeric_limits<float>::max();
    
    // 트랙 정보 사용
    for (int i = 0; i < event.tracks.nTracks; i++) {
        // 시작 위치
        minX = std::min(minX, event.tracks.birthPosX[i]);
        maxX = std::max(maxX, event.tracks.birthPosX[i]);
        minY = std::min(minY, event.tracks.birthPosY[i]);
        maxY = std::max(maxY, event.tracks.birthPosY[i]);
        minZ = std::min(minZ, event.tracks.birthPosZ[i]);
        maxZ = std::max(maxZ, event.tracks.birthPosZ[i]);
        
        // 최종 위치
        minX = std::min(minX, event.tracks.finalPosX[i]);
        maxX = std::max(maxX, event.tracks.finalPosX[i]);
        minY = std::min(minY, event.tracks.finalPosY[i]);
        maxY = std::max(maxY, event.tracks.finalPosY[i]);
        minZ = std::min(minZ, event.tracks.finalPosZ[i]);
        maxZ = std::max(maxZ, event.tracks.finalPosZ[i]);
    }
    
    // 스텝 정보 사용
    for (int i = 0; i < event.steps.nSteps; i++) {
        minX = std::min(minX, event.steps.posX[i]);
        maxX = std::max(maxX, event.steps.posX[i]);
        minY = std::min(minY, event.steps.posY[i]);
        maxY = std::max(maxY, event.steps.posY[i]);
        minZ = std::min(minZ, event.steps.posZ[i]);
        maxZ = std::max(maxZ, event.steps.posZ[i]);
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

// 2D histogram: electron density per location in the pore
void createDensityHistogram(const mcp::Event& event, const std::string& outputFileName) {
    // Calculate data range
    float minX, maxX, minY, maxY, minZ, maxZ;
    calculateDataRange(event, minX, maxX, minY, maxY, minZ, maxZ);
    
    // Calculate pore center
    float poreCenterX = (minX + maxX) / 2;
    float poreCenterY = (minY + maxY) / 2;
    
    TH2F* hDensityXY = new TH2F("hDensityXY", "XY Plane Electron Density (Cumulative)",
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
    
    TH2F* hDensityXZ = new TH2F("hDensityXZ", "XZ Plane Electron Density (Cumulative)",
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
    
    // 시간 범위 계산
    float minTime = std::numeric_limits<float>::max();
    float maxTime = -std::numeric_limits<float>::max();
    
    // 트랙 정보를 사용하여 시간 범위 계산
    for (int i = 0; i < event.tracks.nTracks; i++) {
        minTime = std::min(minTime, event.tracks.birthTime[i]);
        maxTime = std::max(maxTime, event.tracks.finalTime[i]);
    }
    
    std::cout << "Time range: " << minTime << " - " << maxTime << " ps" << std::endl;
    
    // 트랙 정보를 히스토그램에 추가 (시작점과 종료점)
    for (int i = 0; i < event.tracks.nTracks; i++) {
        // 시작점
        float x = event.tracks.birthPosX[i];
        float y = event.tracks.birthPosY[i];
        float z = event.tracks.birthPosZ[i];
        
        if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ) {
            hDensityXY->Fill(x, y);
            hDensityXZ->Fill(x, z);
        }
        
        // 종료점
        x = event.tracks.finalPosX[i];
        y = event.tracks.finalPosY[i];
        z = event.tracks.finalPosZ[i];
        
        if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ) {
            hDensityXY->Fill(x, y);
            hDensityXZ->Fill(x, z);
        }
    }
    
    // 스텝 정보도 히스토그램에 추가 (가중치 낮게)
    for (int i = 0; i < event.steps.nSteps; i++) {
        float x = event.steps.posX[i];
        float y = event.steps.posY[i];
        float z = event.steps.posZ[i];
        
        if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ) {
            // 스텝도 동일한 가중치로 설정 (1.0)
            hDensityXY->Fill(x, y);
            hDensityXZ->Fill(x, z);
        }
    }
    
    // Draw XY plane density graph (separate canvas)
    TCanvas* canvasXY = new TCanvas("canvasDensityXY", "XY Plane Electron Density (Cumulative)", 1000, 800);
    canvasXY->SetLeftMargin(0.15);  
    canvasXY->SetRightMargin(0.15);  
    canvasXY->SetBottomMargin(0.15); 
    canvasXY->SetTopMargin(0.1);     

    gStyle->SetPalette(kViridis);
    
    hDensityXY->Draw("COLZ");
    
    // Draw pore boundary (XY plane)
    drawPoreBoundary(canvasXY, poreCenterX, poreCenterY, minX, maxX, minZ, maxZ);
    
    // 시간 범위 정보 추가
    TPaveText* timeInfo = new TPaveText(0.15, 0.85, 0.55, 0.90, "NDC");
    timeInfo->SetFillColor(0);
    timeInfo->SetTextAlign(12);
    timeInfo->SetTextSize(0.035);
    timeInfo->AddText(Form("Time range: %.2f - %.2f ps", minTime, maxTime));
    timeInfo->Draw();
    
    // 전자 수 정보 추가
    TPaveText* electronInfo = new TPaveText(0.15, 0.80, 0.55, 0.85, "NDC");
    electronInfo->SetFillColor(0);
    electronInfo->SetTextAlign(12);
    electronInfo->SetTextSize(0.035);
    electronInfo->AddText(Form("Total electrons: %d, Anode hits: %d", 
                             event.tracks.nTracks, event.GetAnodeElectronCount()));
    electronInfo->Draw();
    
    // Save XY density plot
    canvasXY->SaveAs((outputFileName + "_density_xy.pdf").c_str());
    delete canvasXY;
    
    // Draw XZ plane density graph (separate canvas)
    TCanvas* canvasXZ = new TCanvas("canvasDensityXZ", "XZ Plane Electron Density (Cumulative)", 1000, 800);
    canvasXZ->SetLeftMargin(0.15);  
    canvasXZ->SetRightMargin(0.15);  
    canvasXZ->SetBottomMargin(0.15); 
    canvasXZ->SetTopMargin(0.1);     

    hDensityXZ->Draw("COLZ");
    
    // 시간 범위 정보 추가
    TPaveText* timeInfoXZ = new TPaveText(0.15, 0.85, 0.55, 0.90, "NDC");
    timeInfoXZ->SetFillColor(0);
    timeInfoXZ->SetTextAlign(12);
    timeInfoXZ->SetTextSize(0.035);
    timeInfoXZ->AddText(Form("Time range: %.2f - %.2f ps", minTime, maxTime));
    timeInfoXZ->Draw();
    
    // 전자 수 정보 추가
    TPaveText* electronInfoXZ = new TPaveText(0.15, 0.80, 0.55, 0.85, "NDC");
    electronInfoXZ->SetFillColor(0);
    electronInfoXZ->SetTextAlign(12);
    electronInfoXZ->SetTextSize(0.035);
    electronInfoXZ->AddText(Form("Total electrons: %d, Anode hits: %d", 
                               event.tracks.nTracks, event.GetAnodeElectronCount()));
    electronInfoXZ->Draw();
    
    canvasXZ->SaveAs((outputFileName + "_density_xz.pdf").c_str());

    delete canvasXZ;
    delete hDensityXY;
    delete hDensityXZ;
    delete timeInfo;
    delete timeInfoXZ;
    delete electronInfo;
    delete electronInfoXZ;
}

// Create single pore amplification process animation
void createCascadeAnimation(const mcp::Event& event, const std::string& outputFileName) {
    // Calculate data range
    float minX, maxX, minY, maxY, minZ, maxZ;
    calculateDataRange(event, minX, maxX, minY, maxY, minZ, maxZ);

    // Calculate pore center
    float poreCenterX = (minX + maxX) / 2;
    float poreCenterY = (minY + maxY) / 2;
    
    // Find time range
    float minTime = std::numeric_limits<float>::max();
    float maxTime = -std::numeric_limits<float>::max();
    float minEnergy = std::numeric_limits<float>::max();
    float maxEnergy = -std::numeric_limits<float>::max();
    
    // 트랙 정보를 사용하여 시간 범위 계산
    for (int i = 0; i < event.tracks.nTracks; i++) {
        minTime = std::min(minTime, event.tracks.birthTime[i]);
        maxTime = std::max(maxTime, event.tracks.finalTime[i]);
        minEnergy = std::min(minEnergy, event.tracks.birthEnergy[i]);
        maxEnergy = std::max(maxEnergy, event.tracks.finalEnergy[i]);
    }
    
    // 스텝 정보를 사용하여 시간 범위 계산
    for (int i = 0; i < event.steps.nSteps; i++) {
        minTime = std::min(minTime, event.steps.time[i]);
        maxTime = std::max(maxTime, event.steps.time[i]);
        minEnergy = std::min(minEnergy, event.steps.energy[i]);
        maxEnergy = std::max(maxEnergy, event.steps.energy[i]);
    }
    
    // Set time interval (reduce memory usage by decreasing frame count)
    const int nFrames = 15; 
    float timeStep = (maxTime - minTime) / nFrames;
    
    gStyle->SetPalette(kViridis);
    
    // Create images for each frame
    for (int frame = 0; frame < nFrames; frame++) {
        float currentTime = minTime + frame * timeStep;
        float timeWindow = 0.1f; // 타임 윈도우를 0.5에서 0.1로 줄임 (더 정확한 시간대의 전자만 표시)
        
        TCanvas* canvas = new TCanvas("canvasAnim", 
                                    Form("MCP Pore Electron Cascade (Time: %.2f ps)", currentTime), 
                                    1300, 1100);  
        canvas->SetLeftMargin(0.15);  
        canvas->SetRightMargin(0.05);  
        canvas->SetBottomMargin(0.15); 
        canvas->SetTopMargin(0.1);     
        canvas->Divide(2, 2, 0.01, 0.01);  
        
        // XY plane
        canvas->cd(1);
        gPad->SetLeftMargin(0.15);  
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.1);     
        
        TH2F* hXY = new TH2F("hXY", Form("XY View (Time: %.2f ps)", currentTime),
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
        drawPoreBoundary(canvas->cd(1), poreCenterX, poreCenterY, minX, maxX, minZ, maxZ);
        
        // XZ plane
        canvas->cd(2);
        gPad->SetLeftMargin(0.15);  
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.1);     
        
        TH2F* hXZ = new TH2F("hXZ", "XZ View",
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
        
        // ZY plane (axis position exchange)
        canvas->cd(3);
        gPad->SetLeftMargin(0.15);  
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.1);     
        
        TH2F* hYZ = new TH2F("hYZ", "ZY View",  
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
        
        // 3D view
        canvas->cd(4);
        gPad->SetLeftMargin(0.15);  
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.18);  
        gPad->SetTopMargin(0.1);     
        gPad->SetTheta(30);
        gPad->SetPhi(30);
        
        TH3F* h3D = new TH3F("h3D", "3D View",
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
        
        // 현재 시간 창에 있는 트랙 표시
        for (int i = 0; i < event.tracks.nTracks; i++) {
            float birthTime = event.tracks.birthTime[i];
            float finalTime = event.tracks.finalTime[i];
            
            // 현재 시간 창에 있는 트랙만 표시
            if ((birthTime <= currentTime + timeWindow && finalTime >= currentTime - timeWindow) ||
                (birthTime >= currentTime - timeWindow && birthTime <= currentTime + timeWindow) ||
                (finalTime >= currentTime - timeWindow && finalTime <= currentTime + timeWindow)) {
                
                // 시작점 또는 종료점 중 현재 시간에 가까운 점 선택
                float time, x, y, z, energy;
                if (std::abs(birthTime - currentTime) < std::abs(finalTime - currentTime)) {
                    // 시작점이 더 가까움
                    time = event.tracks.birthTime[i];
                    x = event.tracks.birthPosX[i];
                    y = event.tracks.birthPosY[i];
                    z = event.tracks.birthPosZ[i];
                    energy = event.tracks.birthEnergy[i];
                } else {
                    // 종료점이 더 가까움
                    time = event.tracks.finalTime[i];
                    x = event.tracks.finalPosX[i];
                    y = event.tracks.finalPosY[i];
                    z = event.tracks.finalPosZ[i];
                    energy = event.tracks.finalEnergy[i];
                }
                
                // 색상 설정
                Int_t color = getEnergyColor(energy, minEnergy, maxEnergy);
                
                // XY 평면에 표시
                canvas->cd(1);
                TMarker* mXY = new TMarker(x, y, 20);
                mXY->SetMarkerColor(color);
                mXY->SetMarkerSize(0.8);
                mXY->Draw();
                
                // XZ 평면에 표시
                canvas->cd(2);
                TMarker* mXZ = new TMarker(x, z, 20);
                mXZ->SetMarkerColor(color);
                mXZ->SetMarkerSize(0.8);
                mXZ->Draw();
                
                // ZY 평면에 표시
                canvas->cd(3);
                TMarker* mYZ = new TMarker(z, y, 20);
                mYZ->SetMarkerColor(color);
                mYZ->SetMarkerSize(0.8);
                mYZ->Draw();
                
                // 3D 뷰에 표시
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
        
        // 현재 시간 창에 있는 스텝 표시 (작은 점으로)
        for (int i = 0; i < event.steps.nSteps; i++) {
            float time = event.steps.time[i];
            
            // 현재 시간 창에 있는 스텝만 표시
            if (time >= currentTime - timeWindow && time <= currentTime + timeWindow) {
                float x = event.steps.posX[i];
                float y = event.steps.posY[i];
                float z = event.steps.posZ[i];
                float energy = event.steps.energy[i];
                
                // 색상 설정
                Int_t color = getEnergyColor(energy, minEnergy, maxEnergy);
                
                // XY 평면에 표시
                canvas->cd(1);
                TMarker* mXY = new TMarker(x, y, 20);
                mXY->SetMarkerColor(color);
                mXY->SetMarkerSize(0.8);  // 스텝은 작게 표시
                mXY->Draw();
                
                // XZ 평면에 표시
                canvas->cd(2);
                TMarker* mXZ = new TMarker(x, z, 20);
                mXZ->SetMarkerColor(color);
                mXZ->SetMarkerSize(0.8);
                mXZ->Draw();
                
                // ZY 평면에 표시
                canvas->cd(3);
                TMarker* mYZ = new TMarker(z, y, 20);
                mYZ->SetMarkerColor(color);
                mYZ->SetMarkerSize(0.8);
                mYZ->Draw();
                
                // 3D 뷰에 표시
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
        
        // 정보 표시
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

        std::string frameFileName = outputFileName + "_frame" + 
                                  std::to_string(frame) + ".png";
        canvas->SaveAs(frameFileName.c_str());
        
        delete hXY;
        delete hXZ;
        delete hYZ;
        delete h3D;
        delete canvas;
        delete info;
        delete stats;
    }
    
    // Create shell script for animation
    std::string scriptName = outputFileName + "_create_gif.sh";
    std::ofstream scriptFile(scriptName);
    scriptFile << "#!/bin/bash\n";
    scriptFile << "convert -delay 20 -loop 0 " << outputFileName << "_frame*.png " 
              << outputFileName << ".gif\n";
    scriptFile.close();
    
    std::string chmodCmd = "chmod +x " + scriptName;
    system(chmodCmd.c_str());
}
