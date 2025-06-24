#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include "include/MCPAnalyzer.h"
#include "include/MCPVisualizer.h"


int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <file> <doAnalysis> <doVisualization> [options]" << std::endl;
        return 1;
    }
    
    std::string inputFileName = argv[1];
    bool doAnalysis = (std::stoi(argv[2]) != 0);
    bool doVisualization = (std::stoi(argv[3]) != 0);
    
    bool anodeOnly = true;  // Default: anode only
    int eventIndex = 0;
    std::string outputPrefix = "output";
    
    // Parse optional arguments
    for (int i = 4; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--allEle") {
            anodeOnly = false;
        } else if (arg.find("--event=") == 0) {
            eventIndex = std::stoi(arg.substr(8));
        } else if (arg.find("--output=") == 0) {
            outputPrefix = arg.substr(9);
        } else {
            std::cout << "Unknown option: " << arg << std::endl;
            return 1;
        }
    }
    
    // Validate arguments
    if (!doAnalysis && !doVisualization) {
        std::cout << "Error: At least one of doAnalysis or doVisualization must be 1" << std::endl;
        return 1;
    }
    
    // Create output directory if needed
    std::string outputDir = "../plots/";
    if (gSystem->AccessPathName(outputDir.c_str())) {
        gSystem->mkdir(outputDir.c_str(), kTRUE);
    }
    
    // Initialize analyzer and load data
    MCPAnalyzer analyzer;
    if (!analyzer.Load(inputFileName)) {
        std::cerr << "Failed to load data from: " << inputFileName << std::endl;
        return 1;
    }
    
    // Select event
    analyzer.SelectEvent(eventIndex);
    
    // Print basic information
    std::cout << "Total electrons: " << analyzer.GetEleCount() << std::endl;
    std::cout << "Electrons reaching anode: " << analyzer.GetEleCount(true) << std::endl;
    
    if (anodeOnly) {
        outputPrefix += "_anode";
    } else {
        outputPrefix += "_all";
    }
    
    // Analysis
    if (doAnalysis) {
        std::cout << "\nStarting analysis..." << std::endl;
        
        // Energy analysis
        std::vector<float> energies = analyzer.GetEnergy(anodeOnly);

        TCanvas* c = new TCanvas("c_energy", "", 800, 600);
        TH1F* energyHist = new TH1F("energy", "Energy;Energy [eV];Count", 100, 0, 2000);

        for (float e : energies) {
            energyHist->Fill(e);
        }
        
        energyHist->SetLineColor(kBlue);
        energyHist->SetLineWidth(2);
        energyHist->Draw("hist");
        c->SaveAs((outputDir + outputPrefix + "_energy.pdf").c_str());
        delete c;
        delete energyHist;
        
        // Time analysis
        std::vector<float> times = analyzer.GetTime(anodeOnly);
        std::vector<float> transitTimes = analyzer.GetTransitTime(anodeOnly);

        // Regular time histogram
        TCanvas* c1 = new TCanvas("c_time", "", 800, 600);
        
        // Calculate median for time range
        std::vector<float> timesSorted = times;
        std::sort(timesSorted.begin(), timesSorted.end());
        float timeMedian = timesSorted[timesSorted.size() / 2];
        float timeMin = timeMedian - 250.0f;
        float timeMax = timeMedian + 250.0f;
        
        TH1F* timeHist = new TH1F("time", "Time;Time [ps];Count", 100, timeMin, timeMax);
        
        for (float t : times) {
            timeHist->Fill(t);
        }
        
        c1->cd();
        timeHist->SetLineColor(kRed);
        timeHist->SetLineWidth(2);
        timeHist->Draw("hist");
        c1->SaveAs((outputDir + outputPrefix + "_time.pdf").c_str());
        delete c1;
        delete timeHist;
        
        // Transit time histogram
        TCanvas* c2 = new TCanvas("c_transit", "", 800, 600);
        
        // Calculate median for transit time range
        std::vector<float> transitTimesSorted = transitTimes;
        std::sort(transitTimesSorted.begin(), transitTimesSorted.end());
        float transitMedian = transitTimesSorted[transitTimesSorted.size() / 2];
        float transitMin = transitMedian - 250.0f;
        float transitMax = transitMedian + 250.0f;
        
        TH1F* transitHist = new TH1F("transit", "Transit Time;Transit Time [ps];Count", 100, transitMin, transitMax);
        
        for (float t : transitTimes) {
            transitHist->Fill(t);
        }
        
        c2->cd();
        transitHist->SetLineColor(kRed);
        transitHist->SetLineWidth(2);
        transitHist->Draw("hist");
        // c2->SaveAs((outputDir + outputPrefix + "_transit_time.pdf").c_str());
        delete c2;
        delete transitHist;

        // Position analysis
        std::vector<float> posX = analyzer.GetPosX(anodeOnly);
        std::vector<float> posY = analyzer.GetPosY(anodeOnly);
        std::vector<float> posZ = analyzer.GetPosZ(anodeOnly);

        TCanvas* c3 = new TCanvas("c_position", "", 800, 600);
        TH2F* posHist = new TH2F("Position", "Position;z [#mum];y [#mum]", 
                                100, -(*std::max_element(posZ.begin(), posZ.end()) * 1.1), 
                                (*std::max_element(posZ.begin(), posZ.end()) * 1.1), 
                                100, -(*std::max_element(posY.begin(), posY.end()) * 1.1), 
                                (*std::max_element(posY.begin(), posY.end()) * 1.1));

        for (int i = 0; i < posX.size(); i++) {
            posHist->Fill(posZ[i], posY[i]);
        }

        posHist->Draw("colz");
        c3->SaveAs((outputDir + outputPrefix + "_position.pdf").c_str());
        delete c3;
        delete posHist;
    }
    
    // Visualization (Trajectory, Density, Animation)
    if (doVisualization) {
        std::cout << "\nStarting visualization..." << std::endl;
        
        MCPVisualizer visualizer(&analyzer);

        // Trajectory visualization - separate XY and XZ
        TCanvas* trajCanvasXY = visualizer.DrawTrajectoryXY();
        trajCanvasXY->SaveAs((outputDir + outputPrefix + "_trajectory_xy.pdf").c_str());
        delete trajCanvasXY;
        
        TCanvas* trajCanvasXZ = visualizer.DrawTrajectoryXZ();
        trajCanvasXZ->SaveAs((outputDir + outputPrefix + "_trajectory_xz.pdf").c_str());
        delete trajCanvasXZ;
        
        // Draw XY density
        TCanvas* densityXY = visualizer.DrawDensityXY();
        densityXY->SaveAs((outputDir + outputPrefix + "_density_xy.pdf").c_str());
        delete densityXY;
        
        // Draw XZ density
        TCanvas* densityXZ = visualizer.DrawDensityXZ();
        densityXZ->SaveAs((outputDir + outputPrefix + "_density_xz.pdf").c_str());
        delete densityXZ;
        
        // Animation
        int frameCount = visualizer.GetAnimationFrameCount();
        
        // Generate each frame
        for (int i = 0; i < frameCount; i++) {
            TCanvas* frame = visualizer.AnimateCascadeFrame(i, frameCount);
            std::string frameFileName = outputDir + outputPrefix + "_animation_frame" + std::to_string(i) + ".png";
            frame->SaveAs(frameFileName.c_str());
            delete frame;
        }
        
        // Create GIF creation script
        std::string scriptName = outputDir + outputPrefix + "_animation_create_gif.sh";
        std::ofstream scriptFile(scriptName);
        scriptFile << "#!/bin/bash\n";
        scriptFile << "convert -delay 20 -loop 0 " << outputDir + outputPrefix + "_animation_frame*.png " 
                   << outputDir + outputPrefix + "_animation.gif\n";
        scriptFile.close();
        
        std::string chmodCmd = "chmod +x " + scriptName;
        system(chmodCmd.c_str());
    }
    
    std::cout << "Analysis completed successfully." << std::endl;

    return 0;
}
