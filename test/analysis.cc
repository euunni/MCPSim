#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


void analysis() {

    // Load files.
    TString path = "/home/jangh/MCPSim/v3_track/output/250604_validation";
    std::vector<std::string> fileLists;

    TSystemDirectory dir("name", path);
    TList* files = dir.GetListOfFiles();
    TIter next(files);
    TSystemFile* file;

    while ((file = (TSystemFile*)next())) {
        if (!file->IsDirectory()) {
            TString fullPath = path + "/" + file->GetName(); 
            fileLists.push_back(fullPath.Data());
            std::cout << "Load file: " << fullPath << std::endl;
        }
    }

    const double mass = 9.1;
    std::vector<std::vector<double>> pos;
    std::vector<std::vector<double>> vel;
    std::vector<double> time;

    for (int i = 0; i < fileLists.size(); i++) {

        // Load data.
        double x, y, z;
        double vx, vy, vz;
        double timeX, timeY, timeZ;

        std::ifstream in;
        in.open((fileLists.at(i)).c_str(), std::ios::in);

        std::string firstLine;
        std::getline(in, firstLine);

        while (true) {
            in >> x >> y >> z >> vx >> vy >> vz >> timeX >> timeY >> timeZ;

            if (!in.good())
                break;

            pos.push_back({x, y, z});
            vel.push_back({vx, vy, vz});
            time.push_back(timeX);
            // std::cout << x << " " << y << " " << z << std::endl;
        }   
        in.close();
    } // End of file loop

    std::cout << "\nAnalyze " << fileLists.size() << " events and " << pos.size() << " secondary electrons.\n" << std::endl;

    // Analyze the data.
    TH1F* vxHist = new TH1F("vx", "vx;Velocity [#sqrt{eV/m}];Counts", 100, -10., 30.);
    TH1F* vyHist = new TH1F("vy", "vy;Velocity [#sqrt{eV/m}];Counts", 100, -10., 30.);
    TH1F* vzHist = new TH1F("vz", "vz;Velocity [#sqrt{eV/m}];Counts", 100, -10., 30.);
    TH1F* energyHist = new TH1F("Energy", "Energy;eV;Counts", 100, 0., 2000.);
    TH1F* TimingHist = new TH1F("Timing", "Timing;ps;Counts", 100, 0., 800.);
    TH2F* positionHist = new TH2F("Position", "Position;z [#mum];y [#mum]", 100, -600., 600., 100, -600., 600.);

    TCanvas* c1 = new TCanvas("", "", 800, 600);

    double timeCut = 220.;
    double timeTolerance = 5.;

    for (int i = 0; i < pos.size(); i++) {

        double energy = 0.5 * mass * (pow(vel.at(i).at(0),2) + pow(vel.at(i).at(1),2) + pow(vel.at(i).at(2),2));
        // if (time.at(i) > timeCut - timeTolerance && time.at(i) < timeCut + timeTolerance) {
        //     energyHist->Fill(energy);
        // }
        energyHist->Fill(energy);
        vxHist->Fill(vel.at(i).at(0));
        vyHist->Fill(vel.at(i).at(1));
        vzHist->Fill(vel.at(i).at(2));
        energyHist->Fill(energy);
        TimingHist->Fill(time.at(i));
        positionHist->Fill(pos.at(i).at(2), pos.at(i).at(1)); // Fill (z, y) values (x: beam axis)
    } // End of electron loop
    
    vxHist->SetLineColor(kBlack);
    vxHist->SetLineWidth(2);
    vxHist->Draw("hist");
    c1->SaveAs("./plot/250604_validation/vx.pdf");

    c1->Update();
    vyHist->SetLineColor(kBlack);
    vyHist->SetLineWidth(2);
    vyHist->Draw("hist");
    c1->SaveAs("./plot/250604_validation/vy.pdf");

    c1->Update();
    vzHist->SetLineColor(kBlack);
    vzHist->SetLineWidth(2);
    vzHist->Draw("hist");
    c1->SaveAs("./plot/250604_validation/vz.pdf");

    c1->Update();
    energyHist->SetLineColor(kBlue);
    energyHist->SetLineWidth(2);
    energyHist->Draw("hist");
    c1->SaveAs("./plot/250604_validation/energy.pdf");

    c1->Update();
    TimingHist->SetLineColor(kRed);
    TimingHist->SetLineWidth(2);
    TimingHist->Draw("hist");
    c1->SaveAs("./plot/250604_validation/Timing.pdf");
    
    c1->Update();
    positionHist->Draw("colz");
    c1->SaveAs("./plot/250604_validation/position.pdf");
}
