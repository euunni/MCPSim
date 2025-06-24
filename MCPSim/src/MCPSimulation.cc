#include "MCPSimulation.h"
#include "MCPConfig.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <cmath>
#include <set>

static int g_nextTrackID = 1; // 전역 trackID는 1부터 시작

namespace MCPSim {

// Process types
const int PROCESS_BACKSCATTER = 0; // Backscattering
const int PROCESS_REDIFFUSED = 1;  // Rediffused electron
const int PROCESS_SECONDARY = 2;   // Secondary emission

Simulation::Simulation()
    : physics(std::make_unique<Physics>())
    , x0(Config::getInstance().get("x0"))
    , x1(Config::getInstance().get("x1"))
    , x2(Config::getInstance().get("x2"))
    , q(Config::getInstance().get("q"))
    , limite(Config::getInstance().get("limite"))
    , I_strip(Config::getInstance().get("I_strip"))
    , c_c(Config::getInstance().get("c_c"))
    , m(Config::getInstance().get("m"))
    , nextTrackID_(0)
{
    Physics::initialize();  // Initialize static members
    tracks_.Reset();        // 트랙 초기화
    steps_.Reset();         // 스텝 초기화
}

// 새 전자 생성 (트랙 추가)
int Simulation::CreateElectron(int parentID, float time,
                             float posX, float posY, float posZ,
                             float velX, float velY, float velZ,
                             float energy, int procType) {
    // 외부 trackID 생성
    int externalTrackID = g_nextTrackID++;
    
    // 내부 trackID 생성 및 매핑
    int internalTrackID = tracks_.AddTrack(parentID, time, 
                                         posX, posY, posZ, 
                                         velX, velY, velZ, 
                                         energy, procType);
    
    // 매핑 저장
    trackIDMap_[externalTrackID] = internalTrackID;
    
    return externalTrackID;
}

// 전자 위치 업데이트 (스텝 추가)
void Simulation::AddElectronStep(int trackID, float time,
                               float posX, float posY, float posZ,
                               float velX, float velY, float velZ,
                               float energy, bool isInteraction) {
    // 외부 trackID에 해당하는 내부 trackID 찾기
    auto it = trackIDMap_.find(trackID);
    if (it != trackIDMap_.end()) {
        int internalTrackID = it->second;
        
        // 스텝 추가
        steps_.AddStep(internalTrackID, time, posX, posY, posZ, velX, velY, velZ, energy, isInteraction);
    }
}

// 전자 종료 (애노드 도달 또는 종료)
void Simulation::FinalizeElectron(int trackID, int status, float time,
                                float posX, float posY, float posZ,
                                float velX, float velY, float velZ,
                                float energy) {
    // 외부 trackID에 해당하는 내부 trackID 찾기
    auto it = trackIDMap_.find(trackID);
    if (it != trackIDMap_.end()) {
        int internalTrackID = it->second;
        
        // 트랙 종료 처리
        tracks_.FinalizeTrack(internalTrackID, status, time, posX, posY, posZ, velX, velY, velZ, energy);
    }
}

// Record the electron's position going to the anode
void Simulation::TrackElectronOutsidePore(const Matrix3x3& start, const Matrix3x3& end, int trackID, double cts) {
    const float timeStep = 1.0f; // Record at 1.0 ps intervals (changed from 0.5 ps)
    
    float startTime = start(2, 0);
    float endTime = end(2, 0);
    float totalTime = endTime - startTime;
    
    // 총 시간이 너무 짧으면 중간 지점 하나만 기록
    if (totalTime < timeStep) {
        // 중간 지점 하나만 기록
        float x = (start(0, 0) + end(0, 0)) * 0.5f;
        float y = (start(0, 1) + end(0, 1)) * 0.5f;
        float z = (start(0, 2) + end(0, 2)) * 0.5f;
        float time = (startTime + endTime) * 0.5f;
        float vx = (start(1, 0) + end(1, 0)) * 0.5f;
        float vy = (start(1, 1) + end(1, 1)) * 0.5f;
        float vz = (start(1, 2) + end(1, 2)) * 0.5f;
        
        // 에너지 계산
        float energy = 0.5f * m * (vx*vx + vy*vy + vz*vz);
        
        // 중간 위치에서의 전자 정보 기록 (스텝)
        AddElectronStep(trackID, time, x, y, z, vx, vy, vz, energy);
        return;
    }
    
    // 중간 단계 수 계산 (메모리 사용량 감소를 위해 단계 수 크게 제한)
    int numSteps = std::min(3, std::max(1, static_cast<int>(totalTime / timeStep)));
    
    for (int step = 1; step < numSteps; step++) {
        float ratio = static_cast<float>(step) / numSteps;
        float t = totalTime * ratio;
        
        // 중간 위치 계산 (물리 방정식 사용)
        float x = start(0, 0) + start(1, 0) * t + 0.5f * cts * t * t;
        float y = start(0, 1) + start(1, 1) * t;
        float z = start(0, 2) + start(1, 2) * t;
        
        // 속도 계산
        float vx = start(1, 0) + cts * t;
        float vy = start(1, 1);
        float vz = start(1, 2);
        
        // 시간 계산
        float time = startTime + t;
        
        // 에너지 계산
        float energy = 0.5f * m * (vx*vx + vy*vy + vz*vz);
        
        // 중간 위치에서의 전자 정보 기록 (스텝)
        AddElectronStep(trackID, time, x, y, z, vx, vy, vz, energy);
    }
}

std::vector<Matrix3x3> Simulation::Run(double initial_energy) {
    auto& config = Config::getInstance();
    double cts = config.get("c_c");

    std::vector<Matrix3x3> Resultat_final;
    
    std::cout << "\nSimulation steps:" << std::endl;
    std::cout << "  1. Generating photoelectron..." << std::endl;
    Matrix3x3 Mat = physics->Pho_ele(initial_energy);  // Generate photoelectron
    
    std::cout << "  2. Moving to MCP entrance..." << std::endl;
    Mat = physics->premiere_arrivee(Mat);  // Photoelectron reaches MCP entrance
    std::cout << "\nAfter premiere_arrivee: x=" << Mat(0,0) << ", y=" << Mat(0,1) << ", z=" << Mat(0,2) << std::endl;
    
    std::cout << "  3. Checking channel entry..." << std::endl;
    auto check = physics->Check_if_hit(Mat);  // Check if the photoelectron enters the channel
    std::cout << "\nChannel entry result: hit=" << check.first << ", channel=" << check.second << std::endl;
    std::cout << "y-position after arrival: " << Mat(0,1) << std::endl;
    
    if (check.first && -config.get("limite") < Mat(0, 1) && Mat(0, 1) < config.get("limite")) {
        std::cout << "  4. Electron entered valid channel region" << std::endl;
        int condition = 0;
        std::vector<Matrix3x3> Ensemble_emi;  // List of electrons that will collide
        
        // Generate 3 initial electrons
        std::cout << "  5. Generating initial 3 electrons..." << std::endl;
        for (int w = 0; w < 3; w++) {
            Matrix3x3 M = Matrix3x3::Zero();
            M(0, 0) = Mat(0, 0) - 0.1 * w;
            M(0, 1) = Mat(0, 1);
            M(0, 2) = Mat(0, 2);
            M(1, 0) = Mat(1, 0);
            M(1, 1) = Mat(1, 1);
            M(1, 2) = Mat(1, 2);
            M(2, 0) = Mat(2, 0);    // current time

            float vx = M(1,0), vy = M(1,1), vz = M(1,2);
            float E  = 0.5f * config.get("m") * (vx*vx + vy*vy + vz*vz);

            // 새 전자 생성 (트랙 추가)
            int trk = CreateElectron(-1, M(2,0),       // parent, time
                                M(0,0), M(0,1), M(0,2),  // pos
                                vx, vy, vz, E,            // vel, energy
                                PROCESS_SECONDARY);    // Initial electron is set to secondary
            
            // 트랙 ID를 행렬에 저장
            M(2, 2) = static_cast<float>(trk);

            double time_1_hit = physics->Point_de_contact2(M, check.second, cts);
            
            // Store the hit time in the matrix
            M(2, 1) = time_1_hit;

            Ensemble_emi.push_back(M);

            std::cout << "\nInitial position (x, y, z) = " << M(0, 0) << ", " << M(0, 1) << ", " << M(0, 2) << std::endl;
            std::cout << "Initial velocity (vx, vy, vz) = " << M(1, 0) << ", " << M(1, 1) << ", " << M(1, 2) << std::endl;
            std::cout << "Electron " << w+1 << " generated with hit time: " << time_1_hit << std::endl;
        }
        
        std::vector<Matrix3x3> Ensemble_non_emi;  // List of electrons that will not collide
        int Etat = 0;  // State of electric field adjustment
        double instant_0 = 0;
        
        std::cout << "\n############################################" << std::endl;
        std::cout << "     Starting main simulation loop..." << std::endl;
        std::cout << "############################################" << std::endl;
        int iteration = 0;
        while (condition == 0) {
            iteration++;
            
            if (!Ensemble_emi.empty()) {
                // Process the electron with the earliest collision
                double I = 0;
                Matrix3x3 lead_M = Ensemble_emi.front();
                double temps = lead_M(2, 1);
                double instant = lead_M(2, 0);
                
                // Move the electron to the collision point
                Matrix3x3 temp = physics->Transporter1(lead_M, temps, cts);
                
                // Preserve trackID
                int parentTrackID = static_cast<int>(lead_M(2, 2));
                temp(2, 2) = parentTrackID;
                
                // Secondary electron emission
                std::vector<ElectronProcess> Elec_secondaire = physics->emi_sec(temp, check.second);
                Ensemble_emi.erase(Ensemble_emi.begin());
                
                // Update non-colliding electrons
                if (!Ensemble_non_emi.empty()) {
                    std::vector<int> indice;
                    for (int j = 0; j < Ensemble_non_emi.size(); j++) {
                        // Modify to work similarly to the original code
                        Matrix3x3 M_passage = physics->Transporter2(Ensemble_non_emi[j], temps, cts);
                        
                        // Preserve trackID
                        int currentTrackID = static_cast<int>(Ensemble_non_emi[j](2, 2));
                        M_passage(2, 2) = currentTrackID;
                        
                        // 포어 밖에서 전자의 위치 기록
                        TrackElectronOutsidePore(Ensemble_non_emi[j], M_passage, currentTrackID, cts);
                        
                        if (M_passage(0, 0) >= config.get("x2")) {
                            // Process the electron that reaches the anode
                            Matrix3x3 M_recuperated = physics->Recuperation(M_passage);
                            
                            // Preserve trackID
                            M_recuperated(2, 2) = currentTrackID;
                            
                            // 애노드에 도달하는 과정의 위치 기록
                            TrackElectronOutsidePore(M_passage, M_recuperated, currentTrackID, cts);
                            
                            // 애노드 도달 처리
                            float finalTime = M_recuperated(2, 0);
                            float finalPosX = M_recuperated(0, 0);
                            float finalPosY = M_recuperated(0, 1);
                            float finalPosZ = M_recuperated(0, 2);
                            float finalVelX = M_recuperated(1, 0);
                            float finalVelY = M_recuperated(1, 1);
                            float finalVelZ = M_recuperated(1, 2);
                            float finalEnergy = 0.5f * m * (finalVelX*finalVelX + finalVelY*finalVelY + finalVelZ*finalVelZ);
                            
                            // 트랙 종료 (애노드 도달)
                            FinalizeElectron(currentTrackID, 1, finalTime,
                                          finalPosX, finalPosY, finalPosZ,
                                          finalVelX, finalVelY, finalVelZ,
                                          finalEnergy);
                            
                            Resultat_final.push_back(M_recuperated);
                            indice.push_back(j);
                        } else {
                            Ensemble_non_emi[j] = M_passage;
                            if (Ensemble_non_emi[j](0, 0) < config.get("x1")) {
                                I += config.get("q") * Ensemble_non_emi[j](1, 0) / (config.get("x1") - config.get("x0"));
                            }
                        }
                    }
                    
                    std::sort(indice.begin(), indice.end(), std::greater<int>());
                    for (int idx : indice) {
                        Ensemble_non_emi.erase(Ensemble_non_emi.begin() + idx);
                    }
                }
                
                // Update electrons that are waiting for collision
                for (int j = 0; j < Ensemble_emi.size(); j++) {
                    // Modify to work similarly to the original code
                    Matrix3x3 M_passage = physics->Transporter1(Ensemble_emi[j], temps, cts);
                    
                    // Preserve trackID
                    int currentTrackID = static_cast<int>(Ensemble_emi[j](2, 2));
                    M_passage(2, 2) = currentTrackID;
                    
                    Ensemble_emi[j] = M_passage;
                    
                    if (Ensemble_emi[j](0, 0) < config.get("x1")) {
                        I += config.get("q") * Ensemble_emi[j](1, 0) / 
                             (config.get("x1") - config.get("x0"));
                    }
                }
                
                // Process new secondary electrons
                for (auto& elec : Elec_secondaire) {
                    // Set time
                    if (elec.matrix(2,0) == 0) {
                        elec.matrix(2,0) = instant; // Set current time
                    }
                    
                    // Kinematic parameters
                    float vx_s = elec.matrix(1,0), vy_s = elec.matrix(1,1), vz_s = elec.matrix(1,2);
                    float en_s = 0.5f * config.get("m") * (vx_s*vx_s + vy_s*vy_s + vz_s*vz_s);
                    
                    // 새 전자 생성 (트랙 추가)
                    int newTrk = CreateElectron(parentTrackID, elec.matrix(2,0),
                                             elec.matrix(0,0), elec.matrix(0,1), elec.matrix(0,2),
                                             vx_s, vy_s, vz_s,
                                             en_s, elec.processType);
                    
                    // 트랙 ID를 행렬에 저장
                    elec.matrix(2,2) = static_cast<float>(newTrk);
                    
                    // Calculate collision time
                    double time = physics->Point_de_contact2(elec.matrix, check.second, cts);
                    
                    if (time == false) {  // MCP exits
                        Ensemble_non_emi.push_back(elec.matrix);
                    } else if (time == true) {  // No collision
                        continue;
                    } else {  // Collision expected
                        Matrix3x3 temp = elec.matrix;
                        temp(2, 1) = time;
                        physics->ajouter_element_trie(Ensemble_emi, temp);
                    }
                }
                
                // Dynamic adjustment of electric field
                if (Etat == 0) {
                    if (I >= 0.05 * config.get("I_strip")) {
                        instant_0 = instant;
                        cts = cts * 0.8;
                        auto result = physics->Rearrangement(Ensemble_emi, Ensemble_non_emi, check.second, cts);
                        Ensemble_emi = result.first;
                        Ensemble_non_emi = result.second;
                        Etat = 1;
                    }
                } else {
                    if (instant >= instant_0 + 5 && I < 0.05 * config.get("I_strip")) {
                        instant_0 = instant;
                        cts = cts / 0.8;
                        auto result = physics->Rearrangement(Ensemble_emi, Ensemble_non_emi, check.second, cts);
                        Ensemble_emi = result.first;
                        Ensemble_non_emi = result.second;
                        
                        if (std::abs(cts - config.get("c_c")) < 1e-10) {
                            Etat = 0;
                        }
                    }
                    
                    if (instant >= instant_0 + 5 && I > 0.05 * config.get("I_strip")) {
                        cts = cts * 0.8;
                        instant_0 = instant;
                        auto result = physics->Rearrangement(Ensemble_emi, Ensemble_non_emi, check.second, cts);
                        Ensemble_emi = result.first;
                        Ensemble_non_emi = result.second;
                    }
                }
            } else {
                // Process remaining non-colliding electrons
                std::vector<int> indice2;
                for (int j = 0; j < Ensemble_non_emi.size(); j++) {
                    Matrix3x3 result = physics->Recuperation(Ensemble_non_emi[j]);
                    
                    // Preserve trackID
                    int trkID = static_cast<int>(Ensemble_non_emi[j](2, 2));
                    result(2, 2) = trkID;
                    
                    // 애노드에 도달하는 과정의 위치 기록
                    TrackElectronOutsidePore(Ensemble_non_emi[j], result, trkID, cts);
                    
                    // 애노드 도달 처리
                    float finalTime = result(2, 0);
                    float finalPosX = result(0, 0);
                    float finalPosY = result(0, 1);
                    float finalPosZ = result(0, 2);
                    float finalVelX = result(1, 0);
                    float finalVelY = result(1, 1);
                    float finalVelZ = result(1, 2);
                    float finalEnergy = 0.5f * m * (finalVelX*finalVelX + finalVelY*finalVelY + finalVelZ*finalVelZ);
                    
                    // 트랙 종료 (애노드 도달)
                    FinalizeElectron(trkID, 1, finalTime,
                                  finalPosX, finalPosY, finalPosZ,
                                  finalVelX, finalVelY, finalVelZ,
                                  finalEnergy);
                    
                    Resultat_final.push_back(result);
                    indice2.push_back(j);
                }

                std::sort(indice2.begin(), indice2.end(), std::greater<int>());
                for (int idx : indice2) {
                    Ensemble_non_emi.erase(Ensemble_non_emi.begin() + idx);
                }
                
                if (Ensemble_emi.empty() && Ensemble_non_emi.empty()) {
                    condition = 1;
                }
            }
        }
    } else {
        std::cout << "Electron did not enter valid channel region" << std::endl;
    }
    
    std::cout << "\nSimulation completed with " << Resultat_final.size() 
              << " electrons reaching the anode" << std::endl;
    
    return Resultat_final;
}

void Simulation::Save(const std::vector<Matrix3x3>& results, double initial_energy, const std::string& rootfile) {
    MCPRootManager rootMgr;
    if (!rootMgr.OpenFile(rootfile)) {
        throw std::runtime_error("Could not create ROOT file: " + rootfile);
    }

    mcp::Event evt = ConvertEvent(results, initial_energy);
    rootMgr.WriteEvt(evt);
    rootMgr.Close();

    std::cout << "Results saved to: " << rootfile << std::endl;
}

mcp::Event Simulation::ConvertEvent(const std::vector<Matrix3x3>& results, double initial_energy) {
    mcp::Event evt;
    
    // Event info
    evt.eventInfo.initialEnergy = initial_energy;
    
    // Config
    evt.config.LoadFromConfig();
    
    evt.tracks = tracks_;
    evt.steps = steps_;
    
    int anodeElectronCount = evt.GetAnodeElectronCount();
    std::cout << "\nTotal number of electrons: " << evt.tracks.nTracks << std::endl;
    std::cout << "Number of electrons reaching the anode: " << anodeElectronCount << std::endl;
    std::cout << "Number of trajectory steps recorded: " << evt.steps.nSteps << std::endl;
    
    return evt;
}

} // namespace MCPSim 