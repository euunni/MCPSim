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

// --- 전자 상태 enum 및 구조체 추가 ---
enum ElectronState { MCP1, GAP1, MCP2, GAP2, ANODE, DEAD };
struct Electron {
    Matrix3x3 state;
    ElectronState status;
    int trackID;
    int parentID;
    int processType;
    int channel; // pore index n that the electron belongs to
};

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
    tracks_.Reset();        // Track initialization
    steps_.Reset();         // Step initialization
}

// Create new electron (add track)
int Simulation::CreateElectron(int parentID, float time,
                             float posX, float posY, float posZ,
                             float velX, float velY, float velZ,
                             float energy, int procType) {
    // Create external trackID
    int externalTrackID = g_nextTrackID++;
    
    // Create internal trackID and mapping
    int internalTrackID = tracks_.AddTrack(parentID, time, 
                                         posX, posY, posZ, 
                                         velX, velY, velZ, 
                                         energy, procType);
    
    // Save mapping
    trackIDMap_[externalTrackID] = internalTrackID;
    
    return externalTrackID;
}

// Update electron position (add step)
void Simulation::AddElectronStep(int trackID, float time,
                               float posX, float posY, float posZ,
                               float velX, float velY, float velZ,
                               float energy) {
    // Find internal trackID corresponding to external trackID
    auto it = trackIDMap_.find(trackID);
    if (it != trackIDMap_.end()) {
        int internalTrackID = it->second;
        
        // Add step
        steps_.AddStep(internalTrackID, time, posX, posY, posZ, velX, velY, velZ, energy);
    }
}

// Electron termination (anode hit or end)
void Simulation::FinalizeElectron(int trackID, int status, float time,
                                float posX, float posY, float posZ,
                                float velX, float velY, float velZ,
                                float energy) {
    // Find internal trackID corresponding to external trackID
    auto it = trackIDMap_.find(trackID);
    if (it != trackIDMap_.end()) {
        int internalTrackID = it->second;
        
        // Finalize track
        tracks_.FinalizeTrack(internalTrackID, status, time, posX, posY, posZ, velX, velY, velZ, energy);
    }
}

// Record the electron's position going to the anode
void Simulation::TrackElectronOutsidePore(const Matrix3x3& start, const Matrix3x3& end, int trackID, double cts) {
    const float timeStep = 1.0f;
    
    float startTime = start(2, 0);
    float endTime = end(2, 0);
    float totalTime = endTime - startTime;
    
    // If totalTime is too short, record only one intermediate point
    if (totalTime < timeStep) {
        // Record only one intermediate point
        float x = (start(0, 0) + end(0, 0)) * 0.5f;
        float y = (start(0, 1) + end(0, 1)) * 0.5f;
        float z = (start(0, 2) + end(0, 2)) * 0.5f;
        float time = (startTime + endTime) * 0.5f;
        float vx = (start(1, 0) + end(1, 0)) * 0.5f;
        float vy = (start(1, 1) + end(1, 1)) * 0.5f;
        float vz = (start(1, 2) + end(1, 2)) * 0.5f;
        
        // Calculate energy
        float energy = 0.5f * m * (vx*vx + vy*vy + vz*vz);
        
        // Record electron information at the intermediate position (step)
        AddElectronStep(trackID, time, x, y, z, vx, vy, vz, energy);
        return;
    }
    
    // Calculate the number of intermediate steps (limit the number of steps to reduce memory usage)
    int numSteps = std::min(3, std::max(1, static_cast<int>(totalTime / timeStep)));
    
    for (int step = 1; step < numSteps; step++) {
        float ratio = static_cast<float>(step) / numSteps;
        float t = totalTime * ratio;
        
        // Calculate the intermediate position (assuming uniform acceleration along X
        // and a possible projection of that acceleration on Y due to the channel tilt).
        float x = start(0, 0) + start(1, 0) * t + 0.5f * cts * t * t;

        float y = start(0, 1) + start(1, 1) * t;

        float z = start(0, 2) + start(1, 2) * t;
        
        // Calculate velocity
        float vx = start(1, 0) + cts * t;
        float vy = start(1, 1);
        float vz = start(1, 2);
        
        // Calculate time
        float time = startTime + t;
        
        // Calculate energy
        float energy = 0.5f * m * (vx*vx + vy*vy + vz*vz);
        
        // Record electron information at the intermediate position (step)
        AddElectronStep(trackID, time, x, y, z, vx, vy, vz, energy);
    }
}

std::vector<Matrix3x3> Simulation::Run(double initial_energy) {

    int num_anode = 100;

    auto& config = Config::getInstance();
    // 파라미터 읽기
    double x0 = config.get("x0");
    double x1 = config.get("x1");
    double x2 = config.get("x2");
    double x3 = config.get("x3");
    double x4 = config.get("x4");
    double alpha1 = config.get("alpha1");
    double alpha2 = config.get("alpha2");
    double R = config.get("R");
    double dia = config.get("dia");
    double pas = config.get("pas");
    double c_c1 = config.get("c_c");
    double c_c2 = config.get("c_c2");
    double m = config.get("m");
    double E0 = config.get("E0");
    double gap_pot1 = config.get("gap_pot1");
    double gap_pot2 = config.get("gap_pot2");
    double q = config.get("q");
    double limite = config.get("limite");
    double I_strip = config.get("I_strip");
    // GAP 구간에서는 현재 MCP 내부에서 사용 중인 가속도 계수( c_c1 / c_c2 )를
    // 그대로 사용한다. 이 값들은 이후 전류(I)에 따라 동적으로 조정될 예정이므로
    // 별도의 고정 상수를 두지 않고 직접 참조한다.

    // -------------------------------
    // Dynamic electric-field coefficients (cts) that will be adjusted
    // according to the instantaneous channel current, following
    // v3_track behaviour.
    double cts1 = c_c1;   // used in MCP1 + GAP1
    double cts2 = c_c2;   // used in MCP2 + GAP2

    // Flags indicating whether the field has been reduced (1) or normal (0)
    int Etat1 = 0;
    int Etat2 = 0;
    // ------------------------------------------------

    std::vector<Matrix3x3> Resultat_final;
    std::vector<Electron> electrons;
    // 1. 광전자 생성 및 MCP1 입구 도달
    Matrix3x3 Mat = physics->Pho_ele(initial_energy);
    Mat = physics->premiere_arrivee(Mat);
    auto check = physics->Check_if_hit(Mat);
    if (!(check.first && -limite < Mat(0, 1) && Mat(0, 1) < limite)) {
        std::cout << "Electron did not enter valid channel region" << std::endl;
        return Resultat_final;
    }

    // Clear ensemble vectors (in case Run() called twice)
    E1_emi_.clear();  E1_nonemi_.clear();  G1_.clear();
    E2_emi_.clear();  E2_nonemi_.clear();  G2_.clear();

    // 2. MCP1 증폭 루프 (초기 3전자)
    for (int w = 0; w < 3; w++) {
        Matrix3x3 M = Matrix3x3::Zero();
        M(0, 0) = Mat(0, 0) - 0.1 * w;
        M(0, 1) = Mat(0, 1);
        M(0, 2) = Mat(0, 2);
        M(1, 0) = Mat(1, 0);
        M(1, 1) = Mat(1, 1);
        M(1, 2) = Mat(1, 2);
        M(2, 0) = Mat(2, 0);

        float vx = M(1,0), vy = M(1,1), vz = M(1,2);
        float E  = 0.5f * m * (vx*vx + vy*vy + vz*vz);
        int trk = CreateElectron(-1, M(2,0), M(0,0), M(0,1), M(0,2), vx, vy, vz, E, PROCESS_SECONDARY);
        M(2, 2) = static_cast<float>(trk);

        // v3_track 방식: 첫 번째 충돌까지 시간 계산
        double t_hit = physics->Point_de_contact2(M, check.second, cts1, alpha1, x0, x1, R, dia, pas);
        M(2, 1) = t_hit;

        // 기존 상태-머신용 electrons 벡터에도 그대로 추가 (단계적 이전을 위해 유지)
        electrons.push_back({M, MCP1, trk, -1, PROCESS_SECONDARY, check.second});

        // --- 새 Ensemble 벡터에도 동일하게 push ---
        if (t_hit == true) {
            // 불가능(벽 안맞고 뒤로 빠짐)이지만 호환 위해 무시
            continue;
        } else if (t_hit == false) {
            E1_nonemi_.push_back(M);
        } else {
            // t_hit > 0 → 벽 충돌 예정, (2,1)에 t 저장 후 시간순 삽입
            physics->ajouter_element_trie(E1_emi_, M);
        }
    }

    // --- 메인 루프: 모든 전자를 동시에 상태 기반으로 관리 ---
    bool all_finished = false; // enter loop at least once
    int iteration = 0;
    while (!all_finished) {
        all_finished = true; // assume finished; will set false if active electron found
        std::cout << "[DEBUG] iteration: " << iteration << std::endl;
        // 상태별 전자 개수 카운트
        int n_mcp1 = 0, n_gap1 = 0, n_mcp2 = 0, n_gap2 = 0, n_anode = 0, n_dead = 0;
        for (const auto& e : electrons) {
            switch (e.status) {
                case MCP1: n_mcp1++; break;
                case GAP1: n_gap1++; break;
                case MCP2: n_mcp2++; break;
                case GAP2: n_gap2++; break;
                case ANODE: n_anode++; break;
                case DEAD: n_dead++; break;
            }
        }
        std::cout << "[iter " << iteration << "] ANODE: " << n_anode << std::endl;

        // --- temporary: process MCP1 ensemble stub ---
        int dummyID = 0;
        ProcessMCP1Ensemble(cts1, alpha1, x0, x1, R, dia, pas, m, E0, dummyID);

        // -------------------------------------------------
        // Current-based electric field adjustment (B)
        // -------------------------------------------------
        double I1 = 0.0, I2 = 0.0;
        for (const auto& ec : electrons) {
            if (ec.status == MCP1 || ec.status == GAP1) {
                // Same denominator (channel length) as v3_track for continuity
                if (ec.state(0,0) < x2) { // up to MCP2 entrance
                    I1 += q * ec.state(1,0) / (x1 - x0);
                }
            } else if (ec.status == MCP2 || ec.status == GAP2) {
                if (ec.state(0,0) < x4) { // up to anode
                    I2 += q * ec.state(1,0) / (x3 - x2);
                }
            }
        }

        // MCP1 coefficient update
        if (Etat1 == 0) {
            if (I1 >= 0.05 * I_strip) {
                cts1 *= 0.8;
                Etat1 = 1;
            }
        } else {
            if (I1 < 0.05 * I_strip) {
                cts1 /= 0.8;
                if (std::abs(cts1 - c_c1) < 1e-12) {
                    Etat1 = 0;
                }
            }
        }

        // MCP2 coefficient update
        if (Etat2 == 0) {
            if (I2 >= 0.05 * I_strip) {
                cts2 *= 0.8;
                Etat2 = 1;
            }
        } else {
            if (I2 < 0.05 * I_strip) {
                cts2 /= 0.8;
                if (std::abs(cts2 - c_c2) < 1e-12) {
                    Etat2 = 0;
                }
            }
        }

        if (n_anode >= num_anode) {
            std::cout << "Early stopping: ANODE electrons reached " << num_anode << std::endl;
            break;
        }
        size_t n_elec = electrons.size();
        std::vector<Electron> new_electrons;
            iteration++;
        for (size_t i = 0; i < n_elec; ++i) {
            auto& e = electrons[i];
            if (e.status == DEAD || e.status == ANODE) continue;
            all_finished = false; // active electron exists
            switch (e.status) {
                case MCP1: {
                    int n_mcp1 = e.channel; // use stored channel index

                    double t_hit = physics->Point_de_contact2(e.state, n_mcp1, cts1, alpha1, x0, x1, R, dia, pas);
                    if (t_hit == false) {
                        // Electron exits MCP1 without hitting the wall – move it exactly to the pore exit (x1) first
                        Matrix3x3 M_exit = physics->RecuperationTo(e.state, x1);
                        TrackElectronOutsidePore(e.state, M_exit, e.trackID, cts1);
                        e.state = M_exit;
                        e.status = GAP1;
                    } else if (t_hit == true) {
                        e.status = DEAD;
                    } else if (std::isnan(t_hit) || t_hit <= 0) {
                        e.status = DEAD;
                        } else {
                        Matrix3x3 temp = physics->Transporter1(e.state, t_hit, cts1, alpha1);
                        temp(2, 2) = e.trackID;
                        AddElectronStep(e.trackID, temp(2,0), temp(0,0), temp(0,1), temp(0,2), temp(1,0), temp(1,1), temp(1,2), 0.5f * m * (temp(1,0)*temp(1,0) + temp(1,1)*temp(1,1) + temp(1,2)*temp(1,2)));
                        auto secondaries = physics->emi_sec(temp, n_mcp1, alpha1, x0, R, dia, pas, m, E0);
                        for (auto& sec : secondaries) {
                            float vx_s = sec.matrix(1,0), vy_s = sec.matrix(1,1), vz_s = sec.matrix(1,2);
                            float en_s = 0.5f * m * (vx_s*vx_s + vy_s*vy_s + vz_s*vz_s);
                            int newTrk = CreateElectron(e.trackID, sec.matrix(2,0), sec.matrix(0,0), sec.matrix(0,1), sec.matrix(0,2), vx_s, vy_s, vz_s, en_s, sec.processType);
                            sec.matrix(2,2) = static_cast<float>(newTrk);
                            double t_hit2 = physics->Point_de_contact2(sec.matrix, e.channel, cts1, alpha1, x0, x1, R, dia, pas);
                            ElectronState next_state;
                            Matrix3x3 state_for_push = sec.matrix;
                            if (t_hit2 == false) {
                                // Transport exactly to pore exit before entering GAP1
                                Matrix3x3 M_exit = physics->RecuperationTo(sec.matrix, x1);
                                TrackElectronOutsidePore(sec.matrix, M_exit, newTrk, cts1);
                                state_for_push = M_exit;
                                next_state = GAP1;
                            } else if (t_hit2 == true) {
                                next_state = DEAD;
                            } else if (std::isnan(t_hit2) || t_hit2 <= 0) {
                                next_state = DEAD;
                            } else {
                                next_state = MCP1;
                            }
                            new_electrons.push_back({state_for_push, next_state, newTrk, e.trackID, sec.processType, e.channel});
                        }
                        e.state = temp;
                    }
                    break;
                }
                case GAP1: {
                    const double dt_gap = 1.0; // 1 ps step like v3_track temps_min
                    double dt = dt_gap;
                    // Predict the x position after dt to check boundary crossing
                    double x_pred = e.state(0, 0) + e.state(1, 0) * dt + 0.5 * cts1 * dt * dt;
                    if (x_pred >= x2) {
                        // compute exact remaining time to reach x2
                        double a = cts1 / 2.0;
                        double b = e.state(1, 0);
                        double c_val = e.state(0, 0) - x2;
                        dt = physics->Resolution(a, b, c_val);
                        if (std::isnan(dt) || dt <= 0) {
                            std::cout << "[Warn] GAP1: invalid dt to boundary -> DEAD" << std::endl;
                            e.status = DEAD;
                            break;
                        }
                    }
                    // Move by dt
                    Matrix3x3 M_step = physics->Transporter2(e.state, dt, cts1, alpha1);
                    TrackElectronOutsidePore(e.state, M_step, e.trackID, cts1);
                    e.state = M_step;

                    // If boundary (x2) reached, attempt to enter MCP2
                    if (e.state(0, 0) >= x2 - 1e-6) {
                        // Snap exactly onto the entry plane
                        Matrix3x3 M_exact = physics->RecuperationTo(e.state, x2);
                        TrackElectronOutsidePore(e.state, M_exact, e.trackID, cts1);
                        e.state = M_exact;

                        // Check if the electron actually hits a pore opening of MCP2
                        auto check2 = physics->Check_if_hit(e.state);

                        if (!(check2.first && -limite < e.state(0, 1) && e.state(0, 1) < limite)) {
                            // Missed the pore entrance of MCP2 → electron is lost.
                            e.status = DEAD;
                        } else {
                            // Valid pore entry; store channel index and move to MCP2.
                            e.channel = check2.second;
                            e.status = MCP2;
                        }
                    }
                    break;
                }
                case MCP2: {
                    int n_mcp2 = e.channel; // stored channel index
                    double t_hit = physics->Point_de_contact2(e.state, n_mcp2, cts2, alpha2, x2, x3, R, dia, pas);
                    if (t_hit == false) {
                        // Electron exits MCP2 – bring it to pore exit (x3) before entering GAP2
                        Matrix3x3 M_exit = physics->RecuperationTo(e.state, x3);
                        TrackElectronOutsidePore(e.state, M_exit, e.trackID, cts2);
                        e.state = M_exit;
                        e.status = GAP2;
                    } else if (t_hit == true) {
                        e.status = DEAD;
                    } else if (std::isnan(t_hit) || t_hit <= 0) {
                        e.status = DEAD;
                    } else {
                        Matrix3x3 temp = physics->Transporter1(e.state, t_hit, cts2, alpha2);
                        temp(2, 2) = e.trackID;
                        AddElectronStep(e.trackID, temp(2,0), temp(0,0), temp(0,1), temp(0,2), temp(1,0), temp(1,1), temp(1,2), 0.5f * m * (temp(1,0)*temp(1,0) + temp(1,1)*temp(1,1) + temp(1,2)*temp(1,2)));
                        auto secondaries = physics->emi_sec(temp, n_mcp2, alpha2, x2, R, dia, pas, m, E0);
                        for (auto& sec : secondaries) {
                            float vx_s = sec.matrix(1,0), vy_s = sec.matrix(1,1), vz_s = sec.matrix(1,2);
                            float en_s = 0.5f * m * (vx_s*vx_s + vy_s*vy_s + vz_s*vz_s);
                            int newTrk = CreateElectron(e.trackID, sec.matrix(2,0), sec.matrix(0,0), sec.matrix(0,1), sec.matrix(0,2), vx_s, vy_s, vz_s, en_s, sec.processType);
                            sec.matrix(2,2) = static_cast<float>(newTrk);
                            double t_hit2 = physics->Point_de_contact2(sec.matrix, e.channel, cts2, alpha2, x2, x3, R, dia, pas);
                            ElectronState next_state;
                            Matrix3x3 state_for_push2 = sec.matrix;
                            if (t_hit2 == false) {
                                // Transport to pore exit (x3) before entering GAP2
                                Matrix3x3 M_exit2 = physics->RecuperationTo(sec.matrix, x3);
                                TrackElectronOutsidePore(sec.matrix, M_exit2, newTrk, cts2);
                                state_for_push2 = M_exit2;
                                next_state = GAP2;
                            } else if (t_hit2 == true) {
                                next_state = DEAD;
                            } else if (std::isnan(t_hit2) || t_hit2 <= 0) {
                                next_state = DEAD;
                            } else {
                                next_state = MCP2;
                            }
                            new_electrons.push_back({state_for_push2, next_state, newTrk, e.trackID, sec.processType, e.channel});
                        }
                        e.state = temp;
                    }
                    break;
                }
                case GAP2: {
                    const double dt_gap = 1.0; // 1 ps step
                    double dt = dt_gap;
                    double x_pred = e.state(0, 0) + e.state(1, 0) * dt + 0.5 * cts2 * dt * dt;
                    if (x_pred >= x4) {
                        double a = cts2 / 2.0;
                        double b = e.state(1, 0);
                        double c_val = e.state(0, 0) - x4;
                        dt = physics->Resolution(a, b, c_val);
                        if (std::isnan(dt) || dt <= 0) {
                            std::cout << "[Warn] GAP2: invalid dt to boundary -> DEAD" << std::endl;
                            e.status = DEAD;
                            break;
                        }
                    }
                    Matrix3x3 M_step = physics->Transporter2(e.state, dt, cts2, alpha2);
                    TrackElectronOutsidePore(e.state, M_step, e.trackID, cts2);
                    e.state = M_step;
                    if (e.state(0, 0) >= x4 - 1e-6) {
                        Matrix3x3 M_exact = physics->RecuperationTo(e.state, x4);
                        TrackElectronOutsidePore(e.state, M_exact, e.trackID, cts2);
                        e.state = M_exact;
                        e.status = ANODE;
                        FinalizeElectron(e.trackID, 1, M_exact(2,0), M_exact(0,0), M_exact(0,1), M_exact(0,2), M_exact(1,0), M_exact(1,1), M_exact(1,2), 0.5f * m * (M_exact(1,0)*M_exact(1,0) + M_exact(1,1)*M_exact(1,1) + M_exact(1,2)*M_exact(1,2)));
                    }
                    break;
                }
                default: break;
            }
        }
        for (auto& ne : new_electrons) electrons.push_back(ne);
    }
    // anode 도달 전자만 Resultat_final에 저장
    for (const auto& e : electrons) {
        if (e.status == ANODE) Resultat_final.push_back(e.state);
    }
    std::cout << "\nSimulation completed with " << Resultat_final.size() << " electrons reaching the anode" << std::endl;
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

// ------------------  Ensemble helpers under migration ------------------
void Simulation::ProcessMCP1Ensemble(double& cts1, double alpha1, double x0, double x1,
                                     double R, double dia, double pas, double m, double E0,
                                     int& nextTrackBaseID) {
    if (E1_emi_.empty() && E1_nonemi_.empty()) return; // nothing to process

    // 1. Determine next collision time (temps)
    double temps = 0.0;
    if (!E1_emi_.empty()) {
        temps = E1_emi_.front()(2, 1); // stored t_hit
    } else {
        // No pending collision; choose small step (1 ps) to propagate non-emi electrons
        temps = 1.0;
    }

    // 2. Move non-colliding electrons by temps
    std::vector<int> eraseIdx;
    for (size_t i = 0; i < E1_nonemi_.size(); ++i) {
        Matrix3x3 &M = E1_nonemi_[i];
        Matrix3x3 M_pass = physics->Transporter2(M, temps, cts1, alpha1);
        // track steps outside pore (optional)
        // If electron exits pore (x >= x1) -> push to GAP1 vector
        if (M_pass(0,0) >= x1) {
            // bring exactly to exit plane for continuity
            Matrix3x3 M_exit = physics->RecuperationTo(M_pass, x1);
            G1_.push_back(M_exit);
            eraseIdx.push_back(i);
        } else {
            E1_nonemi_[i] = M_pass;
        }
    }
    // erase in reverse order
    for (int idx = static_cast<int>(eraseIdx.size()) - 1; idx >= 0; --idx) {
        E1_nonemi_.erase(E1_nonemi_.begin() + eraseIdx[idx]);
    }

    // 3. Process the leading collision if available
    if (!E1_emi_.empty()) {
        Matrix3x3 lead = E1_emi_.front();
        E1_emi_.erase(E1_emi_.begin());

        // Move to collision point
        Matrix3x3 M_coll = physics->Transporter1(lead, temps, cts1, alpha1);

        // NOTE: Secondary emission not yet implemented in 3-A-1; simply discard for now
        // After collision, treat primary as absorbed (dead) so not added back.
    }

    // Debug print
    std::cout << "[Ensemble] after step : E1_emi=" << E1_emi_.size() << " E1_nonemi=" << E1_nonemi_.size() << " G1=" << G1_.size() << std::endl;
}

} // namespace MCPSim 