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
        
        // Calculate the intermediate position (using physical equations)
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
        double time_1_hit = physics->Point_de_contact2(M, check.second, cts1, alpha1, x0, x1, R, dia, pas);
            M(2, 1) = time_1_hit;
        electrons.push_back({M, MCP1, trk, -1, PROCESS_SECONDARY, check.second});
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

        // -------------------------------------------------
        // Current-based electric field adjustment (B)
        // -------------------------------------------------
        double I1 = 0.0, I2 = 0.0;
        for (const auto& ec : electrons) {
            if (ec.status == MCP1) {
                if (ec.state(0,0) < x1) {
                    I1 += q * ec.state(1,0) / (x1 - x0);
                }
            } else if (ec.status == MCP2) {
                if (ec.state(0,0) < x3) {
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
                        // Electron exits MCP1; keep current state and immediately enter GAP1 (v3_track behaviour)
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
                            // Use the same channel as parent for secondary in MCP1
                            double t_hit2 = physics->Point_de_contact2(sec.matrix, e.channel, cts1, alpha1, x0, x1, R, dia, pas);
                            ElectronState next_state;
                            if (t_hit2 == false) { next_state = GAP1; /* std::cout << "[MCP1 secondary→GAP1]" << std::endl; */ }
                            else if (t_hit2 == true) { next_state = DEAD; /* std::cout << "[MCP1 secondary→DEAD] (true)" << std::endl; */ }
                            else if (std::isnan(t_hit2) || t_hit2 <= 0) { next_state = DEAD; /* std::cout << "[Warn] MCP1 secondary: t_hit2 비정상(" << t_hit2 << ") → DEAD 처리" << std::endl; */ }
                            else { next_state = MCP1; /* std::cout << "[MCP1 secondary→MCP1]" << std::endl; */ }
                            new_electrons.push_back({sec.matrix, next_state, newTrk, e.trackID, sec.processType, e.channel});
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
                    // If boundary reached, align exactly and enter MCP2
                    if (e.state(0, 0) >= x2 - 1e-6) {
                        Matrix3x3 M_exact = physics->RecuperationTo(e.state, x2);
                        TrackElectronOutsidePore(e.state, M_exact, e.trackID, cts1);
                        e.state = M_exact;
                        // Determine channel index upon entry to MCP2
                        int n_mcp2_ent;
                        if ((e.state(0, 1) - 1) < 0) {
                            n_mcp2_ent = int((e.state(0, 1) - 1) / (dia + pas)) - 1;
                        } else {
                            n_mcp2_ent = int((e.state(0, 1) - 1) / (dia + pas));
                        }
                        e.channel = n_mcp2_ent;
                        e.status = MCP2;
                    }
                    break;
                }
                case MCP2: {
                    int n_mcp2 = e.channel; // stored channel index
                    double t_hit = physics->Point_de_contact2(e.state, n_mcp2, cts2, alpha2, x2, x3, R, dia, pas);
                    if (t_hit == false) {
                        // Electron exits MCP2; keep current state and immediately enter GAP2 (v3_track behaviour)
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
                            // Same channel as parent
                            double t_hit2 = physics->Point_de_contact2(sec.matrix, e.channel, cts2, alpha2, x2, x3, R, dia, pas);
                            ElectronState next_state;
                            if (t_hit2 == false) { next_state = GAP2; /* std::cout << "[MCP2 secondary→GAP2]" << std::endl; */ }
                            else if (t_hit2 == true) { next_state = DEAD; /* std::cout << "[MCP2 secondary→DEAD] (true)" << std::endl; */ }
                            else if (std::isnan(t_hit2) || t_hit2 <= 0) { next_state = DEAD; /* std::cout << "[Warn] MCP2 secondary: t_hit2 비정상(" << t_hit2 << ") → DEAD 처리" << std::endl; */ }
                            else { next_state = MCP2; /* std::cout << "[MCP2 secondary→MCP2]" << std::endl; */ }
                            new_electrons.push_back({sec.matrix, next_state, newTrk, e.trackID, sec.processType, e.channel});
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

} // namespace MCPSim 