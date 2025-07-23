#include "MCPSimulation.h"
#include "MCPConfig.h"
#include <iostream>
#include <algorithm>
#include <random>
#include <cmath>
#include <exception>

#define DEBUG_MCP2

static int g_nextTrackID = 1;          // Global trackID

namespace MCPSim {

// ---------------------------------------------------------------------
// Output level helper (Track / Node / Step)
// ---------------------------------------------------------------------
OutputLevel GetOutputLevel(){
    auto& cfg = Config::getInstance();
    try{
        double v = cfg.get("outputLevel");
        int iv = static_cast<int>(std::round(v));
        switch(iv){
            case 0: return OutputLevel::kTrack;
            case 2: return OutputLevel::kStep;
            default: return OutputLevel::kNode;   // 1 or any other value
        }
    }catch(const std::exception&){
        // parameter not found → default
        return OutputLevel::kNode;
    }
}

// ───────────────────────────────────────────
// 0. Helpers
// ───────────────────────────────────────────
namespace {
inline double KE(const Matrix3x3& M,double m){
    double vx=M(1,0),vy=M(1,1),vz=M(1,2);
    return 0.5*m*(vx*vx+vy*vy+vz*vz);
}
} // anonymous

// ───────────────────────────────────────────
// 1. C-tors / Track helpers
// ───────────────────────────────────────────
Simulation::Simulation()
    : physics(std::make_unique<Physics>())
{
    Physics::initialize();
    tracks_.Reset();
    steps_.Reset();
    nSecOutMCP1_ = 0;
    // control gap-step recording based on output level (only full-step mode stores every 1 ps)
    recordGapSteps_ = (GetOutputLevel() == OutputLevel::kStep);
}

int  Simulation::CreateElectron(int parent, float t,
                                float x,float y,float z,
                                float vx,float vy,float vz,
                                float E, int proc){
    int ext=g_nextTrackID++;
    int in = tracks_.AddTrack(parent,t,x,y,z,vx,vy,vz,E,proc);
    trackIDMap_[ext]=in;
    // no amplification bookkeeping
    return ext;
}

void Simulation::AddElectronStep(int tid,float t,
                                 float x,float y,float z,
                                 float vx,float vy,float vz,
                                 float E){
    auto it=trackIDMap_.find(tid);
    if(it!=trackIDMap_.end())
        steps_.AddStep(it->second,t,x,y,z,vx,vy,vz,E);
}

void Simulation::FinalizeElectron(int tid,int st,float t,
                                  float x,float y,float z,
                                  float vx,float vy,float vz,
                                  float E){
    auto it=trackIDMap_.find(tid);
    if(it!=trackIDMap_.end())
        tracks_.FinalizeTrack(it->second,st,t,x,y,z,vx,vy,vz,E);
}

// ---------------------------------------------------------------------
// Record a single point (Node) depending on output level
// ---------------------------------------------------------------------
void Simulation::RecordNode(int tid, const Matrix3x3& M){
    AddElectronStep(tid,
                    float(M(2,0)),
                    float(M(0,0)), float(M(0,1)), float(M(0,2)),
                    float(M(1,0)), float(M(1,1)), float(M(1,2)),
                    float(KE(M, Config::getInstance().get("m"))));
}

// ---------------------------------------------------------------------
// Track electron segment outside pore with variable granularity
// ---------------------------------------------------------------------
void Simulation::TrackElectronOutsidePore(const Matrix3x3& A, const Matrix3x3& B,
                                          int tid, double cts)
{
    const auto level = GetOutputLevel();

    if(level == OutputLevel::kTrack){
        return;                                     // no positional recording
    }
    if(level == OutputLevel::kNode){
        RecordNode(tid, A);
        RecordNode(tid, B);
        return;                                     // endpoints only
    }

    // Full step (kStep) – 1 ps sampling along the segment
    const double dt = 1.0;                          // ps
    double T  = A(2,0);
    double T2 = B(2,0);
    double len = T2 - T;
    int n = std::max(1, static_cast<int>(len / dt));

    for(int i = 1; i <= n; ++i){
        double f = static_cast<double>(i) / n;
        double t = len * f;

        Matrix3x3 M = A;
        M.row(0) += (A.row(1) * t);
        M(0,0)  += 0.5 * cts * t * t;
        M(1,0)  += cts * t;
        M(2,0)   = T + t;

        RecordNode(tid, M);
    }
}

// ───────────────────────────────────────────
// 2. Run()
// ───────────────────────────────────────────
std::vector<Matrix3x3> Simulation::Run(double Einit){
    auto& C=Config::getInstance();
    double x0=C.get("x0"),x1=C.get("x1"),x2=C.get("x2"),
           x3=C.get("x3"),x4=C.get("x4");
    double alpha1=C.get("alpha1"),alpha2=C.get("alpha2");
    double R=C.get("R"),dia=C.get("dia"),pas=C.get("pas");
    double cts1=C.get("c_c"),cts2=C.get("c_c2");
    // Gap-specific coefficients (fixed, no dynamic scaling)
    double c_s = C.get("c_s");
    double c_s2 = C.get("c_s2");
    double m=C.get("m"),E0=C.get("E0"),q=C.get("q");
    double limite=C.get("limite"),Istrip=C.get("I_strip");

    anode_hits_.clear();

    // ── 초기 광전자 생성 & 채널 확인 ─────────────────
    Matrix3x3 P=physics->Pho_ele(Einit);
    P = physics->premiere_arrivee(P);
    auto hit=physics->Check_if_hit(P);
    // std::cout << "[DBG] hit=" << hit.first
    //           << ", y=" << P(0,1) << ", limite=" << limite << std::endl;
    if(!(hit.first && std::abs(P(0,1))<limite)){ return anode_hits_;}

    // ── Ensemble 벡터 선언 ──────────────────────────
    std::vector<Matrix3x3> E1_emi,E1_non,G1;
    std::vector<Matrix3x3> E2_emi,E2_non,G2;

    // ── 초기 3전자 MCP-1 벡터에 삽입 ────────────────
    for(int k=0;k<3;k++){
        Matrix3x3 M=P;
        M(0,0)-=0.1*k;
        int trk=CreateElectron(-1,M(2,0),M(0,0),M(0,1),M(0,2),
                               M(1,0),M(1,1),M(1,2),
                               float(KE(M,m)),mcp::PROCESS_SECONDARY);
        M(2,2)=trk;
        // --- Point-de-contact (MCP-1) : work in local plate-1 coordinates ---
        double th;
        {
            Matrix3x3 Mloc = M;                 // copy (keep global M for storage)
            const double pitch = dia + pas;
            int   nz   = int(std::round(Mloc(0,2) / pitch));

            Mloc(0,0) -= x0;                    // local x (plate entrance = 0 µm)
            Mloc(0,2) -= nz * pitch;            // local z (pore centre)

            th = physics->Point_de_contact2(Mloc, hit.second, cts1,
                                              alpha1,
                                              0.0, x1 - x0,   // x0/x1 in local frame
                                              R, dia, pas);
        }
        if (th == false) {
            // Electron exits MCP-1 without wall collision
            E1_non.push_back(M);
        } else if (th == true) {
            // No collision predicted inside the current pore segment → treat like E1_non
            E1_non.push_back(M);
        } else {
            // Positive time means collision scheduled; store in emission queue
            M(2, 1) = th;
            physics->ajouter_element_trie(E1_emi, M);
        }
    }

    // std::cout << "[INIT] after loop  "
    //           << "E1_emi=" << E1_emi.size()
    //           << ", E1_non=" << E1_non.size() << std::endl;

    // ── 상태 변수 ───────────────────────────────────
    int   Etat1=0,Etat2=0;
    double inst1=0,inst2=0;
    int   iter=0;

    // ────────────────────────────────────────────────
    while(true){
        iter++;
        // --- progress print every 100 iterations ---
        if(iter % 100 == 0){
            std::cout << "[Run] iteration " << iter
                      << ", E1_emi = " << E1_emi.size()
                      << ", E1_non = " << E1_non.size()
                      << ", G1 = " << G1.size()
                      << ", E2_emi = " << E2_emi.size()
                      << ", E2_non = " << E2_non.size()
                      << ", G2 = " << G2.size()
                      << ", anode hits = " << anode_hits_.size() << std::endl;
        }
        /*========== 1. MCP-1 Ensemble step =========*/
        if(!E1_emi.empty()){
            double dt=E1_emi.front()(2,1);
            // move non-emi
            std::vector<size_t> del;
            for(size_t i=0;i<E1_non.size();++i){
                auto M=physics->Transporter2(E1_non[i],dt,cts1,alpha1);
                int tid=int(E1_non[i](2,2));
                TrackElectronOutsidePore(E1_non[i],M,tid,cts1);
                if(M(0,0)>=x1){                          // exit -> GAP1 (use dynamic cts1 for final snap)
                    double dt_exit = physics->Resolution(cts1/2.0, M(1,0), M(0,0) - x1);
                    Matrix3x3 M_exit = physics->Transporter2(M, dt_exit, cts1, alpha1);
                    TrackElectronOutsidePore(M, M_exit, tid, cts1);
                    G1.push_back(M_exit);
                    ++nSecOutMCP1_;
                    del.push_back(i);
                }else E1_non[i]=M;
            }
            for(int i=del.size()-1;i>=0;--i) E1_non.erase(E1_non.begin()+del[i]);

            // move emi lead
            Matrix3x3 lead=E1_emi.front();E1_emi.erase(E1_emi.begin());
            Matrix3x3 Mc=physics->Transporter1(lead,dt,cts1,alpha1);
            // record full path inside MCP-1 (lead -> collision point)
            TrackElectronOutsidePore(lead, Mc, int(lead(2,2)), cts1);
            Mc(2,2)=lead(2,2);
            auto secs=physics->emi_sec(Mc,hit.second,alpha1,x0,R,dia,pas,m,E0);

            for(auto& s:secs){
                int tid=CreateElectron(int(Mc(2,2)),s.matrix(2,0),
                                       s.matrix(0,0),s.matrix(0,1),s.matrix(0,2),
                                       s.matrix(1,0),s.matrix(1,1),s.matrix(1,2),
                                       float(KE(s.matrix,m)),s.processType);
                s.matrix(2,2)=tid;
                double t2; // collision time placeholder
                // --- 2-D channel handling: recompute channel index & shift z to local ---
                {
                    // ① 현재 전자가 속한 채널 인덱스(n) 재계산
                    int chanIdx = physics->Check_if_hit(s.matrix).second;

                    // ② z-축 로컬 변환 후 충돌 시간 계산
                    double pitch = dia + pas;
                    int    nz    = int(std::round(s.matrix(0,2) / pitch));
                    double zc    = nz * pitch;

                    Matrix3x3 Mshift = s.matrix;
                    Mshift(0,2) -= zc;      // local z
                    Mshift(0,0) -= x0;      // local x

                    t2 = physics->Point_de_contact2(Mshift, chanIdx, cts1, alpha1,
                                                   0.0, x1 - x0, R, dia, pas);
                }
                if(t2==false){          // Electron leaves channel directly
                    E1_non.push_back(s.matrix);
                }else if(t2==true){     // No collision in this channel segment
                    continue;           // nothing to add
                }else{                  // Collision scheduled
                    s.matrix(2,1)=t2;
                    physics->ajouter_element_trie(E1_emi,s.matrix);
                }
                // std::cout << "[PUSH] t2=" << t2 << std::endl;
            }

            /*--- dynamic field MCP-1 ---*/
            double I1=0;
            auto addI=[&](const Matrix3x3& M){ if(M(0,0)<x1) I1+=q*M(1,0)/(x1-x0);};
            for(auto& M:E1_emi) addI(M);
            for(auto& M:E1_non) addI(M);
            double time=Mc(2,0);
            if(Etat1==0 && I1>=0.05*Istrip){
                cts1*=0.8; Etat1=1; inst1=time;
                auto res=physics->Rearrangement(E1_emi,E1_non,hit.second,cts1,alpha1,x0,x1,R,dia,pas);
                E1_emi=res.first;E1_non=res.second;
            }else if(Etat1==1 && time>=inst1+5.0 && I1<0.05*Istrip){
                cts1/=0.8; Etat1=0; inst1=time;
                auto res=physics->Rearrangement(E1_emi,E1_non,hit.second,cts1,alpha1,x0,x1,R,dia,pas);
                E1_emi=res.first;E1_non=res.second;
            }
        }

        /*--- MCP-1 propagate when no collisions are scheduled ---*/
        if(E1_emi.empty() && !E1_non.empty()){
            const double dtsmall = 1.0;
            std::vector<size_t> delIdx;
            for(size_t i=0;i<E1_non.size();++i){
                auto M = physics->Transporter2(E1_non[i], dtsmall, cts1, alpha1);
                int tid = int(E1_non[i](2,2));
                TrackElectronOutsidePore(E1_non[i], M, tid, cts1);

                if(M(0,0) >= x1){               // exited MCP-1 → GAP-1
                    double dt_exit = physics->Resolution(cts1/2.0, M(1,0), M(0,0) - x1);
                    Matrix3x3 M_exit = physics->Transporter2(M, dt_exit, cts1, alpha1);
                    TrackElectronOutsidePore(E1_non[i], M_exit, tid, cts1);
                    G1.push_back(M_exit);
                    delIdx.push_back(i);
                }else{
                    // r>R check
                    auto hh = physics->Check_if_hit(M);
                    if(!hh.first){               // left pore cylinder – absorb
                        FinalizeElectron(tid, 0, M(2,0), M(0,0),M(0,1),M(0,2),
                                         M(1,0),M(1,1),M(1,2), float(KE(M,m)));
                        delIdx.push_back(i);
                    }else{
                        E1_non[i] = M;          // still inside pore
                    }
                }
            }
            for(int k=delIdx.size()-1; k>=0; --k) E1_non.erase(E1_non.begin()+delIdx[k]);
        }

        /*========== 2. GAP-1 propagate =========*/
        {
            std::vector<size_t> del;
            for(size_t i=0;i<G1.size();++i){
                double dt=1.0;
                double xpred=G1[i](0,0)+G1[i](1,0)*dt+0.5*c_s*dt*dt;
                if(xpred>=x2){
                    dt=physics->Resolution(c_s/2.0,G1[i](1,0),G1[i](0,0)-x2);
                }
                auto M=physics->Transporter2(G1[i],dt,c_s,alpha1);
                TrackElectronOutsidePore(G1[i],M,int(G1[i](2,2)),c_s);
                
                // MCP2 진입 조건: x2 근처 범위에서 검사 (정확한 x2가 아닌 범위 기반)
                const double entrance_tolerance = 0.1;  
                if(M(0,0) >= x2 - entrance_tolerance){
                    M=physics->RecuperationTo(M,x2);
                    TrackElectronOutsidePore(G1[i],M,int(M(2,2)),c_s);
                    auto h=physics->Check_if_hit(M);
                    if(h.first && std::abs(M(0,1))<limite){
                        // 음수 채널 인덱스 허용 (정상 처리)
                        
                        // ★ 디버그: 채널 인덱스와 실제 거리 확인 ★
                        double pitch = dia + pas;
                        // 올바른 거리 계산: Check_if_hit과 동일한 방식 사용
                        double y_prime = M(0,1) - tan(alpha2) * (M(0,0) - x2);
                        double center_y_local = (h.second + 0.5) * pitch;
                        double dy = y_prime - center_y_local;
                        
                        // z축도 포어 중심 기준으로 계산
                        int nz = int(round(M(0,2) / pitch));
                        double center_z = nz * pitch;
                        double dz = M(0,2) - center_z;
                        
                        double r_actual = sqrt(dy*dy + dz*dz);
                        double t2; // collision time placeholder
                        // --- 2-D channel handling: shift z so that local pore axis is at z_c ---
                        {
                            double pitch = dia + pas;
                            int    nz    = int(std::round(M(0,2) / pitch));
                            double zc    = nz * pitch;

                            Matrix3x3 Mshift = M;
                            // Shift global x to local coordinate for MCP-2 (x2 becomes 0)
                            Mshift(0,0) -= x2;
                            Mshift(0,2) -= zc;   // translate into local channel frame (z_c → 0)

                            t2 = physics->Point_de_contact2(Mshift, h.second, cts2, alpha2,
                                                            0.0, x3 - x2, R, dia, pas);
                        }
                        
                        // Record entry Point_de_contact2 result for this electron (internal track index)
                        {
                            int extId = int(M(2,2));
                            auto pitET = trackIDMap_.find(extId);
                            if(pitET != trackIDMap_.end()){
                                double storeVal;
                                if(t2==false) storeVal = -2.0;      // flag: 'false'
                                else if(t2==true) storeVal = -1.0;  // flag: 'true'
                                else storeVal = t2;                 // positive collision time
                                // entryT_MCP2_[pitET->second] = storeVal; // Removed
                            }
                        }

                        if(t2==false){
                            // keep z coordinate; do not overwrite with channel index
                            E2_non.push_back(M);
                        }else if(t2==true){     // No collision yet; keep physical z
                            // Keep z coordinate intact, but store n in a temp list if needed
                            E2_non.push_back(M);
                        }else{                  // Collision scheduled -> add to emi list
                            // keep z coordinate; channel index stored only in variable `h.second`
                            M(2,1)=t2;
                            physics->ajouter_element_trie(E2_emi,M);
                        }
                        del.push_back(i);
                    } else {
                        // 포어에 못 들어간 전자는 끝
                        int tid = int(M(2,2));
                        FinalizeElectron(tid, 0, M(2,0),   // status=0: lost
                                         M(0,0),M(0,1),M(0,2),
                                         M(1,0),M(1,1),M(1,2),
                                         float(KE(M,m)));
                        // Track 끝 위치 한 번만 기록하고 버림
                        TrackElectronOutsidePore(G1[i], M, tid, c_s);
                        del.push_back(i);                  // G1 리스트에서 제거
                    }
                }else G1[i]=M;
            }
            for(int i=del.size()-1;i>=0;--i) G1.erase(G1.begin()+del[i]);
        }

        /*========== 3. MCP-2 Ensemble step =========*/
        if(!E2_emi.empty()){
            double dt=E2_emi.front()(2,1);
            
            std::vector<size_t> delN;
            for(size_t i=0;i<E2_non.size();++i){
                auto M=physics->Transporter2(E2_non[i],dt,cts2,alpha2);
                int tid=int(E2_non[i](2,2));
                TrackElectronOutsidePore(E2_non[i],M,tid,cts2);
                if(M(0,0)>=x3){ 
                    double dt_exit = physics->Resolution(cts2/2.0, M(1,0), M(0,0) - x3);
                    Matrix3x3 M_exit = physics->Transporter2(M, dt_exit, cts2, alpha2);
                    TrackElectronOutsidePore(E2_non[i], M_exit, tid, cts2);
                    G2.push_back(M_exit);
                    delN.push_back(i);} 
                else {
                    // Additional geometry check: has the electron left the pore cylinder?
                    auto hh = physics->Check_if_hit(M);
                    if(!hh.first){
                        // Absorb on silica wall – stop tracking
                        FinalizeElectron(tid, 0, M(2,0),
                                         M(0,0), M(0,1), M(0,2),
                                         M(1,0), M(1,1), M(1,2),
                                         float(KE(M, m)));
                        delN.push_back(i);
                    }else{
                        E2_non[i]=M;
                    }
                }
            }
            for(int i=delN.size()-1;i>=0;--i) E2_non.erase(E2_non.begin()+delN[i]);

            Matrix3x3 lead=E2_emi.front();E2_emi.erase(E2_emi.begin());
            Matrix3x3 Mc=physics->Transporter1(lead,dt,cts2,alpha2);
            // record full path inside MCP-2 (lead -> collision point)
            TrackElectronOutsidePore(lead, Mc, int(lead(2,2)), cts2);
            int ch = physics->Check_if_hit(lead).second;
            auto secs = physics->emi_sec(Mc, ch, alpha2, x2, R, dia, pas, m, E0);

            // no amplification bookkeeping

            for(auto& s:secs){
                int tid=CreateElectron(int(Mc(2,2)),s.matrix(2,0),
                                       s.matrix(0,0),s.matrix(0,1),s.matrix(0,2),
                                       s.matrix(1,0),s.matrix(1,1),s.matrix(1,2),
                                       float(KE(s.matrix,m)),s.processType);
                s.matrix(2,2)=tid;
                double t2; // collision time placeholder
                // --- 2-D channel handling: recompute channel index & shift z to local ---
                {
                    // ① 채널 인덱스 재계산 (글로벌 좌표 기준)
                    int chanIdx = physics->Check_if_hit(s.matrix).second;
                    
                    // 음수 채널 인덱스 허용 (secondary electrons도 정상 처리)

                    // ② z-축 로컬 변환 후 충돌 시간 계산
                    double pitch = dia + pas;
                    int    nz    = int(std::round(s.matrix(0,2) / pitch));
                    double zc    = nz * pitch;

                    Matrix3x3 Mshift = s.matrix;
                    // Shift global x to local coordinate for MCP-2 (x2 becomes 0)
                    Mshift(0,0) -= x2;
                    Mshift(0,2) -= zc;

                    t2 = physics->Point_de_contact2(Mshift, chanIdx, cts2, alpha2,
                                                   0.0, x3 - x2, R, dia, pas);
                }
                if(t2==false){           // already outside MCP-2
                    // do not overwrite z coordinate
                    E2_non.push_back(s.matrix);
                }else if(t2==true){      // no collision yet; keep z
                    // keep z coordinate; don't overwrite
                    E2_non.push_back(s.matrix);
                }else{                   // collision scheduled
                    // keep z coordinate untouched
                    s.matrix(2,1)=t2;
                    physics->ajouter_element_trie(E2_emi,s.matrix);
                }
            }

            /*--- dynamic field MCP-2 ---*/
            double I2=0; auto addI2=[&](const Matrix3x3& M){ if(M(0,0)<x3) I2+=q*M(1,0)/(x3-x2);};
            for(auto& M:E2_emi) addI2(M);for(auto& M:E2_non) addI2(M);
            double time=Mc(2,0);
            if(Etat2==0 && I2>=0.05*Istrip){
                cts2*=0.8; Etat2=1; inst2=time;
                auto res=physics->Rearrangement(E2_emi,E2_non,ch,cts2,alpha2,x2,x3,R,dia,pas);
                E2_emi=res.first;E2_non=res.second;
            }else if(Etat2==1 && time>=inst2+5.0 && I2<0.05*Istrip){
                cts2/=0.8; Etat2=0; inst2=time;
                auto res=physics->Rearrangement(E2_emi,E2_non,ch,cts2,alpha2,x2,x3,R,dia,pas);
                E2_emi=res.first;E2_non=res.second;
            }
        }

        /*--- MCP-2 propagate when no collisions are scheduled ---*/
        if(E2_emi.empty() && !E2_non.empty()){
            const double dtsmall = 1.0; 
            std::vector<size_t> delIdx;
            for(size_t i=0;i<E2_non.size();++i){
                auto M = physics->Transporter2(E2_non[i], dtsmall, cts2, alpha2);
                int tid = int(E2_non[i](2,2));
                TrackElectronOutsidePore(E2_non[i], M, tid, cts2);

                if(M(0,0) >= x3){               // exited MCP-2 → GAP-2
                    double dt_exit = physics->Resolution(cts2/2.0, M(1,0), M(0,0) - x3);
                    Matrix3x3 M_exit = physics->Transporter2(M, dt_exit, cts2, alpha2);
                    TrackElectronOutsidePore(E2_non[i], M_exit, tid, cts2);
                    G2.push_back(M_exit);
                    delIdx.push_back(i);
                }else{
                    // r>R check for small-step propagation as well
                    auto hh = physics->Check_if_hit(M);
                    if(!hh.first){
                        FinalizeElectron(tid, 0, M(2,0), M(0,0),M(0,1),M(0,2),
                                         M(1,0),M(1,1),M(1,2), float(KE(M,m)));
                        delIdx.push_back(i);
                    }else{
                        E2_non[i] = M;          // still inside pore
                    }
                }
            }
            for(int k = delIdx.size()-1; k>=0; --k) E2_non.erase(E2_non.begin()+delIdx[k]);
        }

        /*========== 4. GAP-2 propagate =========*/
        {
            const double dt=0.05;
            std::vector<size_t> del;
            for(size_t i=0;i<G2.size();++i){
                auto M=physics->Transporter2(G2[i],dt,c_s2,alpha2);
                int tid=int(G2[i](2,2));
                TrackElectronOutsidePore(G2[i],M,tid,c_s2);
                if(M(0,0)>=x4){
                    Matrix3x3 Me=physics->RecuperationTo(M,x4);
                    TrackElectronOutsidePore(M,Me,tid,c_s2);
                    anode_hits_.push_back(Me);
                    FinalizeElectron(tid,1,Me(2,0),Me(0,0),Me(0,1),Me(0,2),
                                     Me(1,0),Me(1,1),Me(1,2),float(KE(Me,m)));
                    del.push_back(i);
                }else G2[i]=M;
            }
            for(int i=del.size()-1;i>=0;--i) G2.erase(G2.begin()+del[i]);
        }

        /*========== 5. 종료 조건 =========*/
        // if(anode_hits_.size()>=200) break;
        if(E1_emi.empty()&&E1_non.empty()&&G1.empty()&&
           E2_emi.empty()&&E2_non.empty()&&G2.empty()) break;
    }

    std::cout << "Simulation completed. anode=" << anode_hits_.size()
              << ", MCP1→Gap1 secondaries = " << nSecOutMCP1_ << std::endl;
    return anode_hits_;
}

// ───────────────────────────────────────────
// 3. Save / ConvertEvent (그대로 유지)
// ───────────────────────────────────────────
void Simulation::Save(const std::vector<Matrix3x3>& res,double E,const std::string& f){
    MCPRootManager r;
    if(!r.OpenFile(f)) throw std::runtime_error("root file open error");
    r.WriteEvt(ConvertEvent(res,E)); r.Close();
}
mcp::Event Simulation::ConvertEvent(const std::vector<Matrix3x3>& res,double Ein){
    mcp::Event evt; 
    evt.eventInfo.initialEnergy = Ein;
    evt.config.LoadFromConfig();

    // -----------------------------------------------------------------
    // Filter: if outputLevel == kTrack, keep 애노드에 도달한 Track 만
    // -----------------------------------------------------------------
    if(GetOutputLevel() == OutputLevel::kTrack){
        mcp::Track filt; filt.Reset();
        const auto& src = tracks_;
        for(int i=0;i<src.nTracks;++i){
            if(src.isAnode[i]==1){
                filt.nTracks++;
                filt.trackID.push_back(      src.trackID[i]);
                filt.parentID.push_back(     src.parentID[i]);
                filt.birthTime.push_back(    src.birthTime[i]);
                filt.birthPosX.push_back(    src.birthPosX[i]);
                filt.birthPosY.push_back(    src.birthPosY[i]);
                filt.birthPosZ.push_back(    src.birthPosZ[i]);
                filt.birthVelX.push_back(    src.birthVelX[i]);
                filt.birthVelY.push_back(    src.birthVelY[i]);
                filt.birthVelZ.push_back(    src.birthVelZ[i]);
                filt.birthEnergy.push_back(  src.birthEnergy[i]);
                filt.processType.push_back(  src.processType[i]);
                filt.isAnode.push_back(      src.isAnode[i]);  // always 1
                filt.finalTime.push_back(    src.finalTime[i]);
                filt.finalPosX.push_back(    src.finalPosX[i]);
                filt.finalPosY.push_back(    src.finalPosY[i]);
                filt.finalPosZ.push_back(    src.finalPosZ[i]);
                filt.finalVelX.push_back(    src.finalVelX[i]);
                filt.finalVelY.push_back(    src.finalVelY[i]);
                filt.finalVelZ.push_back(    src.finalVelZ[i]);
                filt.finalEnergy.push_back(  src.finalEnergy[i]);
            }
        }
        evt.tracks = std::move(filt);
    }else{
        evt.tracks = tracks_;
    }

    // steps_: 이미 outputLevel==kTrack 일 때는 기록이 없으므로 그대로 복사해도 비어 있음
    evt.steps = steps_;
    return evt;
}

} // namespace MCPSim 