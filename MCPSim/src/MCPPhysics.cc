#include "MCPPhysics.h"
#include "MCPConfig.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>

namespace MCPSim {

// Static member variable definitions
std::random_device Physics::rd;
std::mt19937 Physics::gen(rd());
double Physics::x0;
double Physics::x1;
double Physics::x2;
double Physics::q;
double Physics::R;
double Physics::limite;
double Physics::P_inf;
double Physics::P_1e;
double Physics::P_1r;
double Physics::Ee;
double Physics::Er;
double Physics::W;
double Physics::p;
double Physics::r;
double Physics::e1;
double Physics::e2;
double Physics::r1;
double Physics::r2;
double Physics::gama_ts0;
double Physics::t1;
double Physics::t2;
double Physics::E_ts;
double Physics::t3;
double Physics::t4;
double Physics::s;
double Physics::m;
double Physics::dia;
double Physics::pas;
double Physics::alpha;
double Physics::E0;
double Physics::diff_pot;
double Physics::c_e;
double Physics::c_s;
double Physics::c_c;
double Physics::Resistance;

// Static initialization function
void Physics::initialize() {
    x0 = Config::getInstance().get("x0");
    x1 = Config::getInstance().get("x1");
    x2 = Config::getInstance().get("x2");
    q = Config::getInstance().get("q");
    R = Config::getInstance().get("R");
    limite = Config::getInstance().get("limite");
    P_inf = Config::getInstance().get("P_inf");
    P_1e = Config::getInstance().get("P_1e");
    P_1r = Config::getInstance().get("P_1r");
    Ee = Config::getInstance().get("Ee");
    Er = Config::getInstance().get("Er");
    W = Config::getInstance().get("W");
    p = Config::getInstance().get("p");
    r = Config::getInstance().get("r");
    e1 = Config::getInstance().get("e1");
    e2 = Config::getInstance().get("e2");
    r1 = Config::getInstance().get("r1");
    r2 = Config::getInstance().get("r2");
    gama_ts0 = Config::getInstance().get("gama_ts0");
    t1 = Config::getInstance().get("t1");
    t2 = Config::getInstance().get("t2");
    E_ts = Config::getInstance().get("E_ts");
    t3 = Config::getInstance().get("t3");
    t4 = Config::getInstance().get("t4");
    s = Config::getInstance().get("s");
    m = Config::getInstance().get("m");
    dia = Config::getInstance().get("dia");
    pas = Config::getInstance().get("pas");
    alpha = Config::getInstance().get("alpha");
    E0 = Config::getInstance().get("E0");
    diff_pot = Config::getInstance().get("diff_pot");
    c_e = Config::getInstance().get("c_e");
    c_s = Config::getInstance().get("c_s");
    c_c = Config::getInstance().get("c_c");
    Resistance = Config::getInstance().get("Resistance");
}

double Physics::comb(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    double c = 1.0;
    for (int i = 0; i < k; i++) {
        c = c * (n - i) / (i + 1);
    }
    return c;
}

// Helper function for solve_quartic
std::vector<double> Physics::solve_quartic(double a4, double a3, double a2, double a1, double a0) {
    std::vector<double> roots;
    
    // a4x^4 + a3x^3 + a2x^2 + a1x + a0 = 0
    // 표준형태로 정규화: x^4 + px^3 + qx^2 + rx + s = 0
    double p = a3 / a4;
    double q = a2 / a4;
    double r = a1 / a4;
    double s = a0 / a4;
    
    // Companion matrix
    Eigen::Matrix4d companion;
    companion << 0, 0, 0, -s,
                 1, 0, 0, -r,
                 0, 1, 0, -q,
                 0, 0, 1, -p;
    
    // Calculate eigenvalues (eigenvalues = roots of the polynomial)
    Eigen::EigenSolver<Eigen::Matrix4d> solver(companion);
    
    // Extract complex roots
    Eigen::Vector4cd eigenvalues = solver.eigenvalues();
    
    // Filter real roots
    for (int i = 0; i < 4; i++) {
        if (std::abs(eigenvalues(i).imag()) < 1e-10) {
            roots.push_back(eigenvalues(i).real());
        }
    }
    
    return roots;
}

double Physics::gama_e0(double E) {
    return P_inf + (P_1e - P_inf) * exp(-pow(fabs(E - Ee) / W, p) / p);
}

double Physics::gama_r0(double E) {
    return P_1r * (1 - exp(-pow(E / Er, r)));
}

double Physics::gama_e(double E, double teta) {
    return gama_e0(E) * (1 + e1 * (1 - pow(cos(teta), e2)));
}

double Physics::gama_r(double E, double teta) {
    return gama_r0(E) * (1 + r1 * (1 - pow(cos(teta), r2)));
}

double Physics::D(double x, double s) {
    return (s * x) / (s - 1 + pow(x, s));
}

double Physics::gama_teta(double teta) {
    return gama_ts0 * (1 + t1 * (1 - pow(cos(teta), t2)));
}

double Physics::E_teta(double teta) {
    return E_ts * (1 + t3 * (1 - pow(cos(teta), t4)));
}

double Physics::gama_ts(double E, double teta) {
    return gama_teta(teta) * D((E / E_teta(teta)), s);
}

double Physics::P_prime(double para, int M, int n) {
    return comb(M, n) * pow(para, n) * pow(1 - para, M - n);
}

double Physics::P(double gamae, double gamar, double para, int M, int n) {
    if (n == 0) {
        return (1 - gamae - gamar) * P_prime(para, M, n);
    } else if (n == 1) {
        return (1 - gamae - gamar) * P_prime(para, M, n) + gamae + gamar;
    } else {
        return (1 - gamae - gamar) * P_prime(para, M, n);
    }
}

int Physics::Tirage(const std::vector<int>& valeurs, const std::vector<double>& probabilites) {
    std::discrete_distribution<int> dist(probabilites.begin(), probabilites.end());
    return valeurs[dist(gen)];
}

int Physics::Generate(double E, double teta) {
    int M = 20;
    double gamae = gama_e(E, teta);
    double gamar = gama_r(E, teta);
    double gamats = gama_ts(E, teta);
    double gama_prime = gamats / (1 - gamae - gamar);
    double para = gama_prime / M;
    
    std::vector<double> probabilities(M + 1);
    std::vector<int> values(M + 1);
    
    for (int i = 0; i <= M; i++) {
        probabilities[i] = P(gamae, gamar, para, M, i);
        values[i] = i;
    }
    
    return Tirage(values, probabilities);
}

double Physics::Delta(double a, double b, double c) {
    return b * b - 4 * a * c;
}

double Physics::Resolution(double a, double b, double c) {
    if (Delta(a, b, c) < 0) {
        return false;
    } else {
        double racine_delta = sqrt(Delta(a, b, c));
        return (-b + racine_delta) / (2 * a);
    }
}

Vector3d Physics::Rot(const Vector3d& v, double teta) {
    Eigen::Matrix3d rotation;
    rotation << cos(teta), -sin(teta), 0,
                sin(teta), cos(teta), 0,
                0, 0, 1;
    return rotation * v;
}

Matrix3x3 Physics::Pho_ele(double Energie_photon) {
    double norme_vitesse = sqrt(2 * Energie_photon / m);
    
    Matrix3x3 Mpe = Matrix3x3::Zero();
    Mpe(0, 1) = 6;
    Mpe(1, 0) = norme_vitesse;
    Mpe(1, 1) = 0;
    
    return Mpe;
}

Matrix3x3 Physics::premiere_arrivee(Matrix3x3 Matrice_photo_electron) {
    double a = c_e / 2;
    double b = Matrice_photo_electron(1, 0);
    double c = -x0;
        
    double t = Resolution(a, b, c);
    double y = Matrice_photo_electron(1, 1) * t + Matrice_photo_electron(0, 1);
    double vx = Matrice_photo_electron(1, 0) + c_e * t;
    
    Matrice_photo_electron(0, 0) = x0;
    Matrice_photo_electron(0, 1) = y;
    Matrice_photo_electron(1, 0) = vx;
    Matrice_photo_electron(2, 0) = t;
    Matrice_photo_electron(2, 1) = 0;
    
    return Matrice_photo_electron;
}

std::pair<bool, int> Physics::Check_if_hit(const Matrix3x3& Matrice_arrive) {
    int n, m;
    
    if ((Matrice_arrive(0, 1) - 1) < 0) {
        n = int((Matrice_arrive(0, 1) - 1) / (dia + pas)) - 1;
    } else {
        n = int((Matrice_arrive(0, 1) - 1) / (dia + pas));
    }
    
    if ((Matrice_arrive(0, 1) + 1) < 0) {
        m = int((Matrice_arrive(0, 1) + 1) / (dia + pas)) - 1;
    } else {
        m = int((Matrice_arrive(0, 1) + 1) / (dia + pas));
    }
    
    if (n == (Matrice_arrive(0, 1) - 1) / (dia + pas)) {
        return {false, n};
    } else if (n != m) {
        return {false, n};
    } else {
        return {true, n};
    }
}

double Physics::Point_de_contact2(const Matrix3x3& Mat, int n, double cts) {
    double Mat00 = Mat(0, 0), Mat01 = Mat(0, 1), Mat02 = Mat(0, 2);
    double Mat10 = Mat(1, 0), Mat11 = Mat(1, 1), Mat12 = Mat(1, 2);
    
    // 4th order polynomial coefficients
    double a4 = ((cts * cts) * (tan(alpha) * tan(alpha)) / 4);
    double a3 = (Mat10 * cts * tan(alpha) * tan(alpha) - Mat11 * cts * tan(alpha));
    double a2 = (Mat00 * cts * tan(alpha) * tan(alpha) - Mat01 * cts * tan(alpha) + 
                (Mat10 * Mat10) * tan(alpha) * tan(alpha) - 2 * Mat10 * Mat11 * tan(alpha) + 
                Mat11 * Mat11 + Mat12 * Mat12 + cts * dia * n * tan(alpha) + 
                cts * dia * tan(alpha) / 2 + cts * n * pas * tan(alpha) + 
                cts * pas * tan(alpha) / 2 - cts * x0 * tan(alpha) * tan(alpha));
    double a1 = (2 * Mat00 * Mat10 * tan(alpha) * tan(alpha) - 2 * Mat00 * Mat11 * tan(alpha) - 
                2 * Mat01 * Mat10 * tan(alpha) + 2 * Mat01 * Mat11 + 2 * Mat02 * Mat12 + 
                2 * Mat10 * dia * n * tan(alpha) + Mat10 * dia * tan(alpha) + 
                2 * Mat10 * n * pas * tan(alpha) + Mat10 * pas * tan(alpha) - 
                2 * Mat10 * x0 * tan(alpha) * tan(alpha) - 2 * Mat11 * dia * n - 
                Mat11 * dia - 2 * Mat11 * n * pas - Mat11 * pas + 
                2 * Mat11 * x0 * tan(alpha));
    double a0 = (Mat00 * Mat00) * tan(alpha) * tan(alpha) - 2 * Mat00 * Mat01 * tan(alpha) + 
                2 * Mat00 * dia * n * tan(alpha) + Mat00 * dia * tan(alpha) + 
                2 * Mat00 * n * pas * tan(alpha) + Mat00 * pas * tan(alpha) - 
                2 * Mat00 * x0 * tan(alpha) * tan(alpha) + Mat01 * Mat01 - 
                2 * Mat01 * dia * n - Mat01 * dia - 2 * Mat01 * n * pas - 
                Mat01 * pas + 2 * Mat01 * x0 * tan(alpha) + Mat02 * Mat02 - 
                R * R + (dia * dia) * n * n + (dia * dia) * n + (dia * dia) / 4 + 
                2 * dia * (n * n) * pas + 2 * dia * n * pas - 
                2 * dia * n * x0 * tan(alpha) + dia * pas / 2 - 
                dia * x0 * tan(alpha) + (n * n) * pas * pas + n * pas * pas - 
                2 * n * pas * x0 * tan(alpha) + (pas * pas) / 4 - 
                pas * x0 * tan(alpha) + (x0 * x0) * tan(alpha) * tan(alpha);
    
    // Solve 4th order polynomial numerically
    // In actual implementation, it is recommended to use Eigen, GSL, etc.
    std::vector<double> roots = solve_quartic(a4, a3, a2, a1, a0);
    std::vector<double> valid_solutions;
    
    for (double root : roots) {
        double yf = Mat01 + root * Mat11;
        double xf = (cts / 2) * root * root + Mat00 + root * Mat10;
        double f = yf - tan(alpha) * xf + tan(alpha) * x0 - (pas + dia) * n - (pas + dia) / 2;
        double zf = Mat02 + root * Mat12;
        double final = sqrt(f * f + zf * zf);
        
        if (fabs(final - R) < 0.001 && root > 0.001) {
            valid_solutions.push_back(root);
        }
    }
    
    if (valid_solutions.empty()) {
        return true;
    } else {
        double t = *std::min_element(valid_solutions.begin(), valid_solutions.end());
        double x = Mat(0, 0) + Mat(1, 0) * t + cts * (t * t) / 2;
        
        if (x >= x1) {
            return false;
        } else if (x <= x0) {
            return true;
        } else {
            return t;
        }
    }
}

std::vector<ElectronProcess> Physics::emi_sec(const Matrix3x3& Mat, int n) {
    std::vector<ElectronProcess> Resultat;
    
    // Calculate energy
    double E = 0.5 * m * (pow(Mat(1, 0), 2) + pow(Mat(1, 1), 2) + pow(Mat(1, 2), 2));
    
    // Calculate teta
    double y_r = Mat(0, 1) - tan(alpha) * Mat(0, 0) + tan(alpha) * x0 - (pas + dia) * n - (pas + dia) / 2;
    double z_r = Mat(0, 2);
    
    double teta = atan2(z_r, y_r);
    if (teta < 0) {
        teta += 2 * M_PI;
    }

    // Create matrix
    Matrix3x3 ModifiedMat = Mat;
    ModifiedMat(0, 1) = R * cos(teta) + tan(alpha) * (Mat(0, 0) - x0) + (pas + dia) * n + (pas + dia) / 2;
    ModifiedMat(0, 2) = R * sin(teta);
    
    // Calculate
    Vector3d e_r(0, cos(teta), sin(teta));
    Vector3d e_o(0, -sin(teta), cos(teta));
    Vector3d z(0, 0, 1);
    Vector3d x(1, 0, 0);
    
    // Apply rotation matrix
    Vector3d normal = Rot(-e_r, alpha);
    Vector3d x_rot = Rot(x, alpha);
    Vector3d mormal = normal.cross(x_rot);
    
    // Calculate velocity vector and angle
    double norme = sqrt(pow(Mat(1, 0), 2) + pow(Mat(1, 1), 2) + pow(Mat(1, 2), 2));
    Vector3d vitesse(Mat(1, 0) / norme, Mat(1, 1) / norme, Mat(1, 2) / norme);
    double scalaire = normal.dot(-vitesse);
    double angle = acos(scalaire);
    
    if (angle > 89 * M_PI / 180) {
        angle = 89 * M_PI / 180;
    }
    
    // Extract parent trackID
    int parentTrackID = static_cast<int>(Mat(2, 2));
    
    // Generate number of electrons
    int nombre_elec = Generate(E, angle);
    
    if (nombre_elec != 0) {
        if (nombre_elec == 1) {
            double g_prime = gama_ts(E, angle) * (1 - gama_e(E, angle) - gama_r(E, angle));
            
            double P_e = gama_e(E, angle);
            double P_r = gama_r(E, angle);
            double P_s = g_prime * exp(-g_prime) * (1 - gama_e(E, angle) - gama_r(E, angle));
            double P_t = P_e + P_r + P_s;
            
            P_e /= P_t;
            P_r /= P_t;
            P_s /= P_t;
            
            std::vector<int> values = {0, 1, 2};
            std::vector<double> probabilities = {P_e, P_r, P_s};
            int process = Tirage(values, probabilities);
            
            if (process == 0) {
                double v_0 = sqrt((2 * E) / m);
                double c_n = cos(angle);
                double c_m = mormal.dot(vitesse);
                double c_x_rot = x_rot.dot(vitesse);
                
                Vector3d u = c_n * normal + c_m * mormal + c_x_rot * x_rot;
                
                Matrix3x3 M = Matrix3x3::Zero();
                M(0, 0) = ModifiedMat(0, 0);
                M(0, 1) = ModifiedMat(0, 1);
                M(0, 2) = ModifiedMat(0, 2);
                M(1, 0) = v_0 * u(0);
                M(1, 1) = v_0 * u(1);
                M(1, 2) = v_0 * u(2);
                M(2, 0) = ModifiedMat(2, 0);
                
                // Preserve trackID (new version)
                M(2, 2) = static_cast<float>(parentTrackID);
                
                ElectronProcess ep;
                ep.matrix = M;
                ep.processType = process; // 0: backscattering
                Resultat.push_back(ep);
                return Resultat;
            }
            
            double Energie;
            double v_0;
            
            if (process == 1) {
                std::uniform_real_distribution<double> dist(E, E0);
                Energie = dist(gen);
                v_0 = sqrt((2 * Energie) / m);
            } else { // process == 2
                v_0 = sqrt((2 * std::min(E, E0)) / m);
            }
            
            std::uniform_real_distribution<double> dist1(0.0174, 1.0);
            std::uniform_real_distribution<double> dist2(0.0, 2.0 * M_PI);
            double u0 = dist1(gen);
            double phi1 = acos(u0);
            double phi2 = dist2(gen);
            
            Vector3d u = cos(phi1) * normal + sin(phi2) * sin(phi1) * mormal + cos(phi2) * sin(phi1) * x_rot;
            
            Matrix3x3 M = Matrix3x3::Zero();
            M(0, 0) = ModifiedMat(0, 0);
            M(0, 1) = ModifiedMat(0, 1);
            M(0, 2) = ModifiedMat(0, 2);
            M(1, 0) = v_0 * u(0);
            M(1, 1) = v_0 * u(1);
            M(1, 2) = v_0 * u(2);
            M(2, 0) = ModifiedMat(2, 0);
            
            // Preserve trackID (new version)
            M(2, 2) = static_cast<float>(parentTrackID);
            
            ElectronProcess ep;
            ep.matrix = M;
            ep.processType = process; // 1: rediffused or 2: secondary
            Resultat.push_back(ep);
        } else {
            double v_0 = sqrt(((2 * std::min((E / nombre_elec), E0))) / m);
            
            std::uniform_real_distribution<double> dist1(0.0174, 1.0);
            std::uniform_real_distribution<double> dist2(0.0, 2.0 * M_PI);
            
            for (int i = 0; i < nombre_elec; i++) {
                double u0 = dist1(gen);
                double phi1 = acos(u0);
                double phi2 = dist2(gen);
                
                Vector3d u = cos(phi1) * normal + sin(phi2) * sin(phi1) * mormal + cos(phi2) * sin(phi1) * x_rot;
                
                Matrix3x3 M = Matrix3x3::Zero();
                M(0, 0) = ModifiedMat(0, 0);
                M(0, 1) = ModifiedMat(0, 1);
                M(0, 2) = ModifiedMat(0, 2);
                M(1, 0) = v_0 * u(0);
                M(1, 1) = v_0 * u(1);
                M(1, 2) = v_0 * u(2);
                M(2, 0) = ModifiedMat(2, 0);
                
                // Preserve trackID (new version)
                M(2, 2) = static_cast<float>(parentTrackID);
                
                ElectronProcess ep;
                ep.matrix = M;
                ep.processType = 2; // All secondary
                Resultat.push_back(ep);
            }
        }
    }
    
    return Resultat;
}

Matrix3x3 Physics::Recuperation(const Matrix3x3& Mat) {
    double a = c_s / 2;
    double b = Mat(1, 0);
    double c = Mat(0, 0) - x2;
    
    double t = Resolution(a, b, c);
    double y2 = Mat(0, 1) + Mat(1, 1) * t;
    double vx = Mat(1, 0) + c_s * t;
    double z2 = Mat(0, 2) + Mat(1, 2) * t;
    
    Matrix3x3 M = Matrix3x3::Zero();
    M(0, 0) = x2;
    M(0, 1) = y2;
    M(0, 2) = z2;
    M(1, 0) = vx;
    M(1, 1) = Mat(1, 1);
    M(1, 2) = Mat(1, 2);
    M(2, 0) = Mat(2, 0) + t;
    
    // Preserve trackID (new version)
    M(2, 2) = Mat(2, 2);
    
    return M;
}

Matrix3x3 Physics::Transporter1(const Matrix3x3& Mat, double t, double cts) {
    double x = Mat(0, 0) + Mat(1, 0) * t + cts * (t * t) / 2;
    double y2 = Mat(0, 1) + Mat(1, 1) * t;
    double z2 = Mat(0, 2) + Mat(1, 2) * t;
    double vx = Mat(1, 0) + cts * t;
    
    Matrix3x3 M = Matrix3x3::Zero();
    M(0, 0) = x;
    M(0, 1) = y2;
    M(0, 2) = z2;
    M(1, 0) = vx;
    M(1, 1) = Mat(1, 1);
    M(1, 2) = Mat(1, 2);
    M(2, 0) = Mat(2, 0) + t;
    M(2, 1) = Mat(2, 1) - t;
    
    // Preserve trackID (new version)
    M(2, 2) = Mat(2, 2);
    
    return M;
}

Matrix3x3 Physics::Transporter2(const Matrix3x3& Mat, double t, double cts) {
    double x = Mat(0, 0) + Mat(1, 0) * t + cts * (t * t) / 2;
    double y2 = Mat(0, 1) + Mat(1, 1) * t;
    double z2 = Mat(0, 2) + Mat(1, 2) * t;
    double vx = Mat(1, 0) + cts * t;
    
    Matrix3x3 M = Matrix3x3::Zero();
    M(0, 0) = x;
    M(0, 1) = y2;
    M(0, 2) = z2;
    M(1, 0) = vx;
    M(1, 1) = Mat(1, 1);
    M(1, 2) = Mat(1, 2);
    M(2, 0) = Mat(2, 0) + t;
    
    // Preserve trackID (new version)
    M(2, 2) = Mat(2, 2);
    
    return M;
}

void Physics::ajouter_element_trie(std::vector<Matrix3x3>& liste, const Matrix3x3& element) {
    if (!liste.empty()) {
        auto it = liste.begin();
        while (it != liste.end() && (*it)(2, 1) < element(2, 1)) {
            ++it;
        }
        liste.insert(it, element);
    } else {
        liste.push_back(element);
    }
}

bool Physics::Erreur(const Matrix3x3& Mat, int n) {
    double y_r = (Mat(0, 1) - tan(alpha) * Mat(0, 0) + tan(alpha) * x0 - (pas + dia) * n - (pas + dia) / 2);
    double z_r = Mat(0, 2);
    if (sqrt(y_r * y_r + z_r * z_r) > R + 0.1) {
        return true;
    } else {
        return false;
    }
}

std::pair<std::vector<Matrix3x3>, std::vector<Matrix3x3>> Physics::Rearrangement(
        std::vector<Matrix3x3> A1, std::vector<Matrix3x3> A2, int n, double cts) {
    
    std::vector<Matrix3x3> Emi;
    std::vector<Matrix3x3> N_Emi;
    
    for (const auto& M : A1) {
        double time = Point_de_contact2(M, n, cts);
        if (time == false) {
            N_Emi.push_back(M);
        } else if (time == true) {
            continue;
        } else {
            Matrix3x3 M_copy = M;
            M_copy(2, 1) = time;
            // Preserve trackID (new version)
            M_copy(2, 2) = M(2, 2);
            ajouter_element_trie(Emi, M_copy);
        }
    }
    
    for (const auto& M : A2) {
        if (M(0, 0) > x1) {
            N_Emi.push_back(M);
        } else {
            double time = Point_de_contact2(M, n, cts);
            if (time == false) {
                N_Emi.push_back(M);
            } else if (time == true) {
                continue;
            } else {
                Matrix3x3 M_copy = M;
                M_copy(2, 1) = time;
                // Preserve trackID (new version)
                M_copy(2, 2) = M(2, 2);
                ajouter_element_trie(Emi, M_copy);
            }
        }
    }
    
    return {Emi, N_Emi};
}

} // namespace MCPSim 