#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

std::random_device rd;
std::mt19937 gen(rd());

const double diff_pot = 600;
const double m = 9.1;
const double x0 = 500;
const double x1 = 900;
const double x2 = 1500;
const double c_c = (1.6 * diff_pot) / ((x1 - x0) * m);
const double c_e = 0.05621;
const double c_s = (1.6 * diff_pot) / ((x1 - x0) * m);
const double pas = 2;
const double dia = 10;
const double R = dia / 2;
const double alpha = 0.13;
const double limite = 599;
const double E0 = 5;
const double Resistance = pow(10, 8);
const double I_strip = diff_pot / Resistance;
const double q = 1.6 * pow(10, -7);

const double pi = M_PI;

const double P_inf = 0.05;
const double P_1e = 0.1;
const double W = 60;
const double Ee = 0;
const double p = 1;
const double P_1r = 0.5;
const double r = 1;
const double Er = 30;
const double e1 = 0.26;
const double e2 = 2;
const double r1 = 0.26;
const double r2 = 2;
const double t1 = 0.25;
const double t2 = 0.8;
const double t3 = 0.25;
const double t4 = 1;
const double gama_ts0 = 3.5;
const double E_ts = 260;
const double s = 1.54;

typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
typedef Eigen::Vector3d Vector3d;

double comb(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    double c = 1.0;
    for (int i = 0; i < k; i++) {
        c = c * (n - i) / (i + 1);
    }
    return c;
}

double gama_e0(double E) {
    return P_inf + (P_1e - P_inf) * exp(-pow(fabs(E - Ee) / W, p) / p);
}

double gama_r0(double E) {
    return P_1r * (1 - exp(-pow(E / Er, r)));
}

double gama_e(double E, double teta) {
    return gama_e0(E) * (1 + e1 * (1 - pow(cos(teta), e2)));
}

double gama_r(double E, double teta) {
    return gama_r0(E) * (1 + r1 * (1 - pow(cos(teta), r2)));
}

double D(double x, double s) {
    return (s * x) / (s - 1 + pow(x, s));
}

double gama_teta(double teta) {
    return gama_ts0 * (1 + t1 * (1 - pow(cos(teta), t2)));
}

double E_teta(double teta) {
    return E_ts * (1 + t3 * (1 - pow(cos(teta), t4)));
}

double gama_ts(double E, double teta) {
    return gama_teta(teta) * D((E / E_teta(teta)), s);
}

double gama_totale(double E) {
    return gama_ts(E, 1.39) / (1 - gama_e(E, 1.39) + gama_r(E, 1.39));
}

// 확률 관련 함수
double P_prime(double para, int M, int n) {
    return comb(M, n) * pow(para, n) * pow(1 - para, M - n);
}

double P(double gamae, double gamar, double para, int M, int n) {
    if (n == 0) {
        return (1 - gamae - gamar) * P_prime(para, M, n);
    } else if (n == 1) {
        return (1 - gamae - gamar) * P_prime(para, M, n) + gamae + gamar;
    } else {
        return (1 - gamae - gamar) * P_prime(para, M, n);
    }
}

// 난수 기반으로 값 선택
int Tirage(const std::vector<int>& valeurs, const std::vector<double>& probabilites) {
    std::discrete_distribution<int> dist(probabilites.begin(), probabilites.end());
    return valeurs[dist(gen)];
}

// 전자 생성 함수
int Generate(double E, double teta) {
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

// 방정식 관련 함수
double Delta(double a, double b, double c) {
    return b * b - 4 * a * c;
}

double Resolution(double a, double b, double c) {
    if (Delta(a, b, c) < 0) {
        return false;  // 실수 해가 없음
    } else {
        double racine_delta = sqrt(Delta(a, b, c));
        return (-b + racine_delta) / (2 * a);
    }
}

std::vector<double> Resolution_bis(double a, double b, double c) {
    std::vector<double> solutions;
    if (Delta(a, b, c) < 0) {
        solutions.push_back(false);  // 실수 해가 없음
        return solutions;
    } else {
        double racine_delta = sqrt(Delta(a, b, c));
        solutions.push_back((-b + racine_delta) / (2 * a));
        solutions.push_back((-b - racine_delta) / (2 * a));
        return solutions;
    }
}

bool Check_iscts(int iscts) {
    return iscts == 1;
}

Vector3d Rot(const Vector3d& v, double teta) {
    Eigen::Matrix3d rotation;
    rotation << cos(teta), -sin(teta), 0,
                sin(teta), cos(teta), 0,
                0, 0, 1;
    return rotation * v;
}

std::vector<double> generate_cosine_angle(int num_samples) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<double> angles(num_samples);
    for (int i = 0; i < num_samples; i++) {
        angles[i] = asin(dist(gen));
    }
    return angles;
}

// 메인 시뮬레이션 함수들
Matrix3x3 Pho_ele(double Energie_photon) {
    double norme_vitesse = sqrt(2 * Energie_photon / m);
    
    Matrix3x3 Mpe = Matrix3x3::Zero();
    Mpe(0, 1) = 6;
    Mpe(1, 0) = norme_vitesse;
    Mpe(1, 1) = 0;
    
    return Mpe;
}

Matrix3x3 premiere_arrivee(Matrix3x3 Matrice_photo_electron) {
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

std::pair<bool, int> Check_if_hit(const Matrix3x3& Matrice_arrive) {
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

// 다항식의 근을 구하는 함수 - 4차 방정식 해법
std::vector<double> solve_quartic(double a4, double a3, double a2, double a1, double a0) {
    std::vector<double> roots;
    
    // a4x^4 + a3x^3 + a2x^2 + a1x + a0 = 0
    // 표준형태로 정규화: x^4 + px^3 + qx^2 + rx + s = 0
    double p = a3 / a4;
    double q = a2 / a4;
    double r = a1 / a4;
    double s = a0 / a4;
    
    // 동반 행렬 구성 (Companion matrix)
    Eigen::Matrix4d companion;
    companion << 0, 0, 0, -s,
                 1, 0, 0, -r,
                 0, 1, 0, -q,
                 0, 0, 1, -p;
    
    // 고유값 계산 (고유값 = 다항식의 근)
    Eigen::EigenSolver<Eigen::Matrix4d> solver(companion);
    
    // 복소수 근 추출
    Eigen::Vector4cd eigenvalues = solver.eigenvalues();
    
    // 실수 근만 필터링
    for (int i = 0; i < 4; i++) {
        if (std::abs(eigenvalues(i).imag()) < 1e-10) {
            roots.push_back(eigenvalues(i).real());
        }
    }
    
    return roots;
}

double Point_de_contact2(const Matrix3x3& Mat, int n, double cts) {
    double Mat00 = Mat(0, 0), Mat01 = Mat(0, 1), Mat02 = Mat(0, 2);
    double Mat10 = Mat(1, 0), Mat11 = Mat(1, 1), Mat12 = Mat(1, 2);
    
    // 4차 다항식 계수
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
    
    // 수치해석으로 4차 방정식 근 구하기
    // 실제 구현에서는 Eigen, GSL 등의 라이브러리 사용 권장
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

std::vector<double> generer_angles(int n, double teta1, double teta2) {
    std::vector<double> angles(n);
    std::uniform_real_distribution<double> dist(teta1, teta2);
    
    for (int j = 0; j < n; j++) {
        angles[j] = dist(gen);
    }
    
    return angles;
}

std::vector<Matrix3x3> emi_sec(const Matrix3x3& Mat, int n) {
    std::vector<Matrix3x3> Resultat;
    
    // 에너지 계산
    double E = 0.5 * m * (pow(Mat(1, 0), 2) + pow(Mat(1, 1), 2) + pow(Mat(1, 2), 2));
    
    // teta 계산
    double y_r = Mat(0, 1) - tan(alpha) * Mat(0, 0) + tan(alpha) * x0 - (pas + dia) * n - (pas + dia) / 2;
    double z_r = Mat(0, 2);
    
    double teta = atan2(z_r, y_r);
    if (teta < 0) {
        teta += 2 * pi;
    }
    
    // 행렬 생성
    Matrix3x3 ModifiedMat = Mat;
    ModifiedMat(0, 1) = R * cos(teta) + tan(alpha) * (Mat(0, 0) - x0) + (pas + dia) * n + (pas + dia) / 2;
    ModifiedMat(0, 2) = R * sin(teta);
    
    // 벡터 계산
    Vector3d e_r(0, cos(teta), sin(teta));
    Vector3d e_o(0, -sin(teta), cos(teta));
    Vector3d z(0, 0, 1);
    Vector3d x(1, 0, 0);
    
    // 회전 행렬 적용
    Vector3d normal = Rot(-e_r, alpha);
    Vector3d x_rot = Rot(x, alpha);
    Vector3d mormal = normal.cross(x_rot);
    
    // 속도 벡터와 각도 계산
    double norme = sqrt(pow(Mat(1, 0), 2) + pow(Mat(1, 1), 2) + pow(Mat(1, 2), 2));
    Vector3d vitesse(Mat(1, 0) / norme, Mat(1, 1) / norme, Mat(1, 2) / norme);
    double scalaire = normal.dot(-vitesse);
    double angle = acos(scalaire);
    
    if (angle > 89 * pi / 180) {
        angle = 89 * pi / 180;
    }
    
    // 전자 수 생성
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
                
                Resultat.push_back(M);
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
            std::uniform_real_distribution<double> dist2(0.0, 2.0 * pi);
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
            
            Resultat.push_back(M);
        } else {
            double v_0 = sqrt(((2 * std::min((E / nombre_elec), E0))) / m);
            
            std::uniform_real_distribution<double> dist1(0.0174, 1.0);
            std::uniform_real_distribution<double> dist2(0.0, 2.0 * pi);
            
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
                
                Resultat.push_back(M);
            }
        }
    }
    
    return Resultat;
}

Matrix3x3 Recuperation(const Matrix3x3& Mat) {
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
    
    return M;
}

Matrix3x3 Transporter1(Matrix3x3 Mat, double t, double cts) {
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
    
    return M;
}

Matrix3x3 Transporter2(Matrix3x3 Mat, double t, double cts) {
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
    
    return M;
}

void ajouter_element_trie(std::vector<Matrix3x3>& liste, const Matrix3x3& element) {
    if (!liste.empty()) {
        auto it = liste.begin();
        while (it != liste.end() && 
               (*it)(2, 1) < element(2, 1)) {
            ++it;
        }
        liste.insert(it, element);
    } else {
        liste.push_back(element);
    }
}

bool Erreur(const Matrix3x3& Mat, int n) {
    double y_r = (Mat(0, 1) - tan(alpha) * Mat(0, 0) + tan(alpha) * x0 - (pas + dia) * n - (pas + dia) / 2);
    double z_r = Mat(0, 2);
    if (sqrt(y_r * y_r + z_r * z_r) > R + 0.1) {
        return true;
    } else {
        return false;
    }
}

std::pair<std::vector<Matrix3x3>, std::vector<Matrix3x3>> Rearrangement(
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
                ajouter_element_trie(Emi, M_copy);
            }
        }
    }
    
    return {Emi, N_Emi};
}

std::vector<Matrix3x3> simulate(double E, double c) {
    double cts = c;
    Matrix3x3 Mat = Pho_ele(E);  // 광전자 생성
    Mat = premiere_arrivee(Mat);  // 광전자 MCP 입구 도달
    
    auto check = Check_if_hit(Mat);  // 포어 진입 여부 확인
    
    // std::cout << "Check_if_hit: " << (check.first ? "Success to arrival" : "Fail to arrival") << ", Position: " << Mat(0, 1) << std::endl;

    std::vector<Matrix3x3> Resultat_final;  // 최종 결과 저장
    
    if (check.first && -limite < Mat(0, 1) && Mat(0, 1) < limite) {
        int condition = 0;
        std::vector<Matrix3x3> Ensemble_emi;  // 충돌할 전자 목록
        
        // 초기 전자 3개 생성
        for (int w = 0; w < 3; w++) {
            Matrix3x3 M = Matrix3x3::Zero();
            M(0, 0) = Mat(0, 0) - 0.1 * w;
            M(0, 1) = Mat(0, 1);
            M(0, 2) = Mat(0, 2);
            M(1, 0) = Mat(1, 0);
            M(1, 1) = Mat(1, 1);
            M(1, 2) = Mat(1, 2);
            M(2, 0) = Mat(2, 0);
            
            double time_1_hit = Point_de_contact2(M, check.second, cts);
            M(2, 1) = time_1_hit;
            Ensemble_emi.push_back(M);
        }

        // std::cout << Ensemble_emi.size() << " initial electrons are generated." << std::endl;
        
        std::vector<Matrix3x3> Ensemble_non_emi;  // 충돌하지 않을 전자 목록
        int Etat = 0;  // 전기장 조절 상태
        double instant_0 = 0;
        
        while (condition == 0) {
            if (!Ensemble_emi.empty()) {
                // 가장 빠른 충돌 처리
                double I = 0;
                Matrix3x3 lead_M = Ensemble_emi.front();
                double temps = lead_M(2, 1);
                double instant = lead_M(2, 0);
                lead_M = Transporter1(lead_M, temps, cts);
                std::vector<Matrix3x3> Elec_secondaire = emi_sec(lead_M, check.second);
                Ensemble_emi.erase(Ensemble_emi.begin());  // 처리한 첫번째 전자 제거
                
                // 비충돌 전자 업데이트
                if (!Ensemble_non_emi.empty()) {
                    std::vector<int> indice;
                    for (int j = 0; j < Ensemble_non_emi.size(); j++) {
                        Matrix3x3 M_passage = Transporter2(Ensemble_non_emi[j], temps, cts);
                        
                        if (M_passage(0, 0) >= x2) {  // 양극 도달
                            Matrix3x3 M_recuperated = Recuperation(Ensemble_non_emi[j]);
                            Resultat_final.push_back(M_recuperated);
                            indice.push_back(j);
                        } else {  // 계속 이동 중
                            Ensemble_non_emi[j] = M_passage;
                            if (Ensemble_non_emi[j](0, 0) < x1) {  // MCP 내부에 있는 경우
                                I += q * Ensemble_non_emi[j](1, 0) / (x1 - x0);
                            }
                        }
                    }
                    
                    // 큰 인덱스부터 제거
                    std::sort(indice.begin(), indice.end(), std::greater<int>());
                    for (int idx : indice) {
                        Ensemble_non_emi.erase(Ensemble_non_emi.begin() + idx);
                    }
                }
                
                // 대기 중인 충돌 전자 업데이트
                for (int j = 0; j < Ensemble_emi.size(); j++) {
                    Ensemble_emi[j] = Transporter1(Ensemble_emi[j], temps, cts);
                    I += q * Ensemble_emi[j](1, 0) / (x1 - x0);
                }
                
                // 새로 생성된 2차 전자 처리
                for (int k = 0; k < Elec_secondaire.size(); k++) {
                    double time = Point_de_contact2(Elec_secondaire[k], check.second, cts);
                    
                    if (time == false) {  // MCP 빠져나감
                        Ensemble_non_emi.push_back(Elec_secondaire[k]);
                    } else if (time == true) {  // 충돌하지 않음
                        continue;
                    } else {  // 충돌할 예정
                        Matrix3x3 temp = Elec_secondaire[k];
                        temp(2, 1) = time;
                        ajouter_element_trie(Ensemble_emi, temp);
                    }
                }
                
                // 전기장 동적 조절 메커니즘
                if (Etat == 0) {  // 정상 상태
                    if (I >= 0.05 * I_strip) {  // 전류 임계값 초과
                        instant_0 = instant;
                        cts = cts * 0.8;  // 전기장 20% 감소
                        auto result = Rearrangement(Ensemble_emi, Ensemble_non_emi, check.second, cts);
                        Ensemble_emi = result.first;
                        Ensemble_non_emi = result.second;
                        Etat = 1;  // 감소 상태로 전환
                    }
                } else {  // 감소 상태
                    if (instant >= instant_0 + 5 && I < 0.05 * I_strip) {  // 5ps 후 전류 안정
                        instant_0 = instant;
                        cts = cts / 0.8;  // 전기장 복원
                        auto result = Rearrangement(Ensemble_emi, Ensemble_non_emi, check.second, cts);
                        Ensemble_emi = result.first;
                        Ensemble_non_emi = result.second;
                        
                        if (std::abs(cts - c) < 1e-10) {  // 원래 값으로 복원 (부동소수점 비교)
                            Etat = 0;
                        }
                    }
                    
                    if (instant >= instant_0 + 5 && I > 0.05 * I_strip) {  // 5ps 후에도 임계값 초과
                        cts = cts * 0.8;  // 전기장 추가 감소
                        instant_0 = instant;
                        auto result = Rearrangement(Ensemble_emi, Ensemble_non_emi, check.second, cts);
                        Ensemble_emi = result.first;
                        Ensemble_non_emi = result.second;
                    }
                }
            } else {  // Ensemble_emi가 비어있는 경우
                std::vector<int> indice2;
                for (int j = 0; j < Ensemble_non_emi.size(); j++) {
                    Matrix3x3 result = Recuperation(Ensemble_non_emi[j]);
                    Resultat_final.push_back(result);
                    indice2.push_back(j);
                }
                
                // 인덱스 큰 것부터 제거
                std::sort(indice2.begin(), indice2.end(), std::greater<int>());
                for (int idx : indice2) {
                    Ensemble_non_emi.erase(Ensemble_non_emi.begin() + idx);
                }
            }
            
            // 종료 조건 체크
            if (Ensemble_emi.empty() && Ensemble_non_emi.empty()) {
                condition = 1;
            }
        }
    }
    
    return Resultat_final;
}

int main(int argc, char* argv[]) {

    std::string job_id = argv[1];
    int index = std::stoi(argv[2]);

    std::vector<Matrix3x3> Rs = simulate(4.21, c_c);
        
    std::ofstream outfile("/home/jangh/MCPSim/pbs_output/FromCpp/400µm600V_3.5_0.5_saturation_0.8_each_5ps_" + job_id + "_" + std::to_string(index) + ".txt");

    if (outfile.is_open()) {
        outfile << "Total number of secondary electrons: " << Rs.size() << std::endl;
        
        for (int i = 0; i < Rs.size(); i++) {
            outfile << Rs[i](0, 0) << " " << Rs[i](0, 1) << " " << Rs[i](0, 2) << " " << Rs[i](1, 0) << " " << Rs[i](1, 1) << " " << Rs[i](1, 2) << " " << Rs[i](2, 0) << " " <<  Rs[i](2, 1) << " " << Rs[i](2, 2) << std::endl;
        }
        outfile.close();
    }
    
    return 0;
}
