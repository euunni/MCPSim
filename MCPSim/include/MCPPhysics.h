#pragma once

#include "MCPConfig.h"
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <utility>
#include <cmath>

namespace MCPSim {

using Matrix3x3 = Eigen::Matrix<double, 3, 3>;
using Vector3d = Eigen::Vector3d;

struct ElectronProcess {
    Matrix3x3 matrix;  
    int processType;   
};

class Physics {
public:
    static std::random_device rd;
    static std::mt19937 gen;

    static void initialize();
    
    static double comb(int n, int k);
    static std::vector<double> solve_quartic(double a4, double a3, double a2, double a1, double a0);
    static double gama_e0(double E);
    static double gama_r0(double E);
    static double gama_e(double E, double teta);
    static double gama_r(double E, double teta);
    static double D(double x, double s);
    static double gama_teta(double teta);
    static double E_teta(double teta);
    static double gama_ts(double E, double teta);
    static double P_prime(double para, int M, int n);
    static double P(double gamae, double gamar, double para, int M, int n);
    static int Tirage(const std::vector<int>& values, const std::vector<double>& probabilities);
    static int Generate(double E, double teta);
    static Vector3d Rot(const Vector3d& v, double teta);
    static Matrix3x3 Pho_ele(double Energie_photon);
    static double Delta(double a, double b, double c);
    static double Resolution(double a, double b, double c);
    static Matrix3x3 premiere_arrivee(Matrix3x3 Matrice_photo_electron);
    static std::pair<bool, int> Check_if_hit(const Matrix3x3& Matrice_arrive);
    static double Point_de_contact2(const Matrix3x3& Mat, int n, double cts, double alpha, double x0, double x1, double R, double dia, double pas);
    static Matrix3x3 Recuperation(const Matrix3x3& Mat);
    static Matrix3x3 RecuperationTo(const Matrix3x3& Mat, double x_target);
    static Matrix3x3 Transporter1(const Matrix3x3& Mat, double t, double cts, double alpha);
    static Matrix3x3 Transporter2(const Matrix3x3& Mat, double t, double cts, double alpha);
    static void ajouter_element_trie(std::vector<Matrix3x3>& liste, const Matrix3x3& element);
    static bool Erreur(const Matrix3x3& Mat, int n, double alpha, double x0, double pas, double dia, double R);
    static std::pair<std::vector<Matrix3x3>, std::vector<Matrix3x3>> Rearrangement(
        std::vector<Matrix3x3> A1, std::vector<Matrix3x3> A2, int n, double cts, double alpha, double x0, double x1, double R, double dia, double pas);
    static std::vector<ElectronProcess> emi_sec(const Matrix3x3& Mat, int n, double alpha, double x0, double R, double dia, double pas, double m, double E0);
    static Matrix3x3 TransportGap(const Matrix3x3& Mat, double t, double cts);

private:
    static double x0; 
    static double x1;
    static double x2;
    static double q;
    static double R;
    static double limite;
    static double P_inf;      // probability of secondary emission lower limit
    static double P_1e;       // probability of elastic scattering
    static double P_1r;       // probability of rear scattering
    static double Ee;         // elastic scattering energy
    static double Er;         // rear scattering energy
    static double W;          // energy width
    static double p;          // probability distribution index
    static double r;          // rear scattering index
    static double e1, e2;     // elastic scattering coefficient
    static double r1, r2;     // rear scattering coefficient
    static double gama_ts0;   // secondary emission coefficient
    static double t1, t2;     // angle dependent coefficient 1
    static double E_ts;       // secondary emission threshold energy
    static double t3, t4;     // angle dependent coefficient 2
    static double s;          // energy distribution index
    static double m;          // electron mass
    static double dia;        // channel diameter
    static double pas;        // channel gap
    static double E0;         // reference energy
    static double diff_pot;   // potential difference
    static double c_e;        // electric field coefficient 1
    static double c_s;        // electric field coefficient 2
    static double c_c;        // electric field coefficient 3 (calculated value)
    static double Resistance; // resistance
};

} // namespace MCPSim 