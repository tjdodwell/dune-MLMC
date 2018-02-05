#ifndef RandomField_Exp2D_h
#define RandomField_Exp2D_h

#include <random>
#include "KLFunctions.h"
#include "general.hh"

class RandomField
{

public:

  int level;

  RandomField(){

    numKLmodes = config.get<int>("RandomField.numKLmodes",100);

    L = config.get<double>("Domain.L",4.0);

    sigKL = config.get<double>("RandomField.sigKL",0.035); // rad  - Standard Deviation of Fibre Misalignment
    ellx = config.get<double>("RandomField.ellx",1.6); // mm    - Correlation Length Scale Parallel to the Fibres
    elly = config.get<double>("RandomField.elly",0.43); // mm    - Correlation Length Scale Parallel to the Fibres


    xi.resize(numKLmodes); // Resize Random Field Container to Stochastic Dimension

    int N = std::sqrt(numKLmodes) + 1;


    freqx.resize(N);  freqy.resize(N);
    lam1Dx.resize(N); lam1Dy.resize(N);

    

    rootFinder(N, ellx / L, freqx); 
    rootFinder(N, elly / L, freqy);

    evaluate_eigenValues(ellx/L,1.0,lam1Dx,freqx);
    evaluate_eigenValues(elly/L,1.0,lam1Dy,freqy);

    lam3D.resize(numKLmodes);
    mode_id_i.resize(numKLmodes); mode_id_j.resize(numKLmodes);

    construct_2D_eigenValues(lam1Dx, lam1Dy, lam3D, mode_id_i, mode_id_j);

    //

    // Derived Properties

    Vf = config.get<double>("Material.Vf",0.5); // Vf     - Fibre Volume Fraction
    d = config.get<double>("Material.d",0.0071); // mm  - Diameter of a  Fibre (mm)

    E1 = config.get<double>("Material.E1",128.0); // GPa - Longitudial Stiffness
    E2 = config.get<double>("Material.E2",9.25); // GPa   - Tranverse Stiffness
    G12 = config.get<double>("Material.G12",5.1); // GPa    - Inplane shear stiffness
    nu12 = config.get<double>("Material.nu12",0.35);//        - Poisson Ratio


    double nu21 = nu12 * (E2 / E1);
    double factor =  1.0 - nu12  * nu21;

    Dune::FieldMatrix<double,3,3> Q(0.0);

    Q[0][0] = E1 / factor;
    Q[0][1] = ( nu12 * E2 ) / factor;
    Q[1][0] = Q[0][1];
    Q[1][1] = E2 / factor;
    Q[2][2] = G12;

    CB = E1 * (d * d) / 16.0;

    D = 0.0;
    for (int i = 0; i < 3; i++){for (int j = 0; j < 3; j++){ D[i][j] = Q[i][j]; }}
    D[2][3] = Q[2][2];
    D[3][2] = Q[2][2];
    D[3][3] = Vf * E1 /(2.0 * (1.0 - 0.2)) + (1.0 - Vf) * G12;
    D[4][4] = CB * (D[3][3] - D[2][2]) / D[3][3];

  }



  void inline setLevel(int i){ level = i;  }

  void inline user_random_field(bool isTest = false){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> randn(0.0,1.0);
    std::fill(xi.begin(), xi.end(), 0.0);
    
    for (int j = 0; j < numKLmodes; j++){
            xi[j] = randn(gen);
    }
  }

   double inline evalPhi(double x, int i, bool isX = true){
      double phi = 0.0;
      double omega;
      x /= L;

      double l = 0.0;

      if (isX){
        l = ellx/L;
        omega =   freqx[i];
      }
      else{
        l = elly/L;
        omega =   freqy[i];
      }

      double norm = std::sin(2.0*omega)*(0.25*(l*l*omega - 1.0/omega)) - 0.5 * l * std::cos(2.0 * omega) + 0.5 * ( 1 + l + l*l*omega*omega);
      norm = 1.0 / std::sqrt(norm);

      phi = norm * (std::sin(x * omega) + l * omega * std::cos(x * omega) );
      return phi;
    }

    Dune::FieldMatrix<double,6,6> evaluateMatrix(Dune::FieldVector<double,2>& x){
        // Compute Misalignment at a point
        double phi = 0.0;
        int R = std::min(numKLmodes,50 + 50 * level);
        for (int j = 0; j < R; j++){
            phi += sigKL * std::sqrt(lam3D[j]) * evalPhi(x[0],mode_id_i[j],true) * evalPhi(x[1],mode_id_j[j],false) * xi[j];
        }
        // Compute Rotation matrix
        Dune::FieldMatrix<double,6,6> T(0.0), TT(0.0);
        double c = std::cos(phi);
        double s = std::sin(phi);

        T[0][0] = c * c;  T[0][1] = s * s;  T[0][2] = c * s;  T[0][3] = c * s;
        T[1][0] = s * s;  T[1][1] = c * c;  T[1][2] = -c * s; T[1][3] = -c * s;
        T[2][0] = -c * s; T[2][1] = c * s;  T[2][2] = c * c;  T[2][3] = -s * s;
        T[3][0] = -c * s; T[3][1] = c * s;  T[3][2] = -s * s; T[3][3] =  c * c;
        T[4][4] = c;  T[4][5] = s;
        T[5][4] = -s; T[5][5] = c;

        for (int i = 0; i < 6; i++){
          for (int j = 0; j < 6;  j++){
            TT[i][j] = T[j][i];
          }
        }

        // Apply rotation to Q  i.e  Qs = tran(T) * Q  * T

        Dune::FieldMatrix<double,6,6> Qs(0.0);

        Qs = D.rightmultiplyany(T);

        Qs = Qs.leftmultiplyany(TT);

        return Qs; // Return Rotated Tensor
    }

    Dune::FieldMatrix<double,6,6> evaluateLocalMatrix(){
      return D;
    }

    double evaluateMisalignment(Dune::FieldVector<double,2>& x){
      double phi = 0.0;
      int R = std::min(numKLmodes,50 + 50 * level);
      for (int j = 0; j < R; j++){
          phi += sigKL * std::sqrt(lam3D[j]) * evalPhi(x[0],mode_id_i[j],true) * evalPhi(x[1],mode_id_j[j],false) * xi[j];
      }
      return phi;
    }

    Dune::FieldMatrix<double,6,6> evaluateRotation(Dune::FieldVector<double,2>& x){
      // Compute Misalignment at a point
      double phi = 0.0;
      int R = std::min(numKLmodes,50 + 50 * level);
      for (int j = 0; j < R; j++){
          phi += std::sqrt(lam3D[j]) * evalPhi(x[0],mode_id_i[j],true) * evalPhi(x[1],mode_id_j[j],false) * xi[j];
      }
      // Compute Rotation matrix
      Dune::FieldMatrix<double,6,6> T(0.0);
      double c = std::cos(phi);
      double s = std::sin(phi);

      T[0][0] = c * c;  T[0][1] = s * s;  T[0][2] = c * s;  T[0][3] = c * s;
      T[1][0] = s * s;  T[1][1] = c * c;  T[1][2] = -c * s; T[1][3] = -c * s;
      T[2][0] = -c * s; T[2][1] = c * s;  T[2][2] = c * c;  T[2][3] = -s * s;
      T[3][0] = -c * s; T[3][1] = c * s;  T[3][2] = -s * s; T[3][3] =  c * c;
      T[4][4] = c;  T[4][5] = s;
      T[5][4] = -s; T[5][5] = c;

      return T; // Return rotation matrix
    }

private:



  double L, sigKL, ellx, elly;

  double Vf, d, E1, E2, G12, nu12,tau_y, PatchSize,CB;

  int numKLmodes;

  std::vector<double> xi, freqx, freqy, lam1Dx, lam1Dy, lam3D;
  std::vector<int> mode_id_i, mode_id_j;

    
  std::vector<Dune::FieldMatrix<double,6,6>> Cijkl;
  Dune::FieldMatrix<double,6,6> D;



};

#endif
