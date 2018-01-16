//
//  user_model.hh
//  
//
//  Created by Tim Dodwell on 30/11/2015.
//
//

#ifndef user_model_h
#define user_model_h

#include <random>
#include "../UQ_functions/KLFunctions.h"
#include "../UQ_functions/general.hh"
#include "linearelasticity.hh"

// ------------------------------------------------------------------------
//             Dirichlet Boundary Conditions Class
// ------------------------------------------------------------------------
template<typename GV, typename RF>
class Dirichlet_BC :
public Dune::PDELab::AnalyticGridFunctionBase<
Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
Dirichlet_BC<GV,RF> >,
public Dune::PDELab::InstationaryFunctionDefaults
{
    
private:
    int dof;
    
public:
    
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, Dirichlet_BC<GV,RF>> BaseT;
    
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    
    // Constructor
    Dirichlet_BC(const GV & gv) : BaseT(gv)
    {
    }
    
    template<typename I>
    bool isDirichlet(const I & ig,
                     const typename Dune::FieldVector<typename I::ctype, I::dimension-1> & x
                     ) const
    {
        Dune::FieldVector<double,2> xg = ig.geometry().global( x );
        bool answer = false;
        if (xg[0] < 1e-6 || xg[1] < 1e-6 || xg[0] > 1.0 - 1e-6 || xg[1] > 1.0 - 1e-6){  answer = true;  }
        return answer;
    }
    
    
    inline void evaluateGlobal(const DomainType & x, RangeType & u) const
    {
        u = 0.0;
    } // end inline function evaluateGlobal
    
    void setDof(int degree_of_freedom){
        dof = degree_of_freedom;
    }
}; // End of Dirichlet BC class


// Class which constructs
class COEFF
{
public:
    // Constructor for COEFF class
    COEFF(int numPly, int numKLmodes, double sigKL, double defProb, double ellKL, const std::vector<double>& mat)
    : numPly_(numPly), numKLmodes_(numKLmodes), sigKL_(sigKL), defectProbability_(defProb),
    ellKL_(ellKL), mat_(mat)
    {
        xi.resize(numPly_ * numKLmodes_);
        
        // Define and Resize vector for KL modes definition
        int N = std::sqrt(numKLmodes) + 1;
        freq.resize(N); lam1D.resize(N); lam2D.resize(numKLmodes);
        mode_id_i.resize(numKLmodes);   mode_id_j.resize(numKLmodes);
        
        
        double E11 = mat[0], E22 = mat[1], nu12 = mat[2], G12 = mat[3];
        
        double nu21 = nu12 * (E22 / E11);
        double factor = 1.0 - nu12 * nu21;
        
        Q[0][0] = E11 / factor; Q[0][1] = nu12 * E22 / factor; Q[0][2] = 0.0;
        Q[1][0] = Q[0][1]; Q[1][1] = E22 / factor; Q[1][2] = 0.0;
        Q[2][0] = 0.0;     Q[2][1] = 0.0;          Q[2][2] = G12;

        Dune::Timer timer;
        
        KLExpansion(N,sigKL_,ellKL_,freq);
        evaluate_eigenValues(sigKL_,ellKL_,lam1D,freq);
        construct_2D_eigenValues(lam1D, lam2D, mode_id_i, mode_id_j);
        
        std::cout << "Karman-Loeve Expansion Computed - " << timer.elapsed() << " secs" << std::endl;

    };
    
    double inline evalPhi(double x, int i){
    
        double omega = freq[i];
        
        double u = 0.0;
        if (i % 2 == 0){
            u = std::cos(omega * (x - 0.5));
            u /= sqrt(0.5 + 0.5 * (std::sin(omega) / omega));
        }
        else{
            u = std::sin(omega * (x - 0.5));
            u /= sqrt(0.5 - 0.5 * (std::sin(omega) / omega));
        }
        return u;
    }
    
    
    Dune::FieldMatrix<double,3,3> inline evaluateTensor(const Dune::FieldVector<double,2> x){
        
        double ply_thickness = 0.2;
    
        Dune::FieldMatrix<double,3,3> A(0.0);
        
        for (int k = 0; k < numPly_; k++)
        {
            double misalignment = 0.0;
            for (int j = 0; j < numKLmodes_; j++){
                misalignment += std::sqrt(lam2D[j]) * evalPhi(x[0],mode_id_i[j]) * evalPhi(x[1],mode_id_j[j]) * xi[k * numKLmodes_ + j];
            }
            
            double phi = 0.0 + misalignment;
            
            // Rotation
            
            
            Dune::FieldMatrix<double,3,3> T(0.0), TT(0.0), Qbar(0.0);
            
            double c = std::cos(phi), s = std::sin(phi);
            
            T[0][0] = c * c; T[0][1] = s * s; T[0][2] = 2.0 * s * c;
            T[1][0] = s * s; T[1][1] = c * c; T[1][2] = -2.0 * s * c;
            T[2][0] = -s * c; T[2][1] = s * c; T[2][2] = c * c - s * s;
            
            
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    TT[i][j] = T[j][i];
                }
            }
            
            Qbar = Q.rightmultiplyany(T);
            
            Qbar = Qbar.leftmultiplyany(TT);
            Qbar *= ply_thickness;
            
            A += Qbar;
        }
        
        return A;
    
    }

    
    void inline user_random_field()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> randu(0.0, 1.0);
        std::normal_distribution<double> randn(0.0,sigKL_);
        
        std::fill(xi.begin(), xi.end(), 0.0);
        
        for (int i = 0; i < numPly_; i++){
            if (randu(gen) > defectProbability_){
                for (int j = 0; j < numKLmodes_; j++){
                    xi[j + numKLmodes_ * i] = randn(gen);
                }
            }
        }
        
    } // end user_random_field
    
    
private:
    int numPly_;
    int numKLmodes_;
    Dune::FieldMatrix<double,3,3> Q;
    double sigKL_, defectProbability_, ellKL_;
    std::vector<double> xi;
    std::vector<double> freq, lam1D, lam2D, mat_;
    std::vector<int> mode_id_i, mode_id_j;
};


template<typename GRID>
class MODEL{
    
public:
    
    int numPly = 8;
    const int numMode1D = 4;
    int numKLmodes = 10;
    double defectProbability = 1.0;
    double sigKL = 1.0;
    double ellKL = 1.0;
    
    std::vector<double> mat;
    
    
    const int numRandomVector = 8 * numKLmodes;


    
    
    // Constructor for MODEL CLASS
    MODEL(GRID& grid):grid_(grid){
    
        mat.resize(4);
        mat[0] = 120.0;
        mat[1] = 9.0;
        mat[2] = 0.25;
        mat[3] = 5.0;
    
    };
    
    double inline user_fem_driver(int l, COEFF& z) const{
        
        typedef typename GRID::LevelGridView LGV;
        LGV gv = grid_.levelGridView(l);
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> randu(0.0, 1.0);
        std::normal_distribution<double> randn(0.0,sigKL);
        
        typedef double RF;
        const int dim = LGV::Grid::dimension;
        typedef typename LGV::Grid::ctype Coord;
        
        typedef Dune::PDELab::QkLocalFiniteElementMap<LGV,Coord,RF,1> FEM;
        FEM fem(gv);
        
        typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
        CON con;
        
        typedef Dune::PDELab::ISTLVectorBackend<> VectorBackend;
        
        typedef Dune::PDELab::GridFunctionSpace<LGV, FEM, CON, VectorBackend> SCALAR_GFS;
        SCALAR_GFS dispU1(gv,fem,con); dispU1.name("U1");
        SCALAR_GFS dispU2(gv,fem,con); dispU2.name("U2");
        
        
        typedef Dune::PDELab::CompositeGridFunctionSpace <VectorBackend,Dune::PDELab::LexicographicOrderingTag, SCALAR_GFS, SCALAR_GFS> GFS;
        GFS gfs(dispU1, dispU2);
        
        // Make constraints map and initialize it from a function
        typedef typename GFS::template ConstraintsContainer<RF>::Type C;
        C cg;
        cg.clear();
        
        typedef Dirichlet_BC<LGV,RF> BC;
        BC U1_cc(gv), U2_cc(gv);
        U1_cc.setDof(1);
        U2_cc.setDof(2);
        
        typedef Dune::PDELab::CompositeConstraintsParameters<BC,BC>
        Constraints;
        Constraints constraints(U1_cc,U2_cc);
        
        Dune::PDELab::constraints(constraints,gfs,cg);
        
        double Q = randn(gen);
        
        return Q;
    }
    
    
private:
    
    GRID& grid_;


    
    
};









#endif /* user_model_h */
