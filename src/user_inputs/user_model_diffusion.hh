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
#include "diffusion.hh"
#include "AdaptGrid.hh"


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

  bool isTest;


public:

    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, Dirichlet_BC<GV,RF>> BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    // Constructor
    Dirichlet_BC(const GV & gv,bool Test = false) : BaseT(gv), isTest(Test)
    {
    }

    template<typename I>
    bool isDirichlet(const I & ig,
                     const typename Dune::FieldVector<typename I::ctype, I::dimension-1> & x
                     ) const
    {
        Dune::FieldVector<double,3> xg = ig.geometry().global( x );
        bool answer = true;
        //if (xg[0] < 1e-6 || xg[1] < 1e-6 || xg[0] > 1.0 - 1e-6 || xg[1] > 1.0 - 1e-6){  answer = true;  }
        return answer;
    }


    inline void evaluateGlobal(const DomainType & x, RangeType & u) const
    {


          u = 0.0;

    } // end inline function evaluateGlobal


}; // End of Dirichlet BC class

#include "Ltwo.hh"

#include "QuantityofInterest.hh"
#include "goallop.hh"

#include<dune/common/fvector.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

//! \brief Parameter class selecting boundary conditions
class BCTypeParam
: public Dune::PDELab::DirichletConstraintsParameters
{
public:
    //! Test whether boundary is Dirichlet-constrained
    template<typename I>
    bool isDirichlet(const I & intersection,
                     const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                     ) const
    {
        Dune::FieldVector<typename I::ctype, I::dimension>
        xg = intersection.geometry().global( coord );

        return true;  // Dirichlet b.c. on all other boundaries
    }

};

// Class which constructs

// Class which constructs
class COEFF
{
public:

    // Constructor for COEFF class
    COEFF(double L,int numKLmodes, double sigKL, double ellKL, const std::vector<double>& mat)
    : numKLmodes_(numKLmodes), sigKL_(sigKL), ellKL_(ellKL), param_(mat)
    {

       const int dim = 2;

        xi.resize(numKLmodes_);

        // Define and Resize vector for KL modes definition
        int N = std::pow(numKLmodes,1.0/dim) + 1;
        freq.resize(N); lam1D.resize(N); lam.resize(numKLmodes);
        mode_id_i.resize(numKLmodes);
        mode_id_j.resize(numKLmodes);
        mode_id_k.resize(numKLmodes);

        KLExpansion(N,0.5 * L,sigKL_,ellKL_,freq);
        evaluate_eigenValues(ellKL_,lam1D,freq);

        construct_2D_eigenValues(lam1D, lam, mode_id_i, mode_id_j);



    };

    double inline evaluateF(Dune::FieldVector<double,2>& x, bool isDual = false) const{


      double f = 0.0;


      if(isDual){

        Dune::FieldVector<double,2> y(0.75); double radius = 0.05;

        double d = 0.0;
        for (int i = 0; i < 2; i++){
          d += (x[i] - y[i]) * (x[i] - y[i]);
        }
        d = std::sqrt(d);

        if (d < radius){ f = 1.0;}

      }
      else{

      if (isTest){

        double alpha = 5.0;

        double u = x[0] * x[1] * (1.0 - x[0]) * (1.0  - x[1]) * std::exp(alpha * x[0]);

        double dudx = (1 - 2 * x[0]) * (x[1]  - x[1] * x[1]) * std::exp(alpha * x[0]) + alpha * u;


        f -=  -2.0 * (x[1] - x[1] * x[1]) * std::exp(alpha * x[0]);
        f -= alpha * (1 - 2.0 * x[0]) * (x[1] - x[1] * x[1]) * std::exp(alpha * x[0]);
        f -= alpha * dudx;
        f -= -2.0 * x[0] * (1  - x[0]) * std::exp(alpha * x[0]);

      }

      else { f = 1.0; }


    }

      return f;

    }

    double inline evalPhi(double x,int i, double L){

        double phi = 0.0;

        double omega = freq[i];

        x -= L;

        double tmp = sqrt(L + std::sin(2 * omega * L) / (2 * omega));

        if (i % 2 == 0){ phi = std::cos(omega * x) / tmp;}
        else { phi = std::sin(omega * x) / tmp; }

        return phi;

    }


    template <typename GV>
    void inline computeTensor(GV& gv, double L,bool Test = false)
    {
      isTest = Test;

        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

        const typename GV::IndexSet& is(gv.indexSet());

        Kij.resize(is.size(0));


        for (ElementIterator eg = gv.template begin<0>(); eg!=gv.template end<0>(); ++eg)
        { // loop through each element
            int id = is.index(*eg); // Get element id

            Kij[id] = evaluateMatrix(eg.geometry().center(),L,isTest);
        }

    }

    Dune::FieldMatrix<double,2,2> evaluateMatrix(const Dune::FieldVector<double,2> x,double L,bool test = false)
    {
        double k = 0.0;

        if(isTest){ k = 1.0; }
        else {
        for (int j = 0; j < numKLmodes_; j++){
            //k += std::sqrt(lam[j]) * evalPhi(x[0],mode_id_i[j],L) * evalPhi(x[1],mode_id_j[j],L) * evalPhi(x[2],mode_id_k[j],L) * xi[j];
            k += std::sqrt(lam[j]) * evalPhi(x[0],mode_id_i[j],L) * evalPhi(x[1],mode_id_j[j],L) * xi[j];
        }
        k = std::exp(k);
        }
        Dune::FieldMatrix<double,2,2> Q(0.0);
        for (int i = 0; i < 2; i++){Q[i][i] = k;}

        return Q;
    }

    void inline evaluateTensor(int id, Dune::FieldMatrix<double,2,2>& Q) const{
        Q = Kij[id];
    }

    double inline evaluateScalar(int id) const {
        return Kij[id][0][0];
    }


    void inline user_random_field()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> randn(0.0,1.0);

        std::fill(xi.begin(), xi.end(), 0.0);

        for (int j = 0; j < numKLmodes_; j++){
            xi[j] = sigKL_ * randn(gen);
        }

    } // end user_random_field


private:


    int numKLmodes_;
    bool isTest = false;
    double sigKL_, ellKL_;
    std::vector<double> xi;
    std::vector<double> freq, lam1D, lam, param_;
    std::vector<int> mode_id_i, mode_id_j, mode_id_k;
    std::vector<Dune::FieldMatrix<double,2,2>> Kij;

}; // End Coeff class

template<int dim, typename GRID>
class MODEL{

public:

  // KL parameters used for 2D tests
    double sigKL = 2.5;
    double ellKL = 1./3.;
    int numKLmodes = std::pow(6,dim);

    // KL Parameters used for 3D tests
    //double sigKL = 1.5;
    //double ellKL = 0.5;
  //  int numKLmodes = std::pow(6,dim);

    Dune::FieldVector<double,dim> L;

    int MaxIt = 5000;
    int Verbosity = 1;
    double tolerance = 1e-6;

    int intOrder = 1;

    std::vector<double> mat;

    int dimension = dim;

    GRID& grid_;

    // Constructor for MODEL CLASS
    MODEL(GRID& grid, Dune::FieldVector<double,dim>& L_):grid_(grid), L(L_){
        mat.resize(1); mat[0] = 0.0;

    };

    void inline getAdaptiveSample(int maxLevel,  COEFF& z, std::vector<double>& output) const{

// Initialisation

    // Initialise output

    double Qc = 0.0; double Qf = 0.0;

    double error = 0.0;

    using Dune::PDELab::Backend::native;

    // Load coarse grid
      std::string gridName = "grids/myGrid.msh";
      GRID grid;
      Dune::GridFactory<GRID> factory(&grid);
      Dune::GmshReader<GRID>::read(factory,gridName,false);
      factory.createGrid();


      //  Define Leaf Grid view
      typedef typename GRID::LeafGridView LGV;
      LGV gv = grid.leafGridView(); // Get finest grid


      typedef double RF;
      typedef typename LGV::Grid::ctype Coord;

      // <<<2>>> Make grid function space

      typedef Dune::PDELab::PkLocalFiniteElementMap<LGV,Coord,RF,1> FEM;
      FEM fem(gv);

      typedef Dune::PDELab::ConformingDirichletConstraints CON;    // constraints class
      typedef Dune::PDELab::istl::VectorBackend<> VBE;
      typedef Dune::PDELab::GridFunctionSpace<LGV,FEM,CON,VBE> GFS;
      GFS gfs(gv,fem);  gfs.name("solution");

      // <<<3>>> assemble constraints on this space
      BCTypeParam bctype; // boundary condition type
      typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
      CC cc;
      Dune::PDELab::constraints( bctype, gfs, cc ); // assemble constraints

      // <<<4>>> make DOF vector
      typedef typename Dune::PDELab::Backend::impl::BackendVectorSelector<GFS,double>::Type U;
      

       // initialise pressure & dual solution

      typedef Dirichlet_BC<LGV,double> G; // boundary value + extension
      G g(gv);
    


      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
        MBE mbe(7);


      // Posterior Error Estimator

        typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,RF,dim> P0FEM;
                P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));

        typedef Dune::PDELab::GridFunctionSpace<LGV,P0FEM,Dune::PDELab::NoConstraints,VBE> P0GFS;
        typedef Dune::PDELab::EmptyTransformation NoTrafo;
                  
        using U0 = Dune::PDELab::Backend::Vector<P0GFS,RF>;

        P0GFS p0gfs(gv,p0fem);  p0gfs.name("error");
        U0 eta(p0gfs,0.0), eta2(p0gfs,0.0);



for (int l = 0; l < maxLevel + 1; l++){

    int K = (l == maxLevel) ? 1 : 2;


    for (int k = 0; k < K; k++){

    

        z.computeTensor<LGV>(gv,0.5 * L[0]); // Compute Tensor on new grid

        // (1) Solve primal problem

        U p(gfs,0.0);

        Dune::PDELab::interpolate(g,gfs,p);

        typedef Dune::PDELab::diffuse<COEFF,LGV,3> LOP;
            LOP lop(gv,z,intOrder,false);

        typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,CC,CC> GO;
            GO go(gfs,cc,gfs,cc,lop,mbe);

        // Select a linear solver backend
        typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
            LS ls(5000,0);

        // Select linear problem solver
        typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
            SLP slp(go,ls,p,1e-10,1e-99,0);

        slp.apply(); // Compute solution to test problem

        Dune::VTKWriter<LGV> vtkwriter(gv,Dune::VTK::conforming);
        Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,p);
        vtkwriter.write("solution",Dune::VTK::appendedraw);


        // (2) Solve Dual Problem

        U w(gfs,0.0);

        Dune::PDELab::interpolate(g,gfs,w);

        typedef Dune::PDELab::diffuse<COEFF,LGV,3> DLOP;
            DLOP dlop(gv,z,intOrder,true);

        typedef Dune::PDELab::GridOperator<GFS,GFS,DLOP,MBE,double,double,double,CC,CC> DGO;
            DGO dgo(gfs,cc,gfs,cc,dlop,mbe);

        // Select a linear solver backend
        typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<DGO> DLS;
            DLS dls(5000,0);

        // Select linear problem solver
        typedef Dune::PDELab::StationaryLinearProblemSolver<DGO,DLS,U> DSLP;
            DSLP dslp(dgo,dls,w,1e-10,1e-99,0);

        dslp.apply(); // Compute solution to dual problem

        Dune::VTKWriter<LGV> vtkwriter2(gv,Dune::VTK::conforming);
        Dune::PDELab::addSolutionToVTKWriter(vtkwriter2,gfs,w);
        vtkwriter2.write("influence",Dune::VTK::appendedraw);

        // (3) Calculate QoI if required

        if (k == 0 && l == maxLevel - 1){  Qc = QuantityofInterest<GFS,CC,MBE,U>(gfs,cc,mbe,p);}
        else if (k == 0 && l == maxLevel){  Qf  = QuantityofInterest<GFS,CC,MBE,U>(gfs,cc,mbe,p);}

        // (4) ==== Adapt Grid if required

        if(l < maxLevel){



            // (A) ==== Compute Cheap Error Estimator

            typedef Dune::PDELab::error_estimator_residual_based<LGV,COEFF,false> ERRLOP;
                ERRLOP errLop(gv,z);

            typedef Dune::PDELab::GridOperator<GFS,P0GFS,ERRLOP,MBE,RF,RF,RF,NoTrafo,NoTrafo> ERRGO;
                ERRGO errgo(gfs,p0gfs,errLop,mbe);
                  
            errgo.residual(p,eta);
            errgo.residual(w,eta2);

            for (int ii = 0; ii < native(eta).size(); ii++){    native(eta)[ii] *= native(eta2)[ii];    }



            // (B) ==== Adapt Grid
    
            AdaptMyGrid<GRID,GFS,CC,BCTypeParam,G,U,U0>(grid,gfs,cc,bctype,g,p,eta);

        }

        if (l == maxLevel){
            // For Primal Problem
            std::vector<double> eta_k_p(gv.size(0));

            std::vector<Dune::FieldVector<double,dim>> flux(gv.size(0));

            computeFlux<GFS,LGV,U,COEFF,dim>(gfs,gv,p,z,flux);

            computeImplicitError<GFS,LGV,U,COEFF,dim,3>(gfs,gv,p,z,flux,eta_k_p);

            // For Dual Problem
            std::vector<double> eta_k_w(gv.size(0));

            std::vector<Dune::FieldVector<double,dim>> flux_w(gv.size(0));

            computeFlux<GFS,LGV,U,COEFF,dim>(gfs,gv,w,z,flux_w);

            computeImplicitError<GFS,LGV,U,COEFF,dim,3>(gfs,gv,w,z,flux_w,eta_k_w);

            // Assemble QoI error estimate
            for (int i = 0; i < eta_k_p.size(); i++){
                error += std::sqrt(eta_k_p[i]) * std::sqrt(eta_k_w[i]);
            }

            error = std::sqrt(error);
        }



    }  // For each interlevel adaptive step


} //  For Each Level




output[0] = Qf - Qc;
output[1] = error;



}


void inline adaptit(int maxLevel, COEFF& z, std::vector<double>& output) const{

    using Dune::PDELab::Backend::native;

    double  Qf = 0.0, Qc = 0.0;

        grid_.globalRefine(1);

        typedef typename GRID::LeafGridView LGV;
        LGV gv = grid_.leafGridView(); // Get finest grid

        typedef double RF;
        typedef typename LGV::Grid::ctype Coord;

        const int dofel = (dim == 2) ? 3 : 4;

        z.computeTensor<LGV>(gv,0.5 * L[0],false); // true indicates test case

        typedef Dune::PDELab::PkLocalFiniteElementMap<LGV,Coord,RF,1> FEM;
        FEM fem(gv);

        typedef Dune::PDELab::ConformingDirichletConstraints CON;

        typedef Dune::PDELab::istl::VectorBackend<> VectorBackend;

        typedef Dune::PDELab::GridFunctionSpace<LGV, FEM, CON, VectorBackend> GFS;
        GFS gfs(gv,fem); gfs.name("pressure");

        // <<<3>>> assemble constraints on this space
        BCTypeParam bctype;
        
        typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
            CC cc;
        
        Dune::PDELab::constraints( bctype, gfs, cc );

        typedef typename Dune::PDELab::Backend::impl::BackendVectorSelector<GFS,RF>::Type U;  
            U p(gfs,0.0), w(gfs,0.0);

        for (int i = 0; i <= maxLevel-1; i++){

            int numAdaptSteps = 3;

            for (int k = 0; k < numAdaptSteps;  k++){

                z.computeTensor<LGV>(gv,0.5 * L[0],false);

                //  Make grid operator
                typedef Dune::PDELab::diffuse<COEFF,LGV,dofel> LOP;
                    LOP lop(gv,z,intOrder,false);

                typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
                    MBE mbe(7); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().


                typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,CC,CC> GO;
                    GO go(gfs,cc,gfs,cc,lop,mbe);

                // Select a linear solver backend
                typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
                    LS ls(5000,0);

                // Select linear problem solver
                typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
                    SLP slp(go,ls,p,1e-10,1e-99,0);
                
                slp.apply();

                // Posterior Error Estimator

                typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,RF,dim> P0FEM;
                        P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));

                typedef Dune::PDELab::GridFunctionSpace<LGV,P0FEM,Dune::PDELab::NoConstraints,VectorBackend> P0GFS;
                typedef Dune::PDELab::EmptyTransformation NoTrafo;
                          
                using U0 = Dune::PDELab::Backend::Vector<P0GFS,RF>;

                P0GFS p0gfs(gv,p0fem);  p0gfs.name("error");
                U0 eta(p0gfs,0.0);

                typedef Dune::PDELab::error_estimator_residual_based<LGV,COEFF,false> ERRLOP;
                    ERRLOP errLop(gv,z);

                typedef Dune::PDELab::GridOperator<GFS,P0GFS,ERRLOP,MBE,RF,RF,RF,NoTrafo,NoTrafo> ERRGO;
                    ERRGO errgo(gfs,p0gfs,errLop,mbe);

                errgo.residual(p,eta);

                if (k == 0 && i == maxLevel-1)
                {
                    Qc = QuantityofInterest<GFS,CC,MBE,U>(gfs,cc,mbe,p);
                }



                //

                typedef Dune::PDELab::diffuse<COEFF,LGV,dofel> DLOP;
                    DLOP dlop(gv,z,intOrder,true);

                typedef Dune::PDELab::GridOperator<GFS,GFS,DLOP,MBE,double,double,double,CC,CC> DGO;
                    DGO dgo(gfs,cc,gfs,cc,dlop,mbe);

                // Select a linear solver backend
                typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<DGO> DLS;
                    DLS dls(5000,0);

                // Select linear problem solver
                typedef Dune::PDELab::StationaryLinearProblemSolver<DGO,DLS,U> DSLP;
                    DSLP dslp(dgo,dls,p,1e-10,1e-99,0);

                dslp.apply(); // Compute solution to test problem


                // Compute Error Estimate

                U0 eta2(p0gfs,0.0);
                errgo.residual(p,eta2);

                for (int i = 0; i < native(eta).size(); i++){
                    native(eta)[i] = std::sqrt(native(eta)[i]) * std::sqrt(native(eta2)[i]);
                }

                // ==== Adaptively refine mesh.

                  double refinement_factor = 0.2;
                  double coursening_factor = 0.0;

                  double alpha(refinement_factor);       // refinement fraction
                  double eta_alpha(0);     // refinement threshold
                  double beta(coursening_factor);        // coarsening fraction
                  double eta_beta(0);      // coarsening threshold
                  int verbose = 0;

                  // <<<10>>> Adapt the grid locally...
                  Dune::PDELab::element_fraction( eta, alpha, beta, eta_alpha, eta_beta, verbose );
                  Dune::PDELab::mark_grid( grid_, eta, eta_alpha, eta_beta ,0 , 100, verbose);
                  Dune::PDELab::adapt_grid( grid_, gfs, p, 2 );

                  Dune::PDELab::constraints(bctype,gfs,cc);
              

              }




        }

        z.computeTensor<LGV>(gv,0.5 * L[0],false);

                //  Make grid operator
                typedef Dune::PDELab::diffuse<COEFF,LGV,dofel> LOP;
                    LOP lop(gv,z,intOrder,false);

                typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
                    MBE mbe(7); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().


                typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,CC,CC> GO;
                    GO go(gfs,cc,gfs,cc,lop,mbe);

                // Select a linear solver backend
                typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
                    LS ls(5000,0);

                // Select linear problem solver
                typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
                    SLP slp(go,ls,p,1e-10,1e-99,0);
                
                slp.apply();

                std::vector<double> eta_k_p(gv.size(0));

                std::vector<Dune::FieldVector<double,dim>> flux(gv.size(0));

                computeFlux<GFS,LGV,U,COEFF,dim>(gfs,gv,p,z,flux);

                computeImplicitError<GFS,LGV,U,COEFF,dim,3>(gfs,gv,p,z,flux,eta_k_p);

                Qf = QuantityofInterest<GFS,CC,MBE,U>(gfs,cc,mbe,p);

                typedef Dune::PDELab::diffuse<COEFF,LGV,dofel> DLOP;
                    DLOP dlop(gv,z,intOrder,true);

                typedef Dune::PDELab::GridOperator<GFS,GFS,DLOP,MBE,double,double,double,CC,CC> DGO;
                    DGO dgo(gfs,cc,gfs,cc,dlop,mbe);

                // Select a linear solver backend
                typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<DGO> DLS;
                    DLS dls(5000,0);

                // Select linear problem solver
                typedef Dune::PDELab::StationaryLinearProblemSolver<DGO,DLS,U> DSLP;
                    DSLP dslp(dgo,dls,p,1e-10,1e-99,0);

                dslp.apply(); // Compute solution to test problem

                double error = 0.0;

                std::vector<double> eta_k_w(gv.size(0));


                std::vector<Dune::FieldVector<double,dim>> flux_w(gv.size(0));

                computeFlux<GFS,LGV,U,COEFF,dim>(gfs,gv,p,z,flux_w);

        

        computeImplicitError<GFS,LGV,U,COEFF,dim,3>(gfs,gv,p,z,flux_w,eta_k_w);

        for (int i = 0; i < eta_k_p.size(); i++){
                    error += std::sqrt(eta_k_p[i]) * std::sqrt(eta_k_w[i]);
        }

        error = std::sqrt(error);

        output[0] = Qf  - Qc; output[1] = error;



        



}


    void inline  testErrorEstimator(int l, COEFF& z,std::vector<double>& output, bool computeBias = true) const{

        using Dune::PDELab::Backend::native;

        typedef typename GRID::LevelGridView LGV;
        LGV gv = grid_.levelGridView(l);

        typedef double RF;
        typedef typename LGV::Grid::ctype Coord;

        const int dofel = (dim == 2) ? 3 : 4;

        z.computeTensor<LGV>(gv,0.5 * L[0],false); // true indicates test case

        typedef Dune::PDELab::PkLocalFiniteElementMap<LGV,Coord,RF,1> FEM;
        FEM fem(gv);

        typedef Dune::PDELab::ConformingDirichletConstraints CON;

        typedef Dune::PDELab::istl::VectorBackend<> VectorBackend;

        typedef Dune::PDELab::GridFunctionSpace<LGV, FEM, CON, VectorBackend> GFS;
        GFS gfs(gv,fem); gfs.name("pressure");

        // <<<3>>> assemble constraints on this space
        BCTypeParam bctype;
            typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
            CC cc;
            Dune::PDELab::constraints( bctype, gfs, cc );

        typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
        MBE mbe(18); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().


        //--- Compute primal solution

        typedef Dune::PDELab::diffuse<COEFF,LGV,dofel> LOP;
        LOP lop(gv,z,intOrder,false);

        typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,CC,CC> GO;
        GO go(gfs,cc,gfs,cc,lop,mbe);

        // Make coefficent vector and initialize it from a function
        typedef typename GO::Traits::Domain Up;
        Up p(gfs,0.0);

        // Select a linear solver backend
        typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
        LS ls(5000,0);

        // Select linear problem solver
        typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,Up> SLP;
        SLP slp(go,ls,p,1e-10,1e-99,0);

        slp.apply(); // Compute solution to test problem

        Dune::VTKWriter<LGV> vtkwriter(gv,Dune::VTK::conforming);
        Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,p);
        vtkwriter.write("solution",Dune::VTK::appendedraw);


        std::vector<double> eta_k_p(gv.size(0));


        if(computeBias){
            std::vector<Dune::FieldVector<double,dim>> flux(gv.size(0));

            computeFlux<GFS,LGV,Up,COEFF,dim>(gfs,gv,p,z,flux);

            computeImplicitError<GFS,LGV,Up,COEFF,dim,3>(gfs,gv,p,z,flux,eta_k_p);
        }

        
        // Compute dual problem

        //--- Get primal solution

        typedef Dune::PDELab::diffuse<COEFF,LGV,dofel> DLOP;
        DLOP dlop(gv,z,intOrder,true);

        typedef Dune::PDELab::GridOperator<GFS,GFS,DLOP,MBE,double,double,double,CC,CC> DGO;
        DGO dgo(gfs,cc,gfs,cc,dlop,mbe);

        // Make coefficent vector and initialize it from a function
        typedef typename DGO::Traits::Domain Uw;
        Uw w(gfs,0.0);

        // Select a linear solver backend
        typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<DGO> DLS;
        DLS dls(5000,0);

        // Select linear problem solver
        typedef Dune::PDELab::StationaryLinearProblemSolver<DGO,DLS,Uw> DSLP;
        DSLP dslp(dgo,dls,w,1e-10,1e-99,0);

        dslp.apply(); // Compute solution to test problem

        Dune::VTKWriter<LGV> vtkwriter2(gv,Dune::VTK::conforming);
        Dune::PDELab::addSolutionToVTKWriter(vtkwriter2,gfs,w);
        vtkwriter2.write("two",Dune::VTK::appendedraw);


        double error = 0.0;

        std::vector<double> eta_k_w(gv.size(0));

        if(computeBias){

        std::vector<Dune::FieldVector<double,dim>> flux_w(gv.size(0));

        computeFlux<GFS,LGV,Uw,COEFF,dim>(gfs,gv,w,z,flux_w);

        

        computeImplicitError<GFS,LGV,Uw,COEFF,dim,3>(gfs,gv,w,z,flux_w,eta_k_w);

        for (int i = 0; i < eta_k_p.size(); i++){
                    error += std::sqrt(eta_k_p[i]) * std::sqrt(eta_k_w[i]);
        }

        error = std::sqrt(error);

        }

       

        // Find QoI

        double Q = QuantityofInterest<GFS,CC,MBE,Up>(gfs,cc,mbe,p);

        // Compute  Bias Error


        output[0] = Q; output[1] = error;


    }





private:




};











#endif /* user_model_h */
