
#include "../PDEs/Cosserat.hh"
#include "user_QuantityofInterest.hh"

// ------------------------------------------------------------------------
//             Dirichlet Boundary Conditions Class
// ------------------------------------------------------------------------
template<typename GV, typename RF>
class Dirichlet_BC_Test_2 :
public Dune::PDELab::AnalyticGridFunctionBase<
Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
Dirichlet_BC<GV,RF> >,
public Dune::PDELab::InstationaryFunctionDefaults
{
private:
  bool isTest;
  int dof;
  double AppliedStrain = 0.1;
  double L;
public:

    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, Dirichlet_BC<GV,RF>> BaseT;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    // Constructor
    Dirichlet_BC_Test_2(const GV & gv,double L_,bool Test = false) : BaseT(gv), L(L_),isTest(Test){}

    template<typename I>
    bool isDirichlet(const I & ig,
                     const typename Dune::FieldVector<typename I::ctype, I::dimension-1> & x
                     ) const
    {
        Dune::FieldVector<double,2> xg = ig.geometry().global(x);

        bool answer = false;

      // if (dof == 3){ answer = true;}
        if ((dof == 1) && (xg[0] < 1e-10 || xg[0]>L-1e-10)){
          answer =  true;
        }
        if ((dof == 2) && (xg[1] < 1e-10 || xg[1]>L-1e-10)){
          answer = true;
        }

        return answer;
    }

    inline void setDof(int d){ dof = d;}

    inline void evaluateGlobal(const DomainType & x, RangeType & u) const
    {
      u = 0.0;

      if  (dof ==  1  && x[0] > L - 1e-10){ u = -AppliedStrain * L;}

    } // end inline function evaluateGlobal


}; // End of Dirichlet BC class




template<int dim, class GRID>
class TestCase2{

public:

  int numKLmodes;

  Dune::FieldVector<double,dim> L;

  double Vf, d, E1, E2, G12, nu12,ellx,elly, sigKL,tau_y, PatchSize,CB, maxLevel;

  Dune::FieldMatrix<double,6,6> D;

  // Solver Paremeters

  int MaxIt = 5000;
  int Verbosity = 1;
  double tolerance = 1e-4;

  GRID& grid;

  TestCase2(GRID& grid_, Dune::FieldVector<double,dim> & L_):grid(grid_), L(L_){

    Dune::FieldMatrix<double,3,3> Q(0.0);

    Vf = config.get<double>("Material.Vf",0.5); // Vf     - Fibre Volume Fraction
    d = config.get<double>("Material.d",0.0071); // mm  - Diameter of a  Fibre (mm)
    sigKL = config.get<double>("RandomField.sigKL",0.035); // rad  - Standard Deviation of Fibre Misalignment
    ellx = config.get<double>("RandomField.ellx",1.6); // mm    - Correlation Length Scale Parallel to the Fibres
    elly = config.get<double>("RandomField.elly",0.43); // mm    - Correlation Length Scale Parallel to the Fibres

    E1 = config.get<double>("Material.E1",128.0); // GPa - Longitudial Stiffness
    E2 = config.get<double>("Material.E2",9.25); // GPa   - Tranverse Stiffness
    G12 = config.get<double>("Material.G12",5.1); // GPa    - Inplane shear stiffness
    nu12 = config.get<double>("Material.nu12",0.35);//        - Poisson Ratio

    tau_y = config.get<double>("Material.tau_y",0.114);// GPa   - Inplane Shear Strength

    PatchSize = config.get<double>("Domain.PatchSize",0.5); //       - Relative Patch Size to Gauge Length

    numKLmodes = config.get<int>("RandomField.numKLmodes",100);; //         - Number of KL Modes

    Qstar = config.get<double>("QoI.Qstar",14.3);

  
    computeM();

  };


  void inline computeM(){

    int maxLevel = config.get<int>("MLMC.maxLevel",5);

    int nel0 = config.get<int>("MLMC.CoarseMeshElements",8);
    
    M_1D.resize(maxLevel);  M.resize(maxLevel);

    M_1D[0] = (double)nel0 + 1.0;
    M[0]    = M_1D[0] * M_1D[0];
    for (int i = 1; i < maxLevel + 1; i++){
      M_1D[i] = 2.0 * M_1D[i-1] - 1.0;
      M[i] = M_1D[i] * M_1D[i];
    }

  }

  double inline getM(int level){ return M[level]; }


  double inline getSample(int l, RandomField& z){

      double alpha_rate = config.get<double>("MLMC.alpha",1.0);

      z.user_random_field(false); // Draw Random Sample from Prior

      double QoI = 0.0;

      double Q0 = 0.0; double Q1 = 10.0 * Qstar;

      if(config.get<bool>("MLMC.selectiveRefinement",false) == false || l < 2){

         Q0 = Solve(l,z);

        if (l > 0){QoI -= Solve(l - 1,z);}

      }
      else{ // level l > 1

        int i = 1;

         Q0 = Solve(0,z);
         Q1 = Solve(1,z);

        double check = std::abs(Q1 - Qstar) - std::abs(Q1 - Q0) / (std::pow(4.0,alpha_rate) - 1.0);

        while (i < l & check < 0.0){
          i += 1;
          Q0 = Q1;
          Q1 = Solve(i,z);
        }

      }

      if (Q0 < Qstar){QoI = 1.0;}

      if (Q1 < Qstar){QoI -= 1.0;}

      return QoI;

  }

  double inline Solve(int l, RandomField& z) const{


    const int numDof = 12;

    // === Define Level Grid view
    typedef typename GRID::LevelGridView LGV;
      LGV gv = grid.levelGridView(l);

      z.setLevel(l);

    typedef typename LGV::Grid::ctype Coord;

    // === Build finite element spaces - Piecewise Constant, Linear and Quadratic


    typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord,double,0,2> FEM0;
      FEM0 fem0;

    typedef Dune::PDELab::QkLocalFiniteElementMap<LGV,Coord,double,1> FEM1;
      FEM1 fem1(gv);
    typedef Dune::PDELab::QkLocalFiniteElementMap<LGV,Coord,double,2> FEM2;
      FEM2 fem2(gv);

    typedef Dune::PDELab::ConformingDirichletConstraints CON;
      CON con;
    typedef Dune::PDELab::istl::VectorBackend<> VectorBackend;

    typedef Dune::PDELab::GridFunctionSpace<LGV, FEM2, CON, VectorBackend> P2_GFS;
    //    P2_GFS dispU1(gv,fem2,con); dispU1.name("U1");
    //    P2_GFS dispU2(gv,fem2,con); dispU2.name("U2");

    typedef Dune::PDELab::GridFunctionSpace<LGV, FEM1, CON, VectorBackend> P1_GFS;
        P1_GFS rot3(gv,fem1,con); rot3.name("ROT3");
        P1_GFS dispU1(gv,fem1,con); dispU1.name("U1");
        P1_GFS dispU2(gv,fem1,con); dispU2.name("U2");

  //  typedef Dune::PDELab::CompositeGridFunctionSpace <VectorBackend,Dune::PDELab::LexicographicOrderingTag, P2_GFS, P2_GFS, P1_GFS> GFS;
    typedef Dune::PDELab::CompositeGridFunctionSpace <VectorBackend,Dune::PDELab::LexicographicOrderingTag, P1_GFS, P1_GFS, P1_GFS> GFS;
          GFS gfs(dispU1, dispU2, rot3);

    // Make constraints map and initialize it from a function
    typedef typename GFS::template ConstraintsContainer<double>::Type C;
    C cg;
    cg.clear();

    // Apply Dirichlet Constraints

     typedef Dirichlet_BC_Test_2<LGV,double> BC;
     BC U1_cc(gv,L[0]), U2_cc(gv,L[0]), ROT3_cc(gv,L[0]);
     U1_cc.setDof(1);
     U2_cc.setDof(2);
     ROT3_cc.setDof(3);

     typedef Dirichlet_BC_Test_2<LGV,double> InitialDisp;
     InitialDisp u1(gv,L[0]), u2(gv,L[0]), rrot3(gv,L[0]);
     u1.setDof(1);
     u2.setDof(2);
     rrot3.setDof(3);

     // Wrap scalar boundary conditions in to vector
     typedef Dune::PDELab::CompositeGridFunction<InitialDisp,InitialDisp,InitialDisp>  InitialSolution;
     InitialSolution initial_solution(u1,u2,rrot3);

     typedef Dune::PDELab::CompositeConstraintsParameters<BC,BC,BC>
        Constraints;
     Constraints constraints(U1_cc,U2_cc,ROT3_cc);

    Dune::PDELab::constraints(constraints,gfs,cg);

    //  Construct Linear Operator on FEM Space
    typedef Dune::PDELab::Cosserat<RandomField,numDof> LOP;
       LOP lop(z);

    typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
      MBE mbe(27); // For all linear shape functions
    //  MBE mbe(59); // Maximal number of nonzeroes per row (25 * 2 displacements + 9 rotations)
    typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,C,C> GO;
      GO go(gfs,cg,gfs,cg,lop,mbe);

    typedef typename GO::Traits::Domain V;
      V x(gfs);
    x = 0.0;
    Dune::PDELab::interpolate(initial_solution,gfs,x);

    Dune::VTKWriter<LGV> vtkwriter2(gv,Dune::VTK::conforming);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter2,gfs,x);
    vtkwriter2.write("initial_solution",Dune::VTK::appendedraw);

    typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> SEQ_CG_AMG_SSOR;
      SEQ_CG_AMG_SSOR ls(MaxIt,0);

    Dune::PDELab::StationaryLinearProblemSolver<GO,SEQ_CG_AMG_SSOR,V> slp(go,ls,x,tolerance,1e-99,0);
      slp.apply();

  /*  Dune::VTKWriter<LGV> vtkwriter(gv,Dune::VTK::conforming);
      Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
      vtkwriter.write("solution",Dune::VTK::appendedraw);*/

    // Need to compute stress over right hand boundary

    using Dune::PDELab::Backend::native;

    //std::cout << "Degrees of Freedom = " << native(x).size() << std::endl;

    typedef Dune::PDELab::GridFunctionSpace<LGV, FEM0, CON, VectorBackend> P0_GFS;
          P0_GFS gfs0(gv,fem0,con); gfs0.name("FC");

    double Q = QuantityofInterest_TestCase2<RandomField,P0_GFS,GFS,MBE,V,LGV,numDof>(z,gfs0,gfs,mbe,x,gv,PatchSize * L[0],E2/G12,tau_y);

    return Q;
  }

  void inline refineGrid(int l = 1){ grid.globalRefine(l);}

private:

  std::vector<double> M_1D, M;

  double Qstar;



};
