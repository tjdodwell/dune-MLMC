


std::vector<double> inline getAdaptiveSample(int L,  COEFF& z) const{

// Initialisation

	// Initialise output

	double Qc = 0.0; double Qf = 0.0;

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
      typedef Dune::PDELab::ISTLVectorBackend<> VBE;
      typedef Dune::PDELab::GridFunctionSpace<LGV,FEM,CON,VBE> GFS;
      GFS gfs(gv,fem);  gfs.name("solution");

      // <<<3>>> assemble constraints on this space
      BCTypeParam bctype; // boundary condition type
      typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
      CC cc;
      Dune::PDELab::constraints( bctype, gfs, cc ); // assemble constraints

      // <<<4>>> make DOF vector
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,double>::Type U;
      U p(gfs,0.0);
      U w(gfs,0.0); // initialise pressure & dual solution

      typedef Dirichlet_BC<LGV,double> G; // boundary value + extension
      G g(gv);
      Dune::PDELab::interpolate(g,gfs,p);


      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
        MBE mbe(7);


      // Posterior Error Estimator

        typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,RF,dim> P0FEM;
            	P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));

            typedef Dune::PDELab::GridFunctionSpace<LGV,P0FEM,Dune::PDELab::NoConstraints,VBE> P0GFS;
            typedef Dune::PDELab::EmptyTransformation NoTrafo;
                  
                  using U0 = Dune::PDELab::Backend::Vector<P0GFS,RF>;

                  P0GFS p0gfs(gv,p0fem);
                    p0gfs.name("error");
                  U0 eta(p0gfs,0.0), eta2(p0gfs,0.0);



for (int l = 0; l < L + 1; l++){

	K = (l == L) ? 1 : numadaptiveSteps;


	for (int k = 0; k < K; k++){

		z.computeTensor<LGV>(gv,0.5 * L[0]); // Compute Tensor on new grid

		// (1) Solve primal problem

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


		// (2) Solve Dual Problem

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

		// (3) Calculate QoI if required

		if (k == 0 && l == L - 1){	Qc = QuantityofInterest<GFS,CC,MBE,Up>(gfs,cc,mbe,p);}
		elseif (k == 0 && l == L){ 	Qf  = QuantityofInterest<GFS,CC,MBE,Up>(gfs,cc,mbe,p);}


		// (4) ====	Adapt Grid if required

		if(l < maxLevel){

			// (A) ====	Compute Cheap Error Estimator

            typedef Dune::PDELab::error_estimator_residual_based<LGV,COEFF,false> ERRLOP;
            	ERRLOP errLop(gv,z);

            typedef Dune::PDELab::GridOperator<GFS,P0GFS,ERRLOP,MBE,RF,RF,RF,NoTrafo,NoTrafo> ERRGO;
            	ERRGO errgo(gfs,p0gfs,errLop,mbe);
                  
            errgo.residual(p,eta);
            errgo.residual(w,eta2);

            for (int ii = 0; ii < native(eta).size(); ii++){	native(eta)[ii] *= native(eta2)[ii];	}

			// (B) ====	Adapt Grid

           	AdaptMyGrid<GRID,GFS,CC,BCTypeParam,G,U,U0>(grid, gfs, cc, bctype, g, p, eta);

		}


	}  // For each interlevel adaptive step


} //  For Each Level

// (5)	Compute Bias Error on Level L


// For Primal Problem
std::vector<double> eta_k_p(gv.size(0));

std::vector<Dune::FieldVector<double,dim>> flux(gv.size(0));

computeFlux<GFS,LGV,Up,COEFF,dim>(gfs,gv,p,z,flux);

computeImplicitError<GFS,LGV,Up,COEFF,dim,3>(gfs,gv,p,z,flux,eta_k_p);

// For Dual Problem
std::vector<double> eta_k_w(gv.size(0));

std::vector<Dune::FieldVector<double,dim>> flux_w(gv.size(0));

computeFlux<GFS,LGV,Uw,COEFF,dim>(gfs,gv,w,z,flux_w);

computeImplicitError<GFS,LGV,Uw,COEFF,dim,3>(gfs,gv,w,z,flux_w,eta_k_w);

// Assemble QoI error estimate
double error = 0.0;
for (int i = 0; i < eta_k_p.size(); i++){
	error += std::sqrt(eta_k_p[i]) * std::sqrt(eta_k_w[i]);
}

error = std::sqrt(error);


std::vector<double> output(2);

output[0] = Qf - Qc;
output[1] = error;

return output

}



