//
//  Cost_Test.hh
//
//
//  Created by Tim Dodwell on 30/11/2015.
//
//

#ifndef STD_MLMC_h
#define STD_MLMC_h


#include "parallel.hh"

using namespace std;


template <class MODEL>
class Cost_Test{
    
public:
    
    STD_MLMC (MODEL& model):model_(model)
    {
    }
    
    
    void inline setParameter(double epsilon_ = 1e-2,int verbosity_ = 1, double initialSamples_Coarsest_ = 100, int adaptive_ = 1){
        epsilon = epsilon_;
        verbosity = verbosity_;
        initialSamples_Coarsest = initialSamples_Coarsest_;
        adaptive = adaptive_;
    }
    
    void inline apply(){
        
        int rank, nproc;
        
        MPI_Comm new_comm;
        
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        if (rank == 0){
            std::cout << "----------------------------------" << std::endl;
            if (adaptive == 1) {
                std::cout << "-    dune-uq : Adaptive MLMC     -" << std::endl;
            }
            else{
                std::cout << "-    dune-uq : Standard MLMC     -" << std::endl;
            }
            std::cout << "----------------------------------" << std::endl;
        }
        
        N = 30;
        L = 5;
        
        std::vector<double> C(L+1), M(L+1);

        
        int numSamples = N / nproc + 1;
        
        N = numSamples * nproc;
        
        COEFF z(model_.L[0],model_.numKLmodes,model_.sigKL,model_.ellKL,model_.mat);
        
        Dune::Timer timer;
        
        for (int l = 0; l < L)
            
            timer.reset(); // Reset timer
            
            if (L > 0 && adaptive == 0){
                model_.grid_.globalRefine(1);
            }
        
            // Compute Samples in Parallel
            
             // Compute number of samples on each processor (rounding up)
        
            
            double * Ytmp = new double[numSamples];
            
            double * Yroot = NULL;
            
            if (rank == 0){ Yroot = new double[numSamples * nproc]; }
            
            double sum = 0.0;
            
            for (int i = 0; i < numSamples; i++){ Ytmp[i] = getSample(L,z,adaptive);
            } // Compute Samples on each processor
        
            MPI_Barrier(MPI_COMM_WORLD);
            
            // Gather all samples from each processor to root processor
            MPI_Gather(Ytmp,numSamples,MPI_DOUBLE,Yroot,numSamples,MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            delete Ytmp;
        
            if (rank == 0){ C[l] = timer.elapsed() / N;
        
        }
    
        
        
        
        
        
        if (rank == 0){
            std::cout << "|| L | " << "M | " << "Cost |" << std::endl;
            
            for (int i = 0; i < L + 1; i++){
                std::cout << "|| " << i << " | " << M[i] << " | "<< C[i] << "|" << std::endl;
            }
            
        }
        
        
    }
    
    
private:
    
    
    template <class Z>
    double getSample(const int l, Z& z, int mode = 0)
    {
        
        // Generate Random Sample
        z.user_random_field();
        
        // Compute Quantities of Interest using FEM
        double Qc = 0.0, Qf = 0.0;
        
        switch (mode) {
                
            case 0: // standard MLMC sampler
                
                Qf = model_.user_fem_driver(l,z);
                
                // If l > 0 we compute coarse grid QoI (Qc), else defaults to Qc = 0.0;
                if (l > 0){ Qc = model_.user_fem_driver(l-1,z);}
                
                break;
                
            case 1: // adaptive MLMC sampler
                
                Qf = model_.user_adaptive_fem_driver(l,z);
                
                break;
                
            default:
                cout << "value of mode unknown - mode = 0 Standard MLMC mode = 1 adaptive MLMC" << endl;
        }
        
        
        
        double Y = Qf - Qc;
        
        return Y;
        
    }
    
    const MODEL model_;
    std::vector<double> Cost;
    std::vector<int> Nopt, N;
    std::vector<double> Bias;
    std::vector<std::vector<double>> Y;
    
    double epsilon = 1e-2;
    int verbosity = 1;
    int initialSamples_Coarsest = 100;
    int adaptive;
    
};






#endif /* STD_MLMC_h */
