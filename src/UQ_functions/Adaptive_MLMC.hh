//
//  STD_MLMC.hh
//
//
//  Created by Tim Dodwell on 30/11/2015.
//
//

#ifndef STD_MLMC_h
#define STD_MLMC_h

#include "../user_inputs/user_model_3Delasticity.hh"
#include "parallel.hh"

using namespace std;
//#include "../user_inputs/linearelasticity3D.hh"

template <class MODEL>
class STD_MLMC{
    
public:
    
    STD_MLMC (MODEL& model):model_(model)
    {
    }
    
    
    void inline setParameter(double epsilon_ = 1e-2,int verbosity_ = 1, double initialSamples_Coarsest_ = 100){
        epsilon = epsilon_;
        verbosity = verbosity_;
        initialSamples_Coarsest = initialSamples_Coarsest_;
    }
    
    void inline apply(){
        
        int rank, nproc;
        
        MPI_Comm new_comm;
        
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        if (rank == 0){
            std::cout << "----------------------------------" << std::endl;
            std::cout << "-    dune-uq : Adaptive MLMC     -" << std::endl;
            std::cout << "----------------------------------" << std::endl;
        }
        
        COEFF z(model_.L[0],model_.numKLmodes,model_.sigKL,model_.ellKL,model_.mat);
        
        
        int L = -1;
        double BiasError = epsilon;
        
        Dune::Timer timer;
        
        while (BiasError > epsilon/(std::sqrt(2.0))){
            
            L += 1; // Increment Number of Levels
            
            Nopt.resize(L + 1);
            N.resize(L + 1);
            
            if (L > 0)model_.grid_.globalRefine(1);
            
            if (rank == 0){
                
                if (L > 0 && verbosity > 0) { std::cout << " Add extra level! ";} // Refine grid once (Uniformly)
                
                std::cout << "Current level : " << L << std::endl;
                
                // << 1 >> Compute Initial samples on level l
                
                // < a > Resize and Initialise  Cost, Nopt, N and Y Vectors
                Cost.resize(L + 1);
                Y.resize(L + 1);
                
                
            }
            
            if (L == 0){    Nopt[L] = initialSamples_Coarsest;}
            else{   Nopt[L] = std::max(N[L-1]/100,5);}
            
            N[L] = Nopt[L];
            
            // < b > Compute Initital Samples in Parallel
            
            int numSamples = N[L] / nproc + 1; // Compute number of samples on each processor (rounding up)
            
            if (rank == 0) {
                std::cout << " Computing " << N[L] <<" Inital Samples" << std::endl;
                timer.reset(); // Start timer on processor 0
            }
            
            double * Ytmp = new double[numSamples];
            
            double * Yroot = NULL;
            
            if (rank == 0){ Yroot = new double[numSamples * nproc]; }
            
            for (int i = 0; i < numSamples; i++){ Ytmp[i] = getSample(L,z); } // Compute Samples on each processor
            
            // Gather all samples from each processor to root processor
            MPI_Gather(Ytmp,numSamples,MPI_DOUBLE,Yroot,numSamples,MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            delete Ytmp;
            
            N[L] = numSamples * nproc;  Nopt[L] = numSamples * nproc;
            
            double k; // Initialise Constant of Proportionality
            
            std::vector<double> V; // Vector to store Variances on each level
            
            
            if (rank == 0){
                
                Y[L].resize(Nopt[L]);
                
                for (int i = 0; i < Nopt[L]; i++){ Y[L][i] = Yroot[i]; }
                
                Cost[L] = timer.elapsed() / Nopt[L]; // Record Average (Parallel) Cost (in secs)
                
                std::cout << " Average Cost =  " << Cost[L] << std::endl;
                
                
                // << 2 >> Compute variance, V, and constant of proportionality for N_l
                
                double tmp = 0.0;
                
                V.resize(L + 1);
                for (int l = 0; l < L + 1; l++){
                    V[l] = var(Y[l]);
                    tmp += std::sqrt(V[l] * Cost[l]);
                }
                
                k = 2.0 * tmp / (epsilon * epsilon);
                
            }
            
            delete Yroot;
            
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            
            // << 3 >> For each level compute extra samples if required
            
            for (int l = 0; l < L + 1; l++)
            {
                
                
                if (rank == 0){
                    // < a > Compute optimal value number of samples on level l
                    Nopt[l] = std::ceil(k * std::sqrt(V[l]/Cost[l]));
                }
                
                MPI_Bcast(&Nopt[l],1,MPI_INT,0,MPI_COMM_WORLD);
                
                
                
                // < b > If Nopt[l] > N[l], then additional samples are evaluated at each
                if (Nopt[l] > N[l]){
                    
                    int moreSamples = Nopt[l] - N[l];
                    
                    if (rank == 0){
                        std::cout << " Computing " << moreSamples <<" more samples on level " << l << std::endl;
                        
                    }
                    
                    int moreSamples_per_processor = moreSamples / nproc + 1;
                    
                    Nopt[l] = N[l] + (moreSamples_per_processor * nproc);
                    
                    double * Ytmp = new double[moreSamples_per_processor];
                    
                    double * Yroot = NULL;
                    
                    if (rank == 0){ Yroot = new double[moreSamples_per_processor * nproc]; }
                    
                    for (int i = 0; i < moreSamples_per_processor; i++){ Ytmp[i] = getSample(l,z); } // Compute Samples on each processor
                    
                    // Gather all samples from each processor to root processor
                    MPI_Gather(Ytmp,moreSamples_per_processor,MPI_DOUBLE,Yroot,moreSamples_per_processor,MPI_DOUBLE,0,MPI_COMM_WORLD);
                    
                    delete Ytmp;
                    
                    if (rank == 0){
                        
                        Y[l].resize(Nopt[l]);
                        
                        for (int i = N[l]; i < Nopt[l]; i++){ Y[l][i] = Yroot[i-N[l]]; }
                        
                    } // end rank
                    
                    delete Yroot;
                    
                    
                    N[l] = Nopt[l];
                    
                    
                } // End if more samples are required
                
                
            } // End for each level
            
            // << 4 >> Compute Bias Error
            
            
            double FEMerror = 0.0;
            
            if (rank == 0){
                
                
                FEMerror = 0.0;
                for (int i = 0; i < Nopt[L]; i++){
                    FEMerror += std::abs(Y[L][i]);
                }
                FEMerror /= N[L];
                
                
                std::cout << " Bias error on level " << L << " = " << FEMerror << std::endl;
                
                Bias.resize(L+1);
                Bias[L] = FEMerror;
            }
            
            
            MPI_Bcast(&FEMerror,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            BiasError = FEMerror;
            
        }
        
        
        if (rank == 0){
            std::cout << "|| L | " << "N | " << "Cost |" << "BiasError ||" << std::endl;
            
            for (int i = 0; i < L + 1; i++){
                std::cout << "|| " << i << " | " << N[i] << " | "<< Cost[i] << " | "<< Bias[i] << " ||"<< std::endl;
            }
        }
        
        
    }
    
    
private:
    
    
    template <class Z>
    double getSample(const int l,Z& z)
    {
        
        // Generate Random Sample
        z.user_random_field();
        
        // Compute Quantities of Interest using FEM
        double Qc = 0.0, Qf = 0.0;
        
        Qf = model_.user_fem_driver(l,z);
        // If l > 0 we compute coarse grid QoI (Qc), else defaults to Qc = 0.0;
        if (l > 0){
            Qc = model_.user_fem_driver(l-1,z);
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
    
};






#endif /* STD_MLMC_h */
