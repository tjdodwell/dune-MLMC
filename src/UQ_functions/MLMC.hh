#include <random>
#include "general.hh"

template <class G>
void MLMC(G& grid)
{
    typedef typename G :: LevelGridView LGV;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> randu(0.0, 1.0);
    
    
    std::vector<int> Nopt(grid.maxLevel());
    std::vector<double> C(grid.maxLevel());
    
    // On level 0 carry out some initial samples
    
    // Extract Grid on level 0
    
    Nopt[0] = 10;
    LGV levelView = grid.levelGridView(0);
    
    std::vector<double> cost(Nopt[0]);
    
    for (int i = 0; i < Nopt[0]; i++){
        
        cost[i] = randu(gen);
        
        C[0] += cost[i];
        
    }
    
    C[0] /= Nopt[0];
    
    
    double error = 1.0;
    
    double FEMerror = error;
    
    while (FEMerror > error / std::sqrt(2.0))
    {
        
        // Refine Grid by one level
        grid.globalrefine(1);
        
        double L = grid.maxLevel();
        
        // Resize Vectors
        
        
        // If not the coarsest level
        if (L > 0){
            
            N[L] = 30;

            // Do some initial samples on level L
            for (int i = 0; i < N[L]; i++){
                
                Y[L][i] = getSample<G,M>(grid,model,l);
            }
            
            C[L] = time / N[L];
            
        }
        
        // Compute Variances
        std::vector<double> V(L);
        double sumVlCl = 0.0;
        
        for (int l = 0; l <= L; l++){
            V[l] = var();
            sumVlCl = std::sqrt(V[l] * C[l]);
        }
        
        double k = 2.0 * sumVlCl / (error * error);
        
        for (int l = 0; l <= L; l++){
            
            Nopt[l] = std::ceil(k * std::sqrt(V[l] / C[l]));
            
            // If extra samples on level l are required do them
            if (Nopt[l] < N[l]){
                
                int moreSamples = Nopt[l] - N[l];
                
                for (int i = N[l] - 1; i < Nopt[l]; i++){
                    
                    Y[l][i] = getSample<G,M>(grid,model,l);
    
                }
             
                N[l] = Nopt[l]; // Update number of samples on level l
            }
            
        
        }
        
        FEMerror = std::abs(mean(Y[l]));
        
        std::cout << "FEMerror " << FEMerror << std::endl;
        
        
    }
    
}

