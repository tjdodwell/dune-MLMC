#include "parallel_Greedy_comms.hh"

template <class M,class Samples>
class GreedyMLMC{

public:

    GreedyMLMC (M& model_):model(model_){ }

    void inline apply(double eps){


      int rank, nproc;

      MPI_Comm new_comm;

      MPI_Comm_size(MPI_COMM_WORLD, &nproc);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      RandomField z;

      ParallelStats myStats;

      Samples Q;

      Dune::Timer watch, overallwatch; // Constructor a timer for estimating cost

      double alpha_rate = config.get<double>("MLMC.alpha",1.0);

      int maxLevel = config.get<int>("MLMC.maxLevel",5);


      //  **** Start of Greedy Algorithm for MLMC ****

      overallwatch.reset(); // Reset Timer for Overall Calculation


      // *** Number of samples on each level

      std::vector<int> Nstar = config.get<std::vector<int>>("MLMC.Nstar",{16,16,16,16,16,16,16,16,16,16});

      for (int i = 0; i < Nstar.size(); i++){
        if(Nstar[i] < 3 * nproc){ Nstar[i] = 3 * nproc; }
      }


      //

      double greedyFactor = config.get<double>("MLMC.greedyFactor",0.1);

      double theta = config.get<double>("MLMC.theta", 0.5);

      double Bias = 10.0 * eps;
      double e2 = eps * eps;

      std::vector<double> PayOff, localEQ, localEQ2;

      int L = -1;

      bool converged = false;


      // Start convergence loop

      while (converged == false && L < maxLevel){

        L += 1; // Increment Level counter

        Q.addLevel(); // Add new level to each local sample container

        if(rank == 0 && config.get<bool>("print.inSim",true)){ 
          std::cout << "Adding Level " << L << " and computing " << Nstar[L] << " initial Samples across " << nproc << " processors." << std::endl; 
        }

        PayOff.resize(L+1);
        localEQ.resize(L + 1); 
        localEQ2.resize(L+1);

        // ** Compute initial samples across all available processors

        for (int i = 0; i < std::ceil(Nstar[L] / nproc); i++){
          watch.reset();
          double QoI = model.getSample(L,z);
          Q.addSample(QoI,L,watch.elapsed());
        }

        Q.UpdateStats(); // Update Local Statistics, now have new samples

        MPI_Barrier(MPI_COMM_WORLD); // Wait for all processors to complete samples

        myStats = parallel_Stats_Update<Samples>(Q,rank,L,nproc,alpha_rate);

        while ( myStats.VoE > (std::sqrt(theta) * e2) ){

          // Compute new greedy samples

          int newSamples = std::ceil(greedyFactor * Q.getN(myStats.greedyLevel) / nproc);
         
          for (int i = 0; i < newSamples; i++){
             watch.reset();
             double QoI = model.getSample(myStats.greedyLevel,z);
             Q.addSample(QoI,myStats.greedyLevel,watch.elapsed());
          }

          Q.UpdateStats(); // Update Local Statistics

          MPI_Barrier(MPI_COMM_WORLD); // Wait for all processors to complete samples

          myStats = parallel_Stats_Update<Samples>(Q,rank,L,nproc,alpha_rate);

        } // end while loop


        if(rank == 0 && config.get<bool>("print.inSim",true)){
            std::cout << "Current Level = " << L << " Bias = " << myStats.Bias << " Variance of Estimator = " << myStats.VoE << std::endl;
        }

        if (myStats.Bias * myStats.Bias < (1.0 - theta) * e2){
            converged = true;
        }


      } // end while

      double totalTime = overallwatch.elapsed();

      // Post Process Results - Write to file

      // Gather all results.



     if (rank == 0 && config.get<bool>("print.Stats", true)){

        std::cout << "****** MLMC Greedy Algorithm Completed ******" << std::endl;
        if (L == maxLevel && converged == false){
          std::cout << "!!! Run not converged before maxLevel reached !!!" << std::endl;
        }
        std::cout << "Total Time Elapsed = " << totalTime << "secs" << std::endl;
        std::cout << "Total Levels = " << L + 1 << std::endl;
        std::cout << "Qhat = " << myStats.expectedValue << std::endl;
        std::cout << "Var(Qhat) = " << myStats.VoE << std::endl;
        double Tol = std::sqrt(myStats.VoE + myStats.Bias * myStats.Bias);
        std::cout << "Tolerance  = " << Tol << " ( As Percentage of target = " << (Tol/eps) << "%)" << std::endl;
        std::cout << "Bias = " << myStats.Bias << std::endl;
        std::cout << "Sample Error = " << std::sqrt(myStats.VoE) << std::endl;

        // Estimate the Rates for this simulation
        std::vector<double> MM(L), EY(L), EC(L), V(L);
        for (int i = 1; i < L+1; i++){
          MM[i-1] = model.getM(i);
          EY[i-1] = std::abs(myStats.getEQ(i));
          EC[i-1] = Q.ExpectedCost(i);
          V[i-1] = myStats.getVQ(i);

        }
        alpha_rate = findRate(MM,EY);
        double beta = findRate(MM,V);
        double gamma = findRate(MM,EC);

        //std::cout << " Final theta Value = " << theta << std::endl;
        std::cout << "alpha = " << alpha_rate << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "gamma = " << gamma << std::endl;
        //std::cout << "*********************************************" << std::endl;
        std::cout << std::endl;
        std::cout << "****** Break Down of Levels ******" << std::endl;
        for (int i = 0; i < L + 1; i++){
          std::cout << "Level " << i << "\t N_" << i << "= " << Q.getN(i) * nproc 
          << "\t E[Y_" << i << "] = " <<  myStats.getEQ(i) << "\t V[Y_" << i << "] = " << myStats.getVQ(i) 
          << "\t V[Y_" << i << "] = " << Q.ExpectedCost(i) << std::endl;
        }


      }



    }

    // line of best fit
double findRate(std::vector<double>& x, std::vector<double>& y){
  std::vector<double> logx(x.size()), logy(y.size());
  for (int i = 0; i < x.size(); i++){
    logx[i] = std::log(x[i]); logy[i] = std::log(y[i]);
  }
  double logx_bar = std::accumulate(logx.begin(),logx.end(),0.0) / (x.size());
  double logy_bar = std::accumulate(logy.begin(),logy.end(),0.0) / (y.size());
  double Cov = 0.0;
  double Var = 0.0;
  for (int i = 0; i < x.size(); i++){
    Cov += (logx[i] - logx_bar) * (logy[i] -logy_bar);
    Var += (logx[i] - logx_bar) * (logx[i] -logx_bar);
  }
  return Cov / Var;
}


private:

    M& model;
    std::vector<double> Cost;
    std::vector<double> Bias;
    std::vector<std::vector<double>> Y;

    double epsilon = 1e-2;
    int verbosity = 1;
    int initialSamples_Coarsest = 100;
    int adaptive;

};
