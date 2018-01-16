

template <class M>
class testSamples{

public:

    testSamples (M& model_):model(model_)
    {

      std::cout << "Testing Class initialised!" << std::endl;
    }


    void inline apply(int l){

      RandomField z;

      int rank;

      MPI_Comm new_comm;

      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      int N = config.get<int>("MLMC.N",32);

      bool isTest = false;

      ofstream myfile;
      myfile.open("Results/results_"+std::to_string(rank)+".txt");

      for (int i = 0; i < N; i++){

        // Generate New Random Field

        z.user_random_field(isTest);

        double Q = model.getSample(0,z);

        if (rank == 0){std::cout << "Sample #" << i << " = " << Q << std::endl;}

      }

      myfile.close();


    
    /*  int N = 1000;
      std::vector<double> Q(N);

      model_.grid.globalRefine(l);


      for (int i = 0; i < N; i++)
      {
        // Generate Random Sample
        z.user_random_field();

        //  if (i > 0){model_.grid.globalRefine(l);}
          Q[i] = model_.getSample(l,z);
      }
      std:cout << "=====" << std::endl;
      for(int i = 0; i < N;i++){
        std::cout << Q[i] << std::endl;
      }*/

   /*   int maxL = 6;
      int N = 1000;

     std::vector<double> Q(maxL), Cost(maxL);

      model_.grid.globalRefine(maxL);

      std::cout << "Grid About to be refined" << std::endl;

      std::cout << "Refined Grid" << std::endl;

      int rank;

      MPI_Comm new_comm;

      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      ofstream myfile;
      myfile.open("results_"+std::to_string(rank)+".txt");
        
        ofstream myfile2;
        myfile2.open("results_costs_"+std::to_string(rank)+".txt");


      for (int j = 0; j < N; j++){

      z.user_random_field();

      std::cout << j << std::endl;

      Dune::Timer timer;

      for (int i = 0; i < maxL; i++){

        timer.reset();
        Q[i] = model_.getSample(i,z);
        Cost[i] += timer.elapsed();

        if (i > 0){myfile << " "; myfile2 << " ";}
        myfile << Q[i];
        myfile2 << Cost[i];

        if (i == maxL-1){myfile << "\n"; myfile2 << "\n";}

      }

    }

    myfile.close();
        myfile2.close();

      std::cout << "Average Costs" << std::endl;

      for (int i = 0 ; i < maxL; i++){
        std::cout << Cost[i]/N << std::endl;
      }




*/



    }


private:

    M& model;
    std::vector<double> Cost;
    std::vector<int> Nopt, N;
    std::vector<double> Bias;
    std::vector<std::vector<double>> Y;

    double epsilon = 1e-2;
    int verbosity = 1;
    int initialSamples_Coarsest = 100;
    int adaptive;

};
