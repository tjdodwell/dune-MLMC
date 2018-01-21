
        MPI_Barrier(MPI_COMM_WORLD); // Wait for all processors to complete samples

         //  *** Compute PayOff, Bias and Variance of Estimator on each level on rank 0

        Q.getEQs(localEQ);  Q.getEQ2s(localEQ2);
        
        double * EQRoot = NULL;
        double * EQ2Root = NULL;

        if (rank == 0){
                EQRoot = new double[nproc];
                EQ2Root = new double[nproc];
        }

        // Gather all samples from each processor to root processor
        
        MPI_Gather(&localEQ[L],1,MPI_DOUBLE,EQRoot,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        

        double Bias = 0.0;
        double Var = 0.0;

        std::vector<double> localPayOff(L + 1), VarEstim_l(L + 1);

        for (int l = 0; l < L + 1; l++){ // For each level

          // Communicate EQ and EQ2 from each processor to rank 0 processor
          MPI_Gather(localEQ[l],1,MPI_DOUBLE,EQRoot,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
          MPI_Gather(localEQ2[l],1,MPI_DOUBLE,EQ2Root,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

          if (rank == 0){

            // * First compute EQ and E[Q^2] - from those collected from each processor
            double EQ_all = std::accumulate(EQRoot.begin(), EQRoot.end(), 0.0) / nproc;
            double Var_all = (std::accumulate(EQ2Root.begin(), EQ2Root.end(), 0.0) / nproc) - EQ_all * EQ_all;
            
            VarEstim_l[l] =  Var_all / (Q.numSamples(l) * nproc); // Variance of Estimate

            PayOff[l] = VarEstim_l[l] / (Q.ExpectedCost(l) *; // Payoff on level l

            if (l == L){ // if level is L then compute Bias;
              Bias = EQ_all / (std::pow(4.0,alpha_rate));
            } // end if for level is L

          } // end if rank == 0

        } // end for each level

        VoE = std::accumulate(Var_Estim_l.begin(), Var_Estim_l.end(),0.0);

        MPI_Bcast(&PayOff, L + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // BroadCast PayOff Vector
        MPI_Bcast(Bias,1,MPI_DOUBLE,0,MPI_COMM_WORLD); // BroadCast Bias
        MPI_Bcast(VoE,1,MPI_DOUBLE,0,MPI_COMM_WORLD); // BroadCast VoE


      }


        while ( VoE > (std::sqrt(theta) * e2)){

          // Find level with maximum payoff
          
          std::vector<double>::iterator maxPayOffLevel = std::max_element(PayOff.begin(),PayOff.end());
          
          int greedyLevel = std::distance(PayOff.begin(),maxPayOffLevel);

          // Compute new greedy samples

          int newSamples = std::ceil(greedyFactor * Q.getN(greedyLevel) / nproc);
         
          for (int i = 0; i < newSamples; i++){
             
             watch.reset();
             
             double QoI = model.getSample(greedyLevel,z);
             
             Q.addSample(QoI,greedyLevel,watch.elapsed());
          }

          Q.UpdateStats(); // Update Local Statistics


          MPI_Barrier(MPI_COMM_WORLD); // Wait for all processors to complete samples

          // Recompute Variance of Estimator and PayOff


                  for (int l = 0; l < L + 1; l++){ // For each level

          // Communicate EQ and EQ2 from each processor to rank 0 processor
          MPI_Gather(localEQ[l],1,MPI_DOUBLE,EQRoot,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
          MPI_Gather(localEQ2[l],1,MPI_DOUBLE,EQ2Root,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

          if (rank == 0){

            // * First compute EQ and E[Q^2] - from those collected from each processor
            double EQ_all = std::accumulate(EQRoot.begin(), EQRoot.end(), 0.0) / nproc;
            double Var_all = (std::accumulate(EQ2Root.begin(), EQ2Root.end(), 0.0) / nproc) - EQ_all * EQ_all;
            
            
            VarEstim_l[l] =  Var_all / (Q.numSamples(l) * nproc); // Variance of Estima

            PayOff[l] = VarEstim_l[l] / (Q.ExpectedCost(l) *; // Payoff on level l

            if (l == L){ // if level is L then compute Bias;
              Bias = EQ_all / (std::pow(4.0,alpha_rate));
            } // end if for level is L

            } // end if rank == 0

          } // end for each level

          VoE = std::accumulate(Var_Estim_l.begin(), Var_Estim_l.end(),0.0);

          MPI_Bcast(&PayOff, L + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // BroadCast PayOff Vector
          MPI_Bcast(VoE,1,MPI_DOUBLE,0,MPI_COMM_WORLD); // BroadCast VoE


        } // End while greedy part of algorithm

          // Recompute Bias

          Bias = Q.computeBias(alpha_rate);

          // Communicate EQ and EQ2 from each processor to rank 0 processor
          MPI_Gather(localEQ[L],1,MPI_DOUBLE,EQRoot,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
          MPI_Gather(localEQ2[L],1,MPI_DOUBLE,EQ2Root,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

          if (rank == 0){

            // * First compute EQ and E[Q^2] - from those collected from each processor
            double EQ_all = std::accumulate(EQRoot.begin(), EQRoot.end(), 0.0) / nproc;
            
            
            Bias = EQ_all / (std::pow(4.0,alpha_rate));
          
          }

          MPI_Bcast(Bias,1,MPI_DOUBLE,0,MPI_COMM_WORLD); // BroadCast Bias

          MPI_Barrier(MPI_COMM_WORLD); // Wait for all processors to complete samples



          if(rank == 0 && config.get<bool>("print.inSim",true)){
            std::cout << "Current Level = " << L << " Bias = " << Bias << " Variance of Estimator = " << VoE << std::endl;
          }

          if (Bias * Bias < (1.0 - theta) * e2){
            converged = true;
          }
