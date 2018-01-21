


class SampleContainer_TestCase2{
	
public:

	SampleContainer_TestCase2(){
		L = -1;
	};

	void inline addLevel(){
		L += 1;
		// Resize Vectors
		QoI.resize(L + 1);	QoI2.resize(L + 1); Cost.resize(L + 1);	EC.resize(L + 1);	
		EQ.resize(L + 1);	EQ2.resize(L + 1);
		Var.resize(L + 1);	N.resize(L + 1);
		pp.resize(L+1); 	pm.resize(L+1);
	}


	void inline addSample(double Q, int level, double time){

		if (Q > 0.0){ pp[level] += 1;}	if (Q < 0.0){ pm[level] += 1;}
	
		QoI[level].push_back(Q);
		QoI2[level].push_back(Q * Q);
		Cost[level].push_back(time);
		N[level] += 1;
	}


	void inline UpdateStats(){

		for (int l = 0; l < L + 1; l++){

			if(config.get<bool>("MLMC.biasEstimator",true)){

				double k = config.get<double>("MLMC.biasEstimator_k",1.0);

				double pmt = ((double)pm[l] + k) / ((double)N[l] + k);
				double ppt = ((double)pp[l] + k) / ((double)N[l] + k);

				EQ[l] = std::abs(ppt - pmt);
				Var[l] = ppt + pmt - (ppt - pmt) * (ppt - pmt);

			}
			else{

				EQ[l] = std::accumulate(QoI[l].begin(),QoI[l].end(),0.0) / QoI[l].size();
				EQ2[l] = std::accumulate(QoI2[l].begin(),QoI2[l].end(),0.0) / QoI2[l].size();
				Var[l] = EQ2[l] - EQ[l] * EQ[l];

			}



			EC[l] = std::accumulate(Cost[l].begin(),Cost[l].end(),0.0) / Cost[l].size();
			
		}
	}

	double inline MLMCestimator(){
		return std::accumulate(EQ.begin(),EQ.end(),0.0);
	}

	double inline computeBias(double alpha_rate){
		return std::abs(EQ[L] / (std::pow(4,alpha_rate) - 1.0));
	}

	double inline VarEstim(int level){ return Var[level] / N[level];}

	double inline computeVarianceOfEstimator(){
		double VoE = 0.0;
		for (int i = 0; i < L+1; i++){
			VoE += Var[i] / N[i];
		}
		return VoE;
	}

	int inline getN(int level){return QoI[level].size(); }
	double inline expectedValue(int level){return EQ[level];}
	double inline Variance(int level){return Var[level];}
	double inline ExpectedCost(int level){ return EC[level]; }
	
private:

	int L;

	std::vector<std::vector<double>> QoI, QoI2, Cost;

	std::vector<int> pp, pm;

	std::vector<int> N;

	std::vector<double> EQ, EQ2, Var, EC;


};


