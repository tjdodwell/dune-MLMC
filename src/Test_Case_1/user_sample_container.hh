


class SampleContainer_TestCase1{
	
public:

	SampleContainer_TestCase1(){
		L = -1;
	};

	void inline addLevel(){
		L += 1;
		// Resize Vectors
		QoI.resize(L + 1);	QoI2.resize(L + 1); Cost.resize(L + 1);	EC.resize(L + 1);	
		EQ.resize(L + 1);	EQ2.resize(L + 1);
		Var.resize(L + 1);	N.resize(L + 1);
	}


	void inline addSample(double Q, int level, double time){
	
		QoI[level].push_back(Q);
		QoI2[level].push_back(Q * Q);
		Cost[level].push_back(time);
		N[level] += 1;
	}


	void inline UpdateStats(){
		for (int l = 0; l < L + 1; l++){
			EC[l] = std::accumulate(Cost[l].begin(),Cost[l].end(),0.0) / Cost[l].size();
			EQ[l] = std::accumulate(QoI[l].begin(),QoI[l].end(),0.0) / QoI[l].size();
			EQ2[l] = std::accumulate(QoI2[l].begin(),QoI2[l].end(),0.0) / QoI2[l].size();
			Var[l] = EQ2[l] - EQ[l] * EQ[l];
		}
	}

	double inline numSamples(int l){return N[l];}

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

	void inline getEQs(std::vector<double>& vec){
		for (int i = 0; i < EQ.size(); i++){
			vec[i] = EQ[i];
		}
	}

	void inline getEQ2s(std::vector<double>& vec){
		for (int i = 0; i < EQ2.size(); i++){
			vec[i] = EQ2[i];
		}
	}
	
private:

	int L;

	std::vector<std::vector<double>> QoI, QoI2, Cost;

	std::vector<int> N;

	std::vector<double> EQ, EQ2, Var, EC;


};


