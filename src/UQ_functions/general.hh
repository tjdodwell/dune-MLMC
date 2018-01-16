//
//  general.hh
//  
//
//  Created by Tim Dodwell on 30/11/2015.
//
//

#ifndef general_h
#define general_h

double mean(std::vector<double>& Y){
    double Yhat = 0.0;
    for (int i = 0; i < Y.size(); i++){
        Yhat += Y[i];
    }
    Yhat /= Y.size();
    return Yhat;
}



double var(std::vector<double>& Y){
    
    double Yhat = mean(Y);
    double variance = 0.0;
    int N = Y.size();
    
    for (int i = 0; i < N; i++){
        variance += (Y[i] - Yhat) * (Y[i] - Yhat);
    }
    
    variance /= (N - 1);
    return variance;
}




#endif /* general_h */
