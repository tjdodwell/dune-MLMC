//
//  KLFunctions.h
//
//
//  Created by Tim Dodwell on 01/12/2015.
//
//

#ifndef KLFunctions_h
#define KLFunctions_h

// ------------------------------------------------------------------------------------
// Define K-L Expansion
// ------------------------------------------------------------------------------------


template<class X>
void inline evaluate_eigenValues(double ell, double sig_KL, X& lambda, X& freq){
  // Checked TD 14/11/2016
    double c = 1 / ell;
    for (int i = 0; i < lambda.size(); i++){
        lambda[i] = 2.0 * c /(freq[i] * freq[i] + c * c);
    }
}

template<class X,class Y>
void bubbleSort(X& index, Y& lam)
{
    const int N = lam.size();

    for (int i = 0; i < N; i++){ index[i] = i; }

    int sw = 1; double tmpreal;
    int tmpint;

    while (sw == 1)
    {
        sw = 0;
        for(int i = 0; i < N-1; i++)
        {
            if (lam[i] < lam[i+1])
            {
                tmpreal = lam[i+1];
                tmpint = index[i+1];
                lam[i+1] = lam[i];
                index[i+1] = index[i];
                lam[i] = tmpreal;
                index[i] = tmpint;
                sw = 1;
            }
        }
    }
}

template<class X, class Y, class Z>
void inline construct_3D_eigenValues(Z& lam1D,X& lambda3D, Y& id2d_i, Y& id2d_j, Y& id2d_k){
    const int N = lam1D.size();

    std::vector<int> index_i(N * N * N), index_j(N * N * N), index_k(N * N * N);
    std::vector<double> lam3D(lam1D.size() * lam1D.size() * lam1D.size());
    std::vector<int> ind(lam1D.size() * lam1D.size() * lam1D.size());

    int counter = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++){

            for (int k = 0; k < N; k++){

                lam3D[counter] = lam1D[i] * lam1D[j] * lam1D[k];

                ind[counter] = counter;

                index_i[counter] = i;
                index_j[counter] = j;
                index_k[counter] = k;

                counter += 1;

            }

        }
    }

    bubbleSort(ind,lam3D);

    for (int i = 0; i < lambda3D.size(); i++){
        id2d_i[i] = index_i[ind[i]];
        id2d_j[i] = index_j[ind[i]];
        id2d_k[i] = index_k[ind[i]];
        lambda3D[i] = lam3D[i];

    }

}




template<class X, class Y, class Z>
void inline construct_2D_eigenValues(Z& lam1Dx, Z& lam1Dy, X& lambda2D, Y& id2d_i, Y& id2d_j){
    const int N = lam1Dx.size();

    std::vector<int> index_i(N * N),    index_j(N * N);
    std::vector<double> lam2D(N * N);
    std::vector<int> ind(N * N);

    int counter = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++){

            lam2D[counter] = lam1Dx[i] * lam1Dy[j];

            ind[counter] = counter;

            index_i[counter] = i;
            index_j[counter] = j;

            counter += 1;

        }
    }

    bubbleSort(ind,lam2D);

    for (int i = 0; i < lambda2D.size(); i++){
        id2d_i[i] = index_i[ind[i]];
        id2d_j[i] = index_j[ind[i]];
        lambda2D[i] = lam2D[i];

    }

}


template<class X, class Y, class Z>
void inline construct_2D_eigenValues(Z& lam1Dx, X& lambda2D, Y& id2d_i, Y& id2d_j){
    const int N = lam1Dx.size();

    Z lam1Dy = lam1Dx;

    std::vector<int> index_i(N * N),    index_j(N * N);
    std::vector<double> lam2D(N * N);
    std::vector<int> ind(N * N);

    int counter = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++){

            lam2D[counter] = lam1Dx[i] * lam1Dy[j];

            ind[counter] = counter;

            index_i[counter] = i;
            index_j[counter] = j;

            counter += 1;

        }
    }

    bubbleSort(ind,lam2D);

    for (int i = 0; i < lambda2D.size(); i++){
        id2d_i[i] = index_i[ind[i]];
        id2d_j[i] = index_j[ind[i]];
        lambda2D[i] = lam2D[i];

    }

}



double res(double x, double ell){
  double c = 1.0 / ell;
  double g = std::tan(x) - 2.0 * c * x / (x * x - c * c);
  return g;
}


double findRoot(double ell, double a, double b){
  double fa, fb, x, fx, error, m;
  error = 1.0;
  while (error > 1e-6){
    fa = res(a,ell);
    fb = res(b,ell);
    m = (fb - fa) / (b - a);
    x = a - fa / m;
    fx = res(x,ell);
    if (((fa < 0) & (fx < 0)) | ((fa > 0) & (fx > 0))) { a = x; }
    else { b = x; }
    error = std::abs(fx);
  }
  return x;
}

void rootFinder(int M, double ell, std::vector<double>& answer){
double c = 1.0 / ell;
std::vector<double> freq(M+2);
// For all intervals

int m = -1;
for (int i = 0; i < M + 1; i++){

//  std::cout << "Root i = " << i << std::endl;

  double w_min = (i - 0.4999) * M_PI;
  double w_max = (i + 0.4999) * M_PI;

//  std::cout << "w_min = " << w_min << std::endl;
  //std::cout << "w_max = " << w_max << std::endl;
  if ((w_min <= c) && (w_max >= c)){

    // If not first interval look for solution near left boundary
    if (w_min > 0.0){
      m += 1;
      freq[m] = findRoot(ell,w_min,0.5*(c+w_min));
    }
    // Always look for solution near right boundary
    m += 1;
    freq[m] = findRoot(ell,0.5*(c + w_max),w_max);
  }
  else{
    m += 1;
    freq[m] = findRoot(ell,w_min,w_max);
  }

}

for (int i = 0; i < M; i++){
  answer[i] = freq[i+1];
}



}



#endif /* KLFunctions_h */
