//
//  parallel.hh
//  dune-uq
//
//  Created by Tim Dodwell on 13/01/2016.
//  Copyright © 2016 Tim Dodwell. All rights reserved.
//

#ifndef parallel_h
#define parallel_h

void block_distribution(std::vector<int>& id, int rank, int nproc, int N){
    
    int Q = N / nproc;
    
    int r = N - Q * nproc;
    
    if (rank < r){
    id.resize(Q + 1);
        for (int i = 0; i < Q + 1; i++){
            id[i] = rank * (Q + 1) + i;
        }
    }
    else{
        id.resize(Q);
        for (int i = 0; i < Q ; i++){
            id[i] = r * (Q + 1) + (rank - r) * Q + i;
        }
    }
    
    
}


#endif /* parallel_h */
