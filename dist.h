//
//  dist.h
//  Base class declaration
//  The base class sets up MPI distribution 
//  
//  Created by Dong Lin on 10/22/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
//

#ifndef DIST_H_
#define DIST_H_
#include "stdcpp.h"
class cDistribute {
 protected:
  int _rank, _size, _root;
  int recvcount, sendcount, stride;
  int *sendbuf, *recvbuf;
  int *sendcounts, *displs;
 public:
  cDistribute(const int rank, const int size, const int root)
    {_rank = rank;_size=size;_root=root;}
  ~cDistribute(){
    if (_root==_rank) {
      delete []sendbuf;
      delete []sendcounts;
      delete []displs;
    }
      delete []recvbuf;
  }
  int compute_count(int,int,int);
  void distribution(int);
  void print_rank(); 

};

#endif // DIST_H_
