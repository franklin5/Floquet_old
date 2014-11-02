//! Base class declaration
/*!

  The base class sets up MPI distribution

  Created by Dong Lin on 10/22/14.
  Copyright (c) 2014 Dong Lin. All rights reserved.
*/

#ifndef DIST_H_
#define DIST_H_
#include "stdcpp.h"
class cDistribute {
 protected:
  int _rank, _size, _root;
  int recvcount, sendcount, stride;
  int *sendbuf, *recvbuf;
  int *sendcounts, *displs;
  int *recvcounts, *displs_r, offset;
  double  *localEig, *TotalEig;
 public:
  cDistribute(const int rank, const int size, const int root)
    {_rank = rank;_size=size;_root=root;}
  ~cDistribute(){
    if (_root==_rank) {
      delete []sendbuf;
      delete []sendcounts;
      delete []displs;
      delete []recvcounts;
	  delete []displs_r;
	  delete []TotalEig;
    }
    delete []recvbuf;
	delete []localEig;
  }
  int compute_count(int rank, int size, int NJOBS);// compute distribution count: the last rank does the remaineder job while the rest do the most even work.
  // NJOBS is the total number of work to be distributed, because each work is indepedent with each other.
  void distribution(int NJOBS);// The paradigm is to send momentum value to different rank
  // and let each processor work on the big matrix independently.
  // --> A two-level parallization is to be investigated in the furture
  //     for a sparse matrix partial eigenvalue spectrum calculation.
  void print_rank(); 

};

#endif // DIST_H_
