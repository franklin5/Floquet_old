//
//  floquet.h
//  Base class declaration
//  
//
//  Created by Dong Lin on 10/22/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
//

#ifndef FLOQUET_H_
#define FLOQUET_H_
#include "stdcpp.h"
class cFloquet {
 protected:
  int  _PMAX, _NMAX, _SMAX, _pblock, _NKX;
  double* _gauss_k, *_gauss_w_k;
  VectorXd  _bdg_E;
  MatrixXcd _bdg_V, _bdg_H;
 public:
  cFloquet(){}
  ~cFloquet(){}
 private:

};

#endif // FLOQUET_H_
