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
  int  _PMAX, _SMAX, _pblock, _NKX;
  double* _gauss_k, *_gauss_w_k;
  VectorXd  _bdg_E;
  MatrixXcd _bdg_V, _bdg_H;
  complex<double> _myI;
 public:
  cFloquet(){_myI = complex<double> (0.0,1.0);}
  ~cFloquet(){
	  delete []_gauss_k;
	  delete []_gauss_w_k;
  }
 private:

};

#endif // FLOQUET_H_
