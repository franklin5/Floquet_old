//! Abstract base class cBdG
/*!
 The base class sets up virtual process of reading input parameters
      and common components of matrix construction
	@author Lin Dong
	@date Oct 2014

  Created by Dong Lin on 10/22/14.
  Copyright (c) 2014 Dong Lin. All rights reserved.
*/

#ifndef BDG_H_
#define BDG_H_
#include "stdcpp.h"
#include "floquet.h"
#include "dist.h"
#include "lgwt.h"
class cBdG: public cFloquet {
 protected:
	int _pblock4, _ibdg;
	double  _h, _mu, _T, _Delta0, _v, _kmax;
 public:
  virtual void file_input()   =0; // this->file_input()
  virtual void file_output()   =0; // this->file_output()
  virtual void construction() =0; // this->construction()
  virtual void compute() =0;		// this->compute()
};

class cBdG_Bulk : public cBdG, public cDistribute{
// Bulk property: Chern number
// Wavefunction is expanded in frequency domain for time-dependent problem;
// p = 0 for time-independent problem.

private:
	int _NKX2;
	int _lowerbound, _upperbound;
	double *curvature_rank, *curvature;
	double _temp_curv;
	complex<double> _chern;
	double _chern_rank, _total_chern;
	char* _chernSolver;

public:
	cBdG_Bulk (const int rank, const int size, const int root) : cDistribute(rank,size,root){
		_chernSolver = new char [100];
	}
	~cBdG_Bulk(){
		delete []_chernSolver;
		if (_root==_rank) {
			delete []curvature;
		}
		delete []curvature_rank;
	}
	void file_input();
	void file_output();
	void construction();
	void update(int nk, double kx, double ky);
	void compute();
	double chern(int nk, double kx, double ky);
	void BrillouinZone();
};

class cBdG_Edge : public cBdG, public cDistribute{
// Edge property: spectrum under hard wall boundary condition.
// Wavefunction is expanded in both y-direction and frequency domain for time-dependent problem;
// p = 0 for time-independent problem.
private:
	int _NMAX;
	double _L;
	double  *localEig, *TotalEig;
public:
	cBdG_Edge (const int rank, const int size, const int root) : cDistribute(rank,size,root){}
	~cBdG_Edge(){}
	void file_input();
	void file_output();
	void construction();
	void update(int nkx);
	void compute();
};
#endif // BDG_H_
