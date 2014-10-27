//                                                                            
//  bdg.h
//  Abstract base class cBdG
//  The base class sets up virtual process of reading input parameters
//      and common components of matrix construction
//  Created by Dong Lin on 10/22/14.                                           
//  Copyright (c) 2014 Dong Lin. All rights reserved.                          
//   

#ifndef BDG_H_
#define BDG_H_
#include "stdcpp.h"
#include "floquet.h"
#include "dist.h"
#include "lgwt.h"
class cBdG: public cFloquet {
 protected:
	int _pblock4, _ibdg;
	double  _h, _mu, _T, _Delta0, _v;
 public:
  virtual void file_input()   =0; // this->file_input()
  virtual void construction() =0; // this->construction()
  virtual void compute() =0;		// this->compute()
};

class cBdG_Bulk : public cBdG, public cDistribute{
// Bulk property: Chern number
// Wavefunction is expanded in frequency domain for time-dependent problem;
// p = 0 for time-independent problem.

private:
	int _NKX2;
	double _kmax, _temp_curv;
	complex<double> _chern;
	char* _chernSolver;
public:
	cBdG_Bulk (const int rank, const int size, const int root) : cDistribute(rank,size,root){
		_chernSolver = new char [100];
	}
	~cBdG_Bulk(){
		delete []_chernSolver;
	}
	void file_input();
	void construction();
	void update(int nk, double kx, double ky);
	void compute();
	double chern(int nk, double kx, double ky);
};

class cBdG_Edge : public cBdG, public cDistribute{
// Edge property: spectrum under hard wall boundary condition.
// Wavefunction is expanded in both y-direction and frequency domain for time-dependent problem;
// p = 0 for time-independent problem.
private:
	int _NMAX;
	double _L, _kmax;
public:
	cBdG_Edge (const int rank, const int size, const int root) : cDistribute(rank,size,root){}
	void file_input();
	void construction();
	void update(int nkx);
	void compute();

};
#endif // BDG_H_
