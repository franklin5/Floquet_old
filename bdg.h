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
class cBdG: public cFloquet {
 protected:
  double  _h, _mu, _T, _Delta0, _v;
 public:
  virtual void file_input()   =0; // this->file_input()
  virtual void construction() =0; // this->construction()

};

class cBdG_Bulk : public cBdG, public cDistribute{
private:
	int _NKX2;
	double _kmax;
public:
	cBdG_Bulk (const int rank, const int size, const int root) : cDistribute(rank,size,root){}
	void file_input();
	void construction();

};

class cBdG_Edge : public cBdG, public cDistribute{
private:
	double _L;
public:
	cBdG_Edge (const int rank, const int size, const int root) : cDistribute(rank,size,root){}
	void file_input();
	void construction();


};
#endif // BDG_H_
