/*
 * prx.h
 * The base class sets up virtual process of reading input parameters
//      and common components of matrix construction
 *  Created on: Oct 24, 2014
 *      Author: ld7
 */

#ifndef PRX_H_
#define PRX_H_
#include "stdcpp.h"
#include "floquet.h"
#include "dist.h"
class cPRX: public cFloquet {
 protected:
  double _J, _b, _a, _omega;
 public:
  virtual void file_input()   =0; // this->file_input()
  virtual void construction() =0; // this->construction()

};

class cPRX_Bulk : public cPRX, public cDistribute{
private:
	int _NKX2;
	double _kmax;
public:
	cPRX_Bulk(const int rank, const int size, const int root) : cDistribute(rank,size,root){}
	void file_input();
	void construction();


};

class cPRX_Edge : public cPRX, public cDistribute{
private:
	double _L;
public:
	cPRX_Edge(const int rank, const int size, const int root) : cDistribute(rank,size,root){}
	void file_input();
	void construction();


};


#endif /* PRX_H_ */
