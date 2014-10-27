/*! prx.h
 *
 * The base class sets up virtual process of reading input parameters
//      and common components of matrix construction
 *  Created on: Oct 24, 2014
 *      @Author Lin Dong
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
// Bulk property: Chern number
// Wavefunction is expanded in frequency domain for time-dependent problem;
// p = 0 for time-independent problem.
private:
	int _NKX2;
	double _kmax;
public:
	cPRX_Bulk(const int rank, const int size, const int root) : cDistribute(rank,size,root){}
	void file_input();
	void construction();


};

class cPRX_Edge : public cPRX, public cDistribute{
// Edge property: spectrum under hard wall boundary condition.
// Wavefunction is expanded in both y-direction and frequency domain for time-dependent problem;
// p = 0 for time-independent problem.
private:
	double _L;
public:
	cPRX_Edge(const int rank, const int size, const int root) : cDistribute(rank,size,root){}
	void file_input();
	void construction();


};


#endif /* PRX_H_ */
