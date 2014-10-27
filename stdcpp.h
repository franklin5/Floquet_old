/*! handy standard header files
  stdcpp.h
//
//
//  Created by Dong Lin on 9/2/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
*/

#ifndef STDCPP_H_
#define STDCPP_H_
#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <new>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <Eigen/Eigenvalues>
#include <time.h>
using namespace std;
using namespace Eigen;
using namespace MPI;
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif
