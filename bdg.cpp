/*! \file
 * bdg.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: ld7
 */
#include "bdg.h"
void cBdG_Bulk::file_input(){
  for (int ig = 0; ig < _size; ++ig) {
    if (ig ==_rank){
      char dummyname[100];
      double dummyvalue;
      int intdummyvalue;
      FILE *input_bdg_bulk;
      input_bdg_bulk = fopen("input_bdg_bulk.txt","r");
      assert(input_bdg_bulk != NULL);
      if (ig == _root)  cout << "Starting to read in parameters from file input_bdg_bulk.txt" << endl;
      fscanf(input_bdg_bulk,"%s %d", dummyname, &intdummyvalue);
      _PMAX = intdummyvalue;	if (ig == _root) cout << dummyname << "=" << _PMAX << endl;
      fscanf(input_bdg_bulk,"%s %d", dummyname, &intdummyvalue);
      _NKX = intdummyvalue;	if (ig == _root) cout << dummyname << "=" << _NKX << endl;
      fscanf(input_bdg_bulk,"%s %lf", dummyname, &dummyvalue);
      _h = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _h << endl;
      fscanf(input_bdg_bulk,"%s %lf", dummyname, &dummyvalue);
      _mu = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _mu << endl;
      fscanf(input_bdg_bulk,"%s %lf", dummyname, &dummyvalue);
      _T = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _T << endl;
      fscanf(input_bdg_bulk,"%s %lf", dummyname, &dummyvalue);
      _Delta0 = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _Delta0 << endl;
      fscanf(input_bdg_bulk,"%s %lf", dummyname, &dummyvalue);
      _v = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _v << endl;
      fscanf(input_bdg_bulk,"%s %lf", dummyname, &dummyvalue);
      _kmax = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _kmax << endl;
      fscanf(input_bdg_bulk,"%s %s", dummyname, _chernSolver);
      if (ig == _root) cout << dummyname << "=" << _chernSolver  <<endl;
      fclose(input_bdg_bulk);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
void cBdG_Bulk::construction(){
  if (_rank == _root) cout << "Construction has started\n";
  _ibdg = 4; // --> It has to be 4 here. (at least currently.)
  _pblock = (2*_PMAX+1);
  _pblock4 = _ibdg*_pblock;
  _SMAX = _pblock4;
  _NKX2 = _NKX*_NKX;
  _bdg_E.resize(_SMAX);
  _bdg_V.resize(_SMAX,_SMAX);
  _bdg_H.resize(_SMAX,_SMAX);
  update(-1,0,0); // null construction that is only done once for arbitrary momentum state.
  _gauss_k = new double [_NKX];  // integration momentum value
  _gauss_w_k = new double [_NKX];// integration momentum weight
  gauss_lgwt(_NKX,-_kmax,_kmax,_gauss_k,_gauss_w_k); // TODO: option of doing 1D integral to be updated.
  // Legendre-Gauss quadrature method of integration:
  // code adapted from a matlab code in FileExchangeMATLAB: (lgwt.m)
  // http://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes/content/lgwt.m
  if (_rank == _root) cout << "Construction is completed." << endl;
}
void cBdG_Bulk::update(int nk, double kx, double ky){
  int p,q;
  if (nk == -1) {
    _bdg_H.setZero(); // This is done only once.
    // The off-diagonal coupling introduced from time-dependent order parameter should be computed only here.
    double dt = 0.0001;int lenT = int(_T/dt);double Delta2 = 0.0;
    double *DELTA_t = new double [lenT];
    for (int it = 0; it < lenT; ++it) {
      DELTA_t[it] = _Delta0 + Delta2*cos(it*dt*2*M_PI/_T); // time-independent problem: replace Delta2 with 0.0;
    }
    complex<double> Gamma2;
    for (int ip = 0; ip < _pblock; ++ip) {
      p = ip - _PMAX;
      for (int iq = 0; iq<=ip;++iq){
	q = iq - _PMAX;
	Gamma2 = complex<double> (0.0,0.0);
	for (int it = 0; it < lenT; ++it) {
	  Gamma2 += DELTA_t[it]*exp(-_myI*2*M_PI*(q-p)/_T*it*dt)/_T*dt;
	}
	_bdg_H(ip*_ibdg+2,iq*_ibdg+1)   =  Gamma2;
	_bdg_H(ip*_ibdg+3,iq*_ibdg) 	= -Gamma2;
      }
    }
    delete []DELTA_t;
  } else {
    double xi = kx*kx+ky*ky-_mu;
    for (int ip = 0; ip < _pblock; ++ip) {
      p = ip - _PMAX;
      _bdg_H(ip*_ibdg,ip*_ibdg)   		=  complex<double> (xi+_h+2*p*M_PI/_T,0.0);
      _bdg_H(ip*_ibdg+1,ip*_ibdg)   	=  complex<double> (_v*kx,_v*ky);
      _bdg_H(ip*_ibdg+1,ip*_ibdg+1)   	=  complex<double> (xi-_h+2*p*M_PI/_T,0.0);
      _bdg_H(ip*_ibdg+2,ip*_ibdg+2)   	=  complex<double> (-(xi+_h)+2*p*M_PI/_T,0.0);
      _bdg_H(ip*_ibdg+3,ip*_ibdg+2)   	=  complex<double> (_v*kx,-_v*ky);
      _bdg_H(ip*_ibdg+3,ip*_ibdg+3)   	=  complex<double> (-(xi-_h)+2*p*M_PI/_T,0.0);
    }
    //cout << _bdg_H << endl;
  }
}
void cBdG_Bulk::BrillouinZone(){
  _lowerbound = -999;
  _upperbound = -999; // ridiculous negative flag
  SelfAdjointEigenSolver<MatrixXcd> ces;
  ces.compute(_bdg_H);
  _bdg_E = ces.eigenvalues();
  _bdg_V = ces.eigenvectors();
  if (_PMAX == 0) {
    _lowerbound = 0;
    _upperbound = _SMAX-1;
  } else {
    for(int ip = 0; ip < _SMAX/2;++ip){
      if (_bdg_E[ip]/(M_PI/_T) >= -1.0) {
	_lowerbound = ip;
	break;
      }
    }
    for(int ip = _SMAX/2; ip < _SMAX;++ip){
      if (_bdg_E[ip]/(M_PI/_T) >= 1.0) {
	_upperbound = ip-1;
	break;
      }
    }
  }
}
double cBdG_Bulk::chern(int nk, double kx, double ky){
  double result;
  complex<double> u,a,b,v,up,ap,bp,vp, Theta1,Theta2, temp;
  if (_lowerbound < 0 || _upperbound < 0){
    _temp_curv = 0.0;
    _chern = complex<double> (0.0,0.0);
    //      cout << "no contribution is added" << endl;
  } else {
    if (string(_chernSolver)=="curvature") {
      //      cout  <<"lower bound = " << lowerbound << " upper bound = " << up	\
      perbound << ", and " << upperbound-lowerbound+1 << " is considered for computation." <<endl;
      double Ediff;
      int nkx = nk % _NKX, nky = int (nk/_NKX);
      _chern = complex<double> (0.0,0.0);
      for(int ih = _lowerbound; ih < _SMAX/2; ++ih) { // hole branch
	for(int ip = _SMAX/2;ip<=_upperbound;++ip){ // particle branch
	  Theta1 = complex<double> (0.0,0.0);
	  Theta2 = complex<double> (0.0,0.0);
	  for(int i = 0; i < _pblock; ++i){ // frequency block adds up
	    u = _bdg_V(i*_ibdg,ih);
	    a = _bdg_V(i*_ibdg+1,ih);
	    b = _bdg_V(i*_ibdg+2,ih);
	    v = _bdg_V(i*_ibdg+3,ih);
	    up = _bdg_V(i*_ibdg,ip);
	    ap = _bdg_V(i*_ibdg+1,ip);
	    bp = _bdg_V(i*_ibdg+2,ip);
	    vp = _bdg_V(i*_ibdg+3,ip);
	    Theta1 += 2*kx*up*conj(u)+_v*ap*conj(u)
	      +_v*up*conj(a)+2*kx*ap*conj(a)
	      -2*kx*bp*conj(b)+_v*vp*conj(b)
	      +_v*bp*conj(v)-2*kx*vp*conj(v);
	    Theta2 += 2*ky*conj(up)*u-_myI*_v*conj(up)*a
	      +_myI*_v*conj(ap)*u+2*ky*conj(ap)*a
	      -2*ky*conj(bp)*b+_myI*_v*conj(bp)*v
	      -_myI*_v*conj(vp)*b-2*ky*conj(vp)*v;
	  }
	  Ediff = _bdg_E[ih]-_bdg_E[ip];
	  _chern += -2.0*Theta1*Theta2/(Ediff*Ediff);
	}
      }
      _temp_curv = _chern.imag();
      _chern = _chern*_gauss_w_k[nkx]*_gauss_w_k[nky]/(2.0*M_PI);
      result = _chern.imag();
    } else if (string(_chernSolver)=="connection") {
      double qx, qy;
      SelfAdjointEigenSolver<MatrixXcd> ces;
      double dk = 0.5 * _kmax * 2.0 /(_NKX-1);
      MatrixXcd loopA,loopB,loopC,loopD;
      loopA.resize(_SMAX,_SMAX);loopB.resize(_SMAX,_SMAX);loopC.resize(_SMAX,_SMAX);loopD.resize(_SMAX,_SMAX);
      for (int i = 1; i < 5; ++i) {
	switch (i) {
	case 1:
	  qx = kx - dk;
	  qy = ky - dk;
	  update(i, qx, qy);
	  ces.compute(_bdg_H);
	  loopA = ces.eigenvectors();
	  //		cout << "loopA is finished." << endl;
	  break;
	case 2:
	  qx = kx + dk;
	  qy = ky - dk;
	  update(i, qx, qy);
	  ces.compute(_bdg_H);
	  loopB = ces.eigenvectors();
	  //		cout << "loopB is finished." << endl;
	  break;
	case 3:
	  qx = kx + dk;
	  qy = ky + dk;
	  update(i, qx, qy);
	  ces.compute(_bdg_H);
	  loopC = ces.eigenvectors();
	  //		cout << "loopC is finished." << endl;
	  break;
	case 4:
	  qx = kx - dk;
	  qy = ky + dk;
	  update(i, qx, qy);
	  ces.compute(_bdg_H);
	  loopD = ces.eigenvectors();
	  //		cout << "loopD is finished." << endl;
	  break;
	default:
    	  break;
	}
      }
      _chern = complex<double> (1.0,0.0);
      complex<double> temp;
      for(int ip = _lowerbound; ip < _SMAX/2; ++ip) {
	temp = loopA.col(ip).adjoint() * loopB.col(ip);
	_chern = _chern * temp/abs(temp);
	temp = loopB.col(ip).adjoint() * loopC.col(ip);
	_chern = _chern * temp/abs(temp);
	temp = loopC.col(ip).adjoint() * loopD.col(ip);
	_chern = _chern * temp/abs(temp);
	temp = loopD.col(ip).adjoint() * loopA.col(ip);
	_chern = _chern * temp/abs(temp);
      }
      _chern = complex<double>(log(std::abs(_chern)),atan(_chern.imag()/_chern.real()));
      _chern = _chern/(2.0*M_PI)*_myI;
      result = _chern.real();
    } else{
      cerr << "Chern number calculation method is not supported. Please check input_bdg_bulk.txt and readme." << endl;
      exit(1);
    }
  }
  return result;
}

void cBdG_Bulk::compute(){
  distribution(_NKX2);
  stride = _SMAX*recvcount;
  localEig = new double[stride];
  double kx, ky;
  int nk, nkx, nky;
  for(int ig = 0; ig< _size; ++ig) {
    if (ig == _rank){
      cout << "rank " << ig << " has started"<< endl;
      _chern_rank = 0.0;
      if (string(_chernSolver)=="curvature") curvature_rank = new double[recvcount];
      for (int i=0; i<recvcount; ++i) { // recvcount is the number of tasks given to each rank.
	clock_t start = clock();
		nk = recvbuf[i]; // distributed momentum value
		if (string(_chernSolver)=="curvature") {
		  // gauss quadrature is good for integration
		  kx = _gauss_k[nk % _NKX];
		  ky = _gauss_k[int (nk/_NKX)];
		} else if (string(_chernSolver)=="connection") {
		  // in connection approach, we have to be wise and cover all the areas that are under consideration
		  // square area with the same size is easier to work with, rather than the quadrature meshes.
		  nkx = nk % _NKX;
		  nky = int (nk/_NKX);
		  kx = -_kmax + nkx * _kmax *2.0 /(_NKX-1);
		  ky = -_kmax + nky * _kmax *2.0 /(_NKX-1);
		}
		update(nk, kx, ky);
		BrillouinZone();
		for(int j = 0; j < _SMAX; ++j){
		  localEig[j+i*_SMAX]=_bdg_E[j];
		}
		_chern_rank += chern(nk, kx,ky);
	clock_t end = clock();
	if (_rank==_root) cout << "task " << recvbuf[i] <<"out of " << recvcount << "used " << double (end-start)/ (double) CLOCKS_PER_SEC  << endl;
	if (string(_chernSolver)=="curvature") curvature_rank[i] = _temp_curv;
      }
      cout << "rank " << _rank << " has finished "<< recvcount << " tasks " << endl;
    }
  }
}
void cBdG_Bulk::file_output(){
  for(int ig = 0; ig<_size; ++ig) {
    if (ig ==_rank){
      cout << "rank " << ig << " has chern number "<< _chern_rank << endl;
    }
    MPI_Barrier(COMM_WORLD);
  }
  MPI_Reduce(&_chern_rank, &_total_chern, 1, MPI_DOUBLE, MPI_SUM, _root, MPI_COMM_WORLD);
  if (_root==_rank) {
    cout << "Total Chern Number is: " << _total_chern << endl;
  }
  if (_root==_rank){
    TotalEig = new double [_NKX2*_SMAX];
    recvcounts = new int[_size];
    displs_r = new int[_size];
    offset = 0;
    for(int ig = 0;ig<_size;++ig){
      recvcounts[ig] = compute_count(ig,_size,_NKX2)*_SMAX; // this is the only difference between TotalEig and curvature
      displs_r[ig] = offset;
      offset += recvcounts[ig];
    }
  }
  MPI_Gatherv(localEig, stride, MPI_DOUBLE, TotalEig, recvcounts, displs_r, MPI_DOUBLE, _root, COMM_WORLD);
  if (_root==_rank) {
    ofstream spectrum_output, akx, aky;
    spectrum_output.open("spectrum.OUT");
    akx.open("AKX.OUT");
    aky.open("AKY.OUT");
    assert(akx.is_open());
    assert(aky.is_open());
    assert(spectrum_output.is_open());
    for(int nky = 0;nky <_NKX;++nky){
      for(int nkx = 0;nkx<_NKX;++nkx){
	akx << _gauss_k[nkx] << '\t';
	aky << _gauss_k[nky] << '\t';
	for (int q = 0; q<_SMAX;++q){
	  spectrum_output << TotalEig[(nkx+nky*_NKX)*_SMAX+q] << '\t';
	}
	spectrum_output << endl;
      }
      akx << endl;
      aky << endl;
      spectrum_output << endl;
    }
    akx.close();
    aky.close();
    spectrum_output.close();
  }
  if (string(_chernSolver)=="curvature"){
    if (_root==_rank) {
      curvature = new double [_NKX2];
      offset = 0;
      for(int ig = 0;ig<_size;++ig){
	recvcounts[ig] = compute_count(ig,_size,_NKX2); // See difference above
	displs_r[ig] = offset;
	offset += recvcounts[ig];
      }
    }
    MPI_Gatherv(curvature_rank,recvcount,MPI_DOUBLE,curvature,recvcounts,displs_r,MPI_DOUBLE,_root,COMM_WORLD);
    
    if (_root==_rank) {
      ofstream curv_output ;
      curv_output.open("curvature.OUT");
      assert(curv_output.is_open());
      for(int nky = 0;nky <_NKX;++nky){
	for(int nkx = 0;nkx<_NKX;++nkx){
	  curv_output << curvature[nkx+nky*_NKX] << '\t';
	}
	curv_output << endl;
      }
      curv_output.close();
    }
  } else if (string(_chernSolver)=="connection")
    if (_root==_rank) cout << "In the Berry connection approach, no curvature mesh function is recorded." << endl;
  
}

void cBdG_Edge::file_input(){
  
  
}

void cBdG_Edge::file_output(){


}
void cBdG_Edge::construction(){
  
  
}
void cBdG_Edge::compute(){


}
