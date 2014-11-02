/*
 * edge.cpp
 *
 *  Created on: Nov 1, 2014
 *      Author: ld7
 */
#include "bdg.h"
void cBdG_Edge::file_input(){
  for (int ig = 0; ig < _size; ++ig) {
    if (ig ==_rank){
      char dummyname[100];
      double dummyvalue;
      int intdummyvalue;
      FILE *input_bdg;
      input_bdg = fopen("input_bdg_edge.txt","r");
      assert(input_bdg != NULL);
      if (ig == _root)  cout << "Starting to read in parameters from file input_bdg_edge.txt" << endl;
      fscanf(input_bdg,"%s %d", dummyname, &intdummyvalue);
      _PMAX = intdummyvalue;	if (ig == _root) cout << dummyname << "=" << _PMAX << endl;
      fscanf(input_bdg,"%s %d", dummyname, &intdummyvalue);
      _NMAX = intdummyvalue;	if (ig == _root) cout << dummyname << "=" << _NMAX << endl;
      fscanf(input_bdg,"%s %d", dummyname, &intdummyvalue);
      _NKX = intdummyvalue;	if (ig == _root) cout << dummyname << "=" << _NKX << endl;
      fscanf(input_bdg,"%s %lf", dummyname, &dummyvalue);
      _h = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _h << endl;
      fscanf(input_bdg,"%s %lf", dummyname, &dummyvalue);
      _mu = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _mu << endl;
      fscanf(input_bdg,"%s %lf", dummyname, &dummyvalue);
      _T = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _T << endl;
      fscanf(input_bdg,"%s %lf", dummyname, &dummyvalue);
      _v = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _v << endl;
      fscanf(input_bdg,"%s %lf", dummyname, &dummyvalue);
      _kmax = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _kmax << endl;
      fscanf(input_bdg,"%s %lf", dummyname, &dummyvalue);
      _L = dummyvalue;	if (ig == _root) cout << dummyname << "=" << _L << endl;
      fclose(input_bdg);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
}


void cBdG_Edge::construction(){
  if (_rank == _root) cout << "Construction has started\n";
  _ibdg = 4; // --> It has to be 4 here. (at least currently.)
  _pblock = (2*_PMAX+1);
  _pblock4 = _ibdg*_pblock;
  _SMAX = _pblock4*_NMAX;
  _bdg_E.resize(_SMAX);
  _bdg_V.resize(_SMAX,_SMAX);
  _bdg_H.resize(_SMAX,_SMAX);
  update(-1); // null construction that is only done once for arbitrary momentum state.
  if (_rank == _root) cout << "Construction is completed." << endl;
  
}

void cBdG_Edge::update(int nkx){
  int p,q,m,n;
  if (nkx == -1) {
    _bdg_H.setZero(); // This is done only once.
    // The off-diagonal coupling introduced from time-dependent order parameter should be computed only here.
    double dt = 0.0005;int lenT = int(_T/dt);
    VectorXcd Delta_t(lenT+1000);Delta_t.setZero();
    FILE *R_Delta, *I_Delta;
    R_Delta = fopen("Rdata_2109.dat","r");
    I_Delta = fopen("Idata_2109.dat","r");
    assert(R_Delta != NULL);assert(I_Delta != NULL);
    int count = 0;
    double reD, imD;
    while (fscanf(R_Delta, "%lf", &reD) != EOF && fscanf(I_Delta, "%lf", &imD) != EOF ){
      Delta_t(count) = complex<double>(reD,imD);
      count++;
    }
    fclose(R_Delta);fclose(I_Delta);
    complex<double> Gamma1, Gamma2, temp;
    double Lambda;
    for (int ip = 0; ip < _pblock; ++ip) {
      p = ip - _PMAX;
      for (int iq = 0; iq<=ip;++iq){
		q = iq - _PMAX;
		Gamma1 = complex<double> (0.0,0.0);Gamma2 = complex<double> (0.0,0.0);
		for (int it = 0; it < count; ++it) {
			temp = Delta_t[it];
			Gamma1 += (temp)		*exp(_myI*2*M_PI*(q-p)/_T*it*dt)/_T*dt;
			Gamma2 += conj(temp)*exp(_myI*2*M_PI*(q-p)/_T*it*dt)/_T*dt;
		}
		//cout << "gamma1=" << Gamma1 << "gamma2="<< Gamma2 << endl;
		for (int im = 0; im < _NMAX; ++im) {
			_bdg_H(ip*_ibdg*_NMAX+im*_ibdg	,iq*_ibdg*_NMAX+im*_ibdg+3) = -Gamma1;
			_bdg_H(ip*_ibdg*_NMAX+im*_ibdg+1,iq*_ibdg*_NMAX+im*_ibdg+2) =  Gamma1;
			_bdg_H(ip*_ibdg*_NMAX+im*_ibdg+2,iq*_ibdg*_NMAX+im*_ibdg+1) =  Gamma2;
			_bdg_H(ip*_ibdg*_NMAX+im*_ibdg+3,iq*_ibdg*_NMAX+im*_ibdg)	= -Gamma2;
		}
      }
      for (int im = 0; im < _NMAX; ++im) {
    	  m = im+1;
    	  for (int in = 0; in < im; ++in) {
    		  n = in+1;
//    		  if (in!=im){
    		  Lambda = 2.0*_v*m*n*(1-pow(-1.0,m+n))/_L/(m*m-n*n);
    		  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg	,ip*_ibdg*_NMAX+in*_ibdg+1) = complex<double>(-Lambda,0.0);
    		  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg+1	,ip*_ibdg*_NMAX+in*_ibdg)	=  complex<double>(Lambda,0.0);
    		  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg+2	,ip*_ibdg*_NMAX+in*_ibdg+3) =  complex<double>(Lambda,0.0);
    		  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg+3	,ip*_ibdg*_NMAX+in*_ibdg+2)	= complex<double>(-Lambda,0.0);
//    		  }
    	  }
      }
    }
  } else {
    double kx = -_kmax +2.0*_kmax/(_NKX-1)*nkx, xi = kx*kx-_mu;
    for (int ip = 0; ip < _pblock; ++ip) {
      p = ip - _PMAX;
      for (int im = 0; im < _NMAX; ++im) {
    	  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg,ip*_ibdg*_NMAX+im*_ibdg) =
    			  complex<double> (xi+pow(m*M_PI/_L,2.0)+_h+2*p*M_PI/_T,0.0);

    	  //_bdg_H(ip*_ibdg*_NMAX+im*_ibdg,ip*_ibdg*_NMAX+im*_ibdg+1) =
    	  //    			  complex<double> (_v*kx,0.0);

    	  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg+1,ip*_ibdg*_NMAX+im*_ibdg) =
    			  complex<double> (_v*kx,0.0);
    	  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg+1,ip*_ibdg*_NMAX+im*_ibdg+1) =
    			  complex<double> (xi+pow(m*M_PI/_L,2.0)-_h+2*p*M_PI/_T,0.0);
    	  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg+2,ip*_ibdg*_NMAX+im*_ibdg+2) =
    			  complex<double> (-(xi+pow(m*M_PI/_L,2.0)+_h)+2*p*M_PI/_T,0.0);

    	  //_bdg_H(ip*_ibdg*_NMAX+im*_ibdg+2,ip*_ibdg*_NMAX+im*_ibdg+3) =
    	  //    			  complex<double> (_v*kx,0.0);

    	  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg+3,ip*_ibdg*_NMAX+im*_ibdg+2) =
    			  complex<double> (_v*kx,0.0);
    	  _bdg_H(ip*_ibdg*_NMAX+im*_ibdg+3,ip*_ibdg*_NMAX+im*_ibdg+3) =
    			  complex<double> (-(xi+pow(m*M_PI/_L,2.0)-_h)+2*p*M_PI/_T,0.0);
      }
	}
  }
}

void cBdG_Edge::compute(){
	distribution(_NKX);
	stride = _SMAX*recvcount;
	localEig = new double[stride];
	SelfAdjointEigenSolver<MatrixXcd> ces;
	for(int ig = 0; ig< _size; ++ig) {
		if (ig == _rank){
		  cout << "rank " << ig << " has started"<< endl;
		  for (int i=0; i<recvcount; ++i) { // recvcount is the number of tasks given to each rank.
			clock_t start = clock();
			update(recvbuf[i]);// distributed momentum value
			ces.compute(_bdg_H,false); // only eigenvalues are computed
			clock_t end = clock();
			if (_rank==_root) cout << "task " << recvbuf[i] <<"out of " << recvcount << "used " << double (end-start)/ (double) CLOCKS_PER_SEC  << endl;
			_bdg_E = ces.eigenvalues();
			for(int j = 0; j < _SMAX; ++j){
			  localEig[j+i*_SMAX]=_bdg_E[j];
			}
			cout << "rank " << _rank << " has finished "<< recvcount << " tasks " << endl;
		  }
		}
	}
}

void cBdG_Edge::file_output(){
  if (_root==_rank){
    TotalEig = new double [_NKX*_SMAX];
    recvcounts = new int[_size];
    displs_r = new int[_size];
    offset = 0;
    for(int ig = 0;ig<_size;++ig){
      recvcounts[ig] = compute_count(ig,_size,_NKX)*_SMAX; // this is the only difference between TotalEig and curvature
      displs_r[ig] = offset;
      offset += recvcounts[ig];
    }
  }
  MPI_Gatherv(localEig, stride, MPI_DOUBLE, TotalEig, recvcounts, displs_r, MPI_DOUBLE, _root, COMM_WORLD);
  if (_root==_rank) {
	//  update(5);
	  ofstream spectrum_output, akx;
    //ofstream bdgR,bdgI;
    spectrum_output.open("spectrum.OUT");
    //bdgR.open("bdgR.OUT");bdgI.open("bdgI.OUT");assert(bdgR.is_open());assert(bdgI.is_open());
    //bdgR << _bdg_H.real() << endl;
    //bdgI << _bdg_H.imag() << endl;
    akx.open("AKX.OUT");
    assert(akx.is_open());
    assert(spectrum_output.is_open());
    for(int nkx = 0;nkx<_NKX;++nkx){
      akx << -_kmax +2.0*_kmax/(_NKX-1)*nkx<<endl;
      for (int q = 0; q<_SMAX;++q){
	spectrum_output << TotalEig[nkx*_SMAX+q] << '\t';
      }
      spectrum_output << endl;
    }
    akx << endl;
    //bdgR.close();bdgI.close();
    akx.close();spectrum_output.close();
  }
}

