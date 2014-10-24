/*
 * main.cpp
 *
 *  Created on: Oct 22, 2014
 *      Author: Lin Dong
 */

#include "stdcpp.h"
#include "dist.h"
#include "floquet.h"
#include "bulk.h"
#include "bdg.h"
#include "edge.h"
#include "prx.h"
#define root 0
#define char_length	100
int main(int argc, char** argv){
  Init(argc,argv);
  int rank, size;
  rank = COMM_WORLD.Get_rank();
  size = COMM_WORLD.Get_size();
  /******************************I/O for main.cpp************************************/
  if (rank == root) {
    cout << "======================================================================\n"
	 << "The purpose of this program is to study\n"
	 << "the Floquet Topological property of a periodically driven Hamiltonian.\n"
	 << "Motivated, proposed, designed, implemented and researched \n"
	 << "by Lin Dong at Rice University. \n"
	 << "at " << __TIME__ << ", on " << __DATE__ << endl
	 << "MPI is initialized and program starts from " << __FILE__ << endl
	 << "======================================================================\n";  
  }
  char dummy_input_name[char_length], hamiltonian[char_length], topology[char_length];
  FILE *input_open;
  input_open = fopen("input.txt","r");
  assert(input_open != NULL);
  /* example file of input.txt with comment:
     Hamiltonian             bdg  || prx
     topology                bulk || edge
  */
  fscanf(input_open,"%s %s", dummy_input_name, hamiltonian);
  if (rank == root) {
    cout << dummy_input_name << "=" << hamiltonian  <<endl;
    if (string(hamiltonian) != "bdg" && string(hamiltonian) != "prx" ) {
      cerr << "Hamiltonian is not supported. Please check input.txt and readme." << endl;
      exit(1);
    }
  }
  fscanf(input_open,"%s %s", dummy_input_name, topology);
  if (rank == root) {
    cout << dummy_input_name << "=" << topology << endl;
    if (string(topology) != "bulk" && string(topology) != "edge") {
      cerr << "Topological property is not supported. Please check input.txt and readme." << endl;
      exit(1);
    }
  }
  fclose(input_open);
  /************************************************************************************/
  //cDistribute dist(rank,size,root);
  cout.precision(16);
  switch (string(hamiltonian)=="bdg") {
	case 1:
		cBdG* bdg;
		switch (string(topology)=="bulk") {
			case 1:
				bdg = new cBdG_Bulk(rank,size,root);
				break;
			case 0:
				bdg = new cBdG_Edge(rank,size,root);
				break;
			default:
				break;
		}
		delete bdg;
		break;
	case 0:
		cPRX* prx;
		switch (string(topology)=="bulk") {
			case 1:
				prx = new cPRX_Bulk(rank,size,root);
				break;
			case 0:
				prx = new cPRX_Edge(rank,size,root);
				break;
			default:
				break;
		}
		delete prx;
		break;
	default:
		break;
  }
  Finalize();
  return 0;
}

