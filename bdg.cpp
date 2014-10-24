/*
 * bdg.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: ld7
 */
#include "bdg.h"

void cBdG_Bulk::file_input(){
char dummyname[100];
double dummyvalue;
int intdummyvalue;
  FILE *input_open;
  input_open = fopen("input_bdg_bulk.txt","r");
  assert(input_open != NULL);
  if (_rank == _root)  cout << "Starting to read in parameters from file input_bdg_bulk.txt" << endl;
  fscanf(input_open,"%s %d", dummyname, &intdummyvalue);
  _PMAX = intdummyvalue;	if (_rank == _root) cout << dummyname << "=" << _PMAX << endl;

  fclose(input_open);
}

void cBdG_Bulk::construction(){

	cout << "construction has started\n";
}


void cBdG_Edge::file_input(){
char dummyname[100];
double dummyvalue;
int intdummyvalue;
  FILE *input_open;
  input_open = fopen("input_bdg_bulk.txt","r");
  assert(input_open != NULL);
  if (_rank == _root)  cout << "Starting to read in parameters from file input_bdg_bulk.txt" << endl;
  fscanf(input_open,"%s %d", dummyname, &intdummyvalue);
  _PMAX = intdummyvalue;	if (_rank == _root) cout << dummyname << "=" << _PMAX << endl;

  fclose(input_open);
}

void cBdG_Edge::construction(){

	cout << "construction has started\n";
}
