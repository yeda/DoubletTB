/*
 * Created by Eda Yildirim
 *  (2016 JGU)
 *
 *  email:eda.yildirim@cern.ch
 */
#ifndef Analyze_h
#define Analyze_h

#include <iostream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObject.h"

using namespace std;

void createHistos();
void fillHistos();

void readAlignmentParameters(TString alignmenttxt, double alignmentpar[]);
void calculateEfficiency();
void fillHitMap2D();

void printResolution(TString runnum);
void printEfficiency(TString runnum);

bool isGoodEvent();
void setBranches(TTree* input_tree);
void createOutputFile(TFile *fout);

double getAngleABC(double A[] , double B[], double C[]);
void findRefLayers(int layerid, int* arr);
double getExpectedHit(double x1, double y1, double x2, double y2, double at_y);
#endif



