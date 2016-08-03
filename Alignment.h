/*
 * Created by Eda Yildirim
 *  (2016 JGU)
 *
 *  email:eda.yildirim@cern.ch
 */
#ifndef Alignment_h
#define Alignment_h

#include <iostream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObject.h"

using namespace std;

void fillHistos();
void alignLayers();
void readAlignmentParameters(double alignmentpar[]);
void fillHistosAfterAlignment(double alignmentpar[]);

bool isGoodEvent();

void setBranches(TTree* input_tree);

void createHistos();

void createOutputFile(TFile *fout);
double getExpectedHit(double x1, double y1, double x2, double y2, double at_y);


#endif



