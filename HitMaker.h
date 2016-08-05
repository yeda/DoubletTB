/*
 * Created by Eda Yildirim
 *  (2016 JGU)
 *
 *  email:eda.yildirim@cern.ch
 */
#ifndef HitMaker_h
#define HitMaker_h

#include <iostream>
#include <vector>
#include <map>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObject.h"

using namespace std;

struct Point {
    double x;
    double y;
    double z;
};

struct Cluster {
    unsigned short _layerID;
    unsigned short _cluster_size;
    float _hit_amplitude;
    float _hit_position;
    
};



// Declaration of leaf types
// in the root file saved by mmdaq

void setBranches(TTree* rawtree, TTree* datatree);

void createHitMapHistos();

void createOutputFile(TFile *fout);

void plotEvent(TString prefix);
void processEvent(vector<Cluster> clusters);

void createOtherHistos();
bool isGoodEvent(vector<Cluster> clus);
bool isNoHitOnLayer(vector<Cluster> clus);
void readHitAmplitudeCuts(TString runnum);


#endif

/*
 
 UInt_t                  apv_evt;
 Int_t                   time_s;
 Int_t                   time_us;
 vector<UInt_t>    *apv_fecNo;
 vector<UInt_t>    *apv_id;
 vector<UInt_t>    *apv_ch;
 vector<string>          *mm_id;
 vector<UInt_t>    *mm_readout;
 vector<UInt_t>    *mm_strip;
 vector<vector<short> >  *apv_q;
 UInt_t                  apv_presamples;
 vector<short>           *apv_qmax;
 vector<short>           *apv_tbqmax;
 
 */




