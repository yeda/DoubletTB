/*
 * Created by Eda Yildirim
 *  (2016 JGU)
 *
 *  email:eda.yildirim@cern.ch
 */
#ifndef Settings_h
#define Settings_h

#include <iostream>
#include <map>
#include <string>
#include <stdlib.h>
#include <sstream>


#include "TString.h"

using namespace std;

TString output_resolution_txtfile= TString("./results/resolution.txt");
TString output_efficiency_txtfile= TString("./results/efficiency.txt");


// should end with /
TString inputPath = TString("./input/");
TString outputPath = TString("./output/");

TString inputRawTreeName = TString("raw");
TString inputDataTreeName = TString("data");

const double pitchsize = 0.25; // in mm

const int NLayer=6;

//////////// HitMaker /////////
TString hitAmplitudeCut_file= TString("./HitAmplitudeCuts.txt");
int SignalCut = 1;




/*
std::map<TString, double> HitAmplitudeCut={
    {TString("DownX"),1},{TString("DownY"),1},
    {TString("UpX"),1},{TString("UpY"),1},
    {TString("RefX"),1},{TString("RefY"),1}
};
*/

// local coordinates
int CluSizeCutX = 1;
int CluSizeCutY = 1;

const unsigned int MaxEventNum = 10000000;
//const unsigned int MaxEventNum = 1000000;
const int MAXEVENTTOPLOT = 10;


//////////// Alignment /////////
bool PlotHistosAfterAlignment =  true;

// for desyTB.config

// TMM1 -> Dublett down
// mmdaq.Chamber.TMM1.Strips.X.Chips: APV6, APV7, APV0
// mmdaq.Chamber.TMM1.Strips.Y.Chips: APV1, APV2, APV3
// TMM2 -> Dublett up
//mmdaq.Chamber.TMM2.Strips.X.Chips: APV8, APV9, APV4
//mmdaq.Chamber.TMM2.Strips.Y.Chips: APV5, APV10, APV11
// TMM3 -> Reference
//mmdaq.Chamber.TMM3.Strips.X.Chips: APV12, APV13
//mmdaq.Chamber.TMM3.Strips.Y.Chips: APV15, APV14

// Note that "X" and "Y" is reserved to define axis

std::map<int,TString> apvIDmap2 =
    {
        {6,TString("B-L2")}, {7,TString("B-L2")}, {0,TString("B-L2")},
        {1,TString("B-L1")}, {2,TString("B-L1")}, {3,TString("B-L1")},
        {8,TString("A-L2")}, {9,TString("A-L2")}, {4,TString("A-L2")},
        {5,TString("A-L1")}, {10,TString("A-L1")}, {11,TString("A-L1")},
        {12,TString("R-L2")}, {13,TString("R-L2")},
        {15,TString("R-L1")}, {14,TString("R-L1")}
    };

/*
 DownX->B-L2, DownY->B-L1, UpX->A-L2, UpY->A-L1, RefX->R-L2, RefY->R-L1
*/


// for desyTBUntil11_05_2016_15_30.config

// TMM1 -> Dublett down
//mmdaq.Chamber.TMM1.Strips.X.Chips: APV6, APV7, APV4
//mmdaq.Chamber.TMM1.Strips.Y.Chips: APV5, APV2, APV3
// TMM2 -> Dublett up
//mmdaq.Chamber.TMM2.Strips.X.Chips: APV8, APV9, APV0
//mmdaq.Chamber.TMM2.Strips.Y.Chips: APV1, APV10, APV11
// TMM3 -> Reference
//mmdaq.Chamber.TMM3.Strips.X.Chips: APV12, APV13
//mmdaq.Chamber.TMM3.Strips.Y.Chips: APV15, APV14

// Note that "X" and "Y" is reserved to define axis

std::map<int,TString> apvIDmap1 =
{
    {6,TString("B-L2")}, {7,TString("B-L2")}, {4,TString("B-L2")},
    {5,TString("B-L1")}, {2,TString("B-L1")}, {3,TString("B-L1")},
    {8,TString("A-L2")}, {9,TString("A-L2")}, {0,TString("A-L2")},
    {1,TString("A-L1")}, {10,TString("A-L1")}, {11,TString("A-L1")},
    {12,TString("R-L2")}, {13,TString("R-L2")},
    {15,TString("R-L1")}, {14,TString("R-L1")}
};

// Define Layer IDs, you can define ids as you like
// For the alignment step : By default
// layerID 0 should belong to the fixed X axis and
// layerID 1 should belong to the fixed Y axis
std::map<TString, unsigned short> layerIDmap={
    {TString("R-L2"),0},{TString("R-L1"),1},
    {TString("A-L1"),2},{TString("A-L2"),3},
    {TString("B-L2"),4},{TString("B-L1"),5}
};
/*
 DownX->B-L2, DownY->B-L1, UpX->A-L2, UpY->A-L1, RefX->R-L2, RefY->R-L1
 */

std::map<unsigned short,TString> IDlayermap={
    {0,TString("R-L2")},{1,TString("R-L1")},
    {2,TString("A-L1")},{3,TString("A-L2")},
    {4,TString("B-L2")},{5,TString("B-L1")}
};

// Ref, Up, Down --> R, A, B
int xlayers[3] ={0,2,4};
int ylayers[3] = {1,3,5};
int zpos[3] = {0,138,182};

std::map<unsigned short,double> layerZposition={
    {0,0.0},{1,0.0},
    {2,138.0},{3,138.0},
    {4,182.0},{5,182.0}
};

std::map<unsigned short,double> layerResolution={
    {0,10.0},{1,10.0},
    {2,10.0},{3,10.0},
    {4,10.0},{5,10.0}
};


vector<string> splitstring(string s, char delim) {
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

#endif

