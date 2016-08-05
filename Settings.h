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
int CluSizeCutY = 5;

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
        {6,TString("DownX")}, {7,TString("DownX")}, {0,TString("DownX")},
        {1,TString("DownY")}, {2,TString("DownY")}, {3,TString("DownY")},
        {8,TString("UpX")}, {9,TString("UpX")}, {4,TString("UpX")},
        {5,TString("UpY")}, {10,TString("UpY")}, {11,TString("UpY")},
        {12,TString("RefX")}, {13,TString("RefX")},
        {15,TString("RefY")}, {14,TString("RefY")}
    };


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
    {6,TString("DownX")}, {7,TString("DownX")}, {4,TString("DownX")},
    {5,TString("DownY")}, {2,TString("DownY")}, {3,TString("DownY")},
    {8,TString("UpX")}, {9,TString("UpX")}, {0,TString("UpX")},
    {1,TString("UpY")}, {10,TString("UpY")}, {11,TString("UpY")},
    {12,TString("RefX")}, {13,TString("RefX")},
    {15,TString("RefY")}, {14,TString("RefY")}
};

// Define Layer IDs, you can define ids as you like
// For the alignment step : By default
// layerID 0 should belong to the fixed X axis and
// layerID 1 should belong to the fixed Y axis
std::map<TString, unsigned short> layerIDmap={
    {TString("DownX"),0},{TString("DownY"),1},
    {TString("UpX"),2},{TString("UpY"),3},
    {TString("RefX"),4},{TString("RefY"),5}
};

std::map<unsigned short,TString> IDlayermap={
    {0,TString("DownX")},{1,TString("DownY")},
    {2,TString("UpX")},{3,TString("UpY")},
    {4,TString("RefX")},{5,TString("RefY")}
};

// Ref, Up, Down
int xlayers[3] ={4,2,0};
int ylayers[3] = {5,3,1};
int zpos[3] = {0,138,182};

std::map<unsigned short,double> layerZposition={
    {0,182.0},{1,182.0},
    {2,138.0},{3,138.0},
    {4,0.0},{5,0.0}
};

std::map<unsigned short,double> layerResolution={
    {0,2000.0},{1,2000.0},
    {2,2000.0},{3,2000.0},
    {4,2000.0},{5,2000.0}
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

