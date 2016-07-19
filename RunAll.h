/*
 * Created by Eda Yildirim
 *  (2016 JGU)
 *
 *  email:eda.yildirim@cern.ch
 */


#ifndef RunAll_h
#define RunAll_h

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "TString.h"
#include "TColor.h"

using namespace std;

TString runlistfile = TString("./MM_DESYTB_May_2016.txt");

Color_t color[9] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange-3, kPink-9, kViolet+1, kCyan+2};


void readRunList(TString filename);

struct Meas {
    int RunNum;
    double Angle1;
    double Angle2;
    double BeamEnergy;
    
    double eff_Down[2];
    double eff_Up[2];
    double eff_Ref[2];
    
    double res_Down[2];
    double res_Down_err[2];
    double res_Up[2];
    double res_Up_err[2];
    
    double ang_res;
    double ang_res_err;

    Meas(){};
    ~Meas(){};

    void set_eff_DownX(double avalue) {eff_Down[0] = avalue;}
    void set_eff_DownY(double avalue) {eff_Down[1] = avalue;}
    void set_eff_UpX(double avalue) {eff_Up[0] = avalue;}
    void set_eff_UpY(double avalue) {eff_Up[1] = avalue;}
    void set_eff_RefX(double avalue) {eff_Ref[0] = avalue;}
    void set_eff_RefY(double avalue) {eff_Ref[1] = avalue;}

    double get_eff_DownX() {return eff_Down[0];}
    double get_eff_DownY() {return eff_Down[1];}
    double get_eff_UpX() {return eff_Up[0];}
    double get_eff_UpY() {return eff_Up[1];}
    double get_eff_RefX() {return eff_Ref[0];}
    double get_eff_RefY() {return eff_Ref[1];}

    void set_res_DownX(double avalue) {res_Down[0] = avalue;}
    void set_res_DownY(double avalue) {res_Down[1] = avalue;}
    void set_res_DownX_err(double avalue) {res_Down_err[0] = avalue;}
    void set_res_DownY_err(double avalue) {res_Down_err[1] = avalue;}
    void set_res_UpX(double avalue) {res_Up[0] = avalue;}
    void set_res_UpY(double avalue) {res_Up[1] = avalue;}
    void set_res_UpX_err(double avalue) {res_Up_err[0] = avalue;}
    void set_res_UpY_err(double avalue) {res_Up_err[1] = avalue;}
    void set_ang_res(double avalue) {ang_res = avalue;}
    void set_ang_res_err(double avalue) {ang_res_err = avalue;}

    double get_res_DownX() {return res_Down[0];}
    double get_res_DownY() {return res_Down[1];}
    double get_res_DownX_err() {return res_Down_err[0];}
    double get_res_DownY_err() {return res_Down_err[1];}
    double get_res_UpX() {return res_Up[0];}
    double get_res_UpY() {return res_Up[1];}
    double get_res_UpX_err() {return res_Up_err[0];}
    double get_res_UpY_err() {return res_Up_err[1];}
    double get_ang_res() {return ang_res;}
    double get_ang_res_err() {return ang_res_err;}

    void set_RunNum (int avalue) { RunNum = avalue;}
    void set_Angle1(double avalue) {Angle1 = avalue;}
    void set_Angle2(double avalue) {Angle2 = avalue;}
    void set_BeamEnergy(double avalue) { BeamEnergy= avalue;}
    
    int get_RunNum() {return RunNum;}
    double get_Angle1() {return Angle1;}
    double get_Angle2() {return Angle2;}
    double get_BeamEnergy() {return BeamEnergy;}


};

struct RunInfo{
    int PedRun;
    int DataRun;
    double BeamEnergy;
    double Distance1;
    double Angle1;
    double Distance2;
    double Angle2;
    string config;
    string map;
    string Comment;

    void print(){
        cout<<"DataRun "<< DataRun <<" PedRun "<<PedRun<<" Angle1 "<<Angle1<< " Angle2 "<<Angle2<<" config "<<config<< " BeamEnergy "<< BeamEnergy <<endl;
    }
};

void runHitMaker();
void runAlignment();
void runAnalyze();


void readMeasEff();
void readMeasRes();
void makePlots();
void makeMultiGraps();

void createOutputFile(TFile *fout);

#endif
