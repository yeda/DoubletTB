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
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TAxis.h"

using namespace std;

TString runlistfile = TString("./MM_DESYTB_May_2016.txt");

Color_t color[9] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange-3, kPink-9, kViolet+1, kCyan+2};

    // Doublet B (B) correction --> doublet b res = m*(real res)+b 
    //	m = 2.54091, b = 0.00527 mm
    // Doublet A (A) correction --> doublet a res = m*(real res)+b 
    //	m = 1.31425, b = 0.00640 mm

double m_B = 2.54091;
double b_B = 0.00527;
double m_A = 1.31425;
double b_A = 0.00640;

void readRunList(TString filename);

struct Meas {
    int RunNum;
    double Angle1;
    double Angle2;
    double BeamEnergy;
    
    double eff_B[2];
    double eff_A[2];
    double eff_R[2];
    
    double eff_B_err[2];
    double eff_A_err[2];
    double eff_R_err[2];
    
    double res_R[2];
    double res_R_err[2];
    double res_B[2];
    double res_B_err[2];
    double res_A[2];
    double res_A_err[2];
    
    double ang_res;
    double ang_res_err;
    
    Meas(){};
    ~Meas(){};
    void set_eff_BL1(double avalue) {eff_B[0] = avalue;}
    void set_eff_BL2(double avalue) {eff_B[1] = avalue;}
    void set_eff_AL1(double avalue) {eff_A[0] = avalue;}
    void set_eff_AL2(double avalue) {eff_A[1] = avalue;}
    void set_eff_RL1(double avalue) {eff_R[0] = avalue;}
    void set_eff_RL2(double avalue) {eff_R[1] = avalue;}
    
    double get_eff_BL1() {return eff_B[0];}
    double get_eff_BL2() {return eff_B[1];}
    double get_eff_AL1() {return eff_A[0];}
    double get_eff_AL2() {return eff_A[1];}
    double get_eff_RL1() {return eff_R[0];}
    double get_eff_RL2() {return eff_R[1];}
    
    void set_eff_BL1_err(double avalue) {eff_B_err[0] = avalue;}
    void set_eff_BL2_err(double avalue) {eff_B_err[1] = avalue;}
    void set_eff_AL1_err(double avalue) {eff_A_err[0] = avalue;}
    void set_eff_AL2_err(double avalue) {eff_A_err[1] = avalue;}
    void set_eff_RL1_err(double avalue) {eff_R_err[0] = avalue;}
    void set_eff_RL2_err(double avalue) {eff_R_err[1] = avalue;}
    
    double get_eff_BL1_err() {return eff_B_err[0];}
    double get_eff_BL2_err() {return eff_B_err[1];}
    double get_eff_AL1_err() {return eff_A_err[0];}
    double get_eff_AL2_err() {return eff_A_err[1];}
    double get_eff_RL1_err() {return eff_R_err[0];}
    double get_eff_RL2_err() {return eff_R_err[1];}
    
    void set_res_RL1(double avalue) {res_R[0] = avalue;}
    void set_res_RL2(double avalue) {res_R[1] = avalue;}
    void set_res_RL1_err(double avalue) {res_R_err[0] = avalue;}
    void set_res_RL2_err(double avalue) {res_R_err[1] = avalue;}
    void set_res_BL1(double avalue) {res_B[0] = avalue;}
    void set_res_BL2(double avalue) {res_B[1] = avalue;}
    void set_res_BL1_err(double avalue) {res_B_err[0] = avalue;}
    void set_res_BL2_err(double avalue) {res_B_err[1] = avalue;}
    void set_res_AL1(double avalue) {res_A[0] = avalue;}
    void set_res_AL2(double avalue) {res_A[1] = avalue;}
    void set_res_AL1_err(double avalue) {res_A_err[0] = avalue;}
    void set_res_AL2_err(double avalue) {res_A_err[1] = avalue;}
    void set_ang_res(double avalue) {ang_res = avalue;}
    void set_ang_res_err(double avalue) {ang_res_err = avalue;}
    
    double get_res_RL1() {return res_R[0];}
    double get_res_RL2() {return res_R[1];}
    double get_res_RL1_err() {return res_R_err[0];}
    double get_res_RL2_err() {return res_R_err[1];}
    double get_res_BL1() {return res_B[0];}
    double get_res_BL2() {return res_B[1];}
    double get_res_BL1_err() {return res_B_err[0];}
    double get_res_BL2_err() {return res_B_err[1];}
    double get_res_AL1() {return res_A[0];}
    double get_res_AL2() {return res_A[1];}
    double get_res_AL1_err() {return res_A_err[0];}
    double get_res_AL2_err() {return res_A_err[1];}

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


void formatCanvas1D(TCanvas *cc){
    cc->SetRightMargin(0.125);
    cc->SetLeftMargin(0.125);
    cc->SetBottomMargin(0.13);
    cc->SetTopMargin(0.07);
}


void formatMultiGraph(TMultiGraph *gr){
    gr->GetXaxis()->SetLabelSize(0.045);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetTitleOffset(1.2);
    
    gr->GetYaxis()->SetLabelSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleOffset(1.2);
    
    //gr->SetTitle("");
}

#endif
