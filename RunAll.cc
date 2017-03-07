/*
 * Created by Eda Yildirim
 *  (2016 JGU)
 *
 *  email:eda.yildirim@cern.ch
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <map>
#include <utility>
#include <algorithm>

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TObject.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLatex.h"

#include "RunAll.h"
#include "Settings.h"


using namespace std;

const double NOTSET = 100000;

map<int,bool> is_tiltX ={ {71, false},{74, false},{76, false},{79, false},{84, false},{86, false},{88, false},{90, false},{92, false},{139, true}, {137, true}, {102, true}, {104, true}, {106, true}, {110, true}, {112, true}, {118, true}, {124, true}, {127, true}};


TString output_root_filename = TString("./results/measurementResults.root");

TString eff_histoname= TString("eff_");
TString res_histoname= TString("res_");
TString angres_histoname= TString("angres");

vector<RunInfo> runlist;
map<int, Meas* > measurements;

map<TString,TObject*> rootobjects;

int main(int argc, char *argv[]){
    
    readRunList(runlistfile);
    
    //runHitMaker();
  //runAlignment();
    //system("./alignment 134");
     
    //system("rm ./results/efficiency.txt"); system("rm ./results/resolution.txt");
    //runAnalyze();
    readMeasEff();
    readMeasRes();
    
    TFile *fout = new TFile(output_root_filename.Data(),"RECREATE");
    
    makePlots();
    makeMultiGraps();
    createOutputFile(fout);
    
    return 0;
}



void makePlots(){
    TString histname, title;
    TGraphErrors* gre;
    int runnum;
    
    TString tiltname;
    TString basehistname;
    // comb_result[tiltx/y][B/A/R][L1/L2][angle,eff,eff_err,spa_res,spa_res_err][meas]
    double comb_result[2][3][2][5][30]={0};
    // comb_ang_res[tiltx/y][angle,angres,angreserr][meas]
    double comb_ang_res[2][3][30]={0};
    double angle_err[30];
    fill(angle_err, angle_err+30, 0.5);
    double zero_err[30]={0};
    
    int i_point[2] = {0};
    for (map< int, Meas* >::iterator i_meas=measurements.begin(); i_meas != measurements.end(); i_meas++) {
        runnum = i_meas->first;
        Meas* ameas = i_meas->second;
        
        if (is_tiltX.find(runnum) == is_tiltX.end()) continue;
        
        if(is_tiltX[runnum]){
            // comb_result[tiltx/y][B/A/R][L1/L2][angle,eff,eff_err,spa_res,spa_res_err][meas]
            // comb_ang_res[tiltx/y][angle,angres,angreserr][meas]
            
            // tiltX, B, L1
            comb_result[0][0][0][0][i_point[0]]=ameas->get_Angle2();
            comb_result[0][0][0][1][i_point[0]]=ameas->get_eff_BL1();
            comb_result[0][0][0][2][i_point[0]]=ameas->get_eff_BL1_err();
            comb_result[0][0][0][3][i_point[0]]=ameas->get_res_BL1();
            comb_result[0][0][0][4][i_point[0]]=ameas->get_res_BL1_err();

            // tiltX, B, L2
            comb_result[0][0][1][0][i_point[0]]=ameas->get_Angle2();
            comb_result[0][0][1][1][i_point[0]]=ameas->get_eff_BL2();
            comb_result[0][0][1][2][i_point[0]]=ameas->get_eff_BL2_err();
            comb_result[0][0][1][3][i_point[0]]=ameas->get_res_BL2();
            comb_result[0][0][1][4][i_point[0]]=ameas->get_res_BL2_err();
            
            // tiltX, A, L1
            comb_result[0][1][0][0][i_point[0]]=ameas->get_Angle2();
            comb_result[0][1][0][1][i_point[0]]=ameas->get_eff_AL1();
            comb_result[0][1][0][2][i_point[0]]=ameas->get_eff_AL1_err();
            comb_result[0][1][0][3][i_point[0]]=ameas->get_res_AL1();
            comb_result[0][1][0][4][i_point[0]]=ameas->get_res_AL1_err();

            // tiltX, A, L2
            comb_result[0][1][1][0][i_point[0]]=ameas->get_Angle2();
            comb_result[0][1][1][1][i_point[0]]=ameas->get_eff_AL2();
            comb_result[0][1][1][2][i_point[0]]=ameas->get_eff_AL2_err();
            comb_result[0][1][1][3][i_point[0]]=ameas->get_res_AL2();
            comb_result[0][1][1][4][i_point[0]]=ameas->get_res_AL2_err();
            
            // tiltX, R, L1
            comb_result[0][2][0][0][i_point[0]]=ameas->get_Angle2();
            comb_result[0][2][0][1][i_point[0]]=ameas->get_eff_RL1();
            comb_result[0][2][0][2][i_point[0]]=ameas->get_eff_RL1_err();

            // tiltX, R, L2
            comb_result[0][2][1][0][i_point[0]]=ameas->get_Angle2();
            comb_result[0][2][1][1][i_point[0]]=ameas->get_eff_RL2();
            comb_result[0][2][1][2][i_point[0]]=ameas->get_eff_RL2_err();
            
            comb_ang_res[0][0][i_point[0]]=ameas->get_Angle2();
            comb_ang_res[0][1][i_point[0]]=ameas->get_ang_res();
            comb_ang_res[0][2][i_point[0]]=ameas->get_ang_res_err();
            i_point[0]++;
            
        }
        else {
            // comb_result[tiltx/y][B/A/R][L1/L2][angle,eff,eff_err,spa_res,spa_res_err][meas]
            // comb_ang_res[tiltx/y][angle,angres,angreserr][meas]
            
            // tiltY, B, L1
            comb_result[1][0][0][0][i_point[1]]=ameas->get_Angle1();
            comb_result[1][0][0][1][i_point[1]]=ameas->get_eff_BL1();
            comb_result[1][0][0][2][i_point[1]]=ameas->get_eff_BL1_err();
            comb_result[1][0][0][3][i_point[1]]=ameas->get_res_BL1();
            comb_result[1][0][0][4][i_point[1]]=ameas->get_res_BL1_err();
            
            // tiltY, B, L2
            comb_result[1][0][1][0][i_point[1]]=ameas->get_Angle1();
            comb_result[1][0][1][1][i_point[1]]=ameas->get_eff_BL2();
            comb_result[1][0][1][2][i_point[1]]=ameas->get_eff_BL2_err();
            comb_result[1][0][1][3][i_point[1]]=ameas->get_res_BL2();
            comb_result[1][0][1][4][i_point[1]]=ameas->get_res_BL2_err();
            
            // tiltY, A, L1
            comb_result[1][1][0][0][i_point[1]]=ameas->get_Angle1();
            comb_result[1][1][0][1][i_point[1]]=ameas->get_eff_AL1();
            comb_result[1][1][0][2][i_point[1]]=ameas->get_eff_AL1_err();
            comb_result[1][1][0][3][i_point[1]]=ameas->get_res_AL1();
            comb_result[1][1][0][4][i_point[1]]=ameas->get_res_AL1_err();
            
            // tiltY, A, L2
            comb_result[1][1][1][0][i_point[1]]=ameas->get_Angle1();
            comb_result[1][1][1][1][i_point[1]]=ameas->get_eff_AL2();
            comb_result[1][1][1][2][i_point[1]]=ameas->get_eff_AL2_err();
            comb_result[1][1][1][3][i_point[1]]=ameas->get_res_AL2();
            comb_result[1][1][1][4][i_point[1]]=ameas->get_res_AL2_err();
           
            // tiltY, R, L1
            comb_result[1][2][0][0][i_point[1]]=ameas->get_Angle1();
            comb_result[1][2][0][1][i_point[1]]=ameas->get_eff_RL1();
            comb_result[1][2][0][2][i_point[1]]=ameas->get_eff_RL1_err();
            
            // tiltY, R, L2
            comb_result[1][2][1][0][i_point[1]]=ameas->get_Angle1();
            comb_result[1][2][1][1][i_point[1]]=ameas->get_eff_RL2();
            comb_result[1][2][1][2][i_point[1]]=ameas->get_eff_RL2_err();
            
            comb_ang_res[1][0][i_point[1]]=ameas->get_Angle1();
            comb_ang_res[1][1][i_point[1]]=ameas->get_ang_res();
            comb_ang_res[1][2][i_point[1]]=ameas->get_ang_res_err();
            i_point[1]++;
            
        }
        
    }
    // comb_result[tiltx/y][down/up/ref][x/y][angle,eff,eff_err,spa_res,spa_res_err][meas]
    
    cout<< "i_point[0] "<<i_point[0]<<" i_point[1] "<<i_point[1]<<endl;
    
    for (int itilt=0; itilt<2; itilt++) {
        for (int idet=0; idet<3; idet++) {
            for (int ixy=0; ixy<2; ixy++) {
                for (int imeastype=1; imeastype<3; imeastype++) {
                    histname = TString("MEAS_tilt111_22233");
                    
                    if (itilt==0) histname.ReplaceAll("111", "X");
                    else histname.ReplaceAll("111", "Y");
                    
                    if(idet==0) histname.ReplaceAll("222", "B");
                    else if(idet==1) histname.ReplaceAll("222", "A");
                    else histname.ReplaceAll("222", "R");
                    
                    if (ixy==0) histname.ReplaceAll("33", "-L1");
                    else histname.ReplaceAll("33", "-L2");
                    
                    if(imeastype==0) continue;
                    else if(imeastype==1) histname.ReplaceAll("MEAS", eff_histoname);
                    else if(imeastype==2) histname.ReplaceAll("MEAS", res_histoname);
                    
                    if (imeastype==1) {
                        gre = new TGraphErrors(i_point[itilt], comb_result[itilt][idet][ixy][0], comb_result[itilt][idet][ixy][1], angle_err, comb_result[itilt][idet][ixy][2]);
                        gre->SetName(histname.Data());
                        title= histname + TString(";Tilt angle (degrees);Efficiency");
                        gre->SetTitle(title.Data());
                        gre->SetMarkerStyle(20);
                        gre->SetMarkerSize(1.6);
                        gre->Sort();
                        rootobjects.insert(pair<TString,TObject*>(histname,gre));
                    } else if (imeastype==2){
                        if (idet!=2){
                            gre = new TGraphErrors(i_point[itilt], comb_result[itilt][idet][ixy][0], comb_result[itilt][idet][ixy][3],angle_err,comb_result[itilt][idet][ixy][4]);
                            gre->SetName(histname.Data());
                            gre->SetMarkerStyle(20);
                            gre->SetMarkerSize(1.6);
                            
                            title= histname + TString(";Tilt angle (degrees);Spatial Resolution (mm)");
                            gre->Sort();
                            gre->SetTitle(title.Data());
                            rootobjects.insert(pair<TString,TObject*>(histname,gre));
                        }
                    }
                    
                }
            }
        }
    }
    
    // comb_ang_res[tiltx/y][angle,angres,angreserr][meas]
    
    for (int itilt=0; itilt<2; itilt++) {
        histname = TString("MEAS_tilt111_B");
        if (itilt==0) histname.ReplaceAll("111", "X");
        else histname.ReplaceAll("111", "Y");
        histname.ReplaceAll("MEAS", angres_histoname);
        
        gre = new TGraphErrors(i_point[itilt], comb_ang_res[itilt][0], comb_ang_res[itilt][1],angle_err,comb_ang_res[itilt][2]);
        gre->SetName(histname.Data());
        title= histname + TString(";Tilt angle (degrees);Angular Resolution (degrees)");
        gre->SetTitle(title.Data());
        gre->Sort();
        
        rootobjects.insert(pair<TString,TObject*>(histname,gre));
        
    }
    
    
    
}




void makeMultiGraps(){
    TGraphErrors *gre;
    TLatex *latex =new TLatex();
    latex->SetNDC();
    latex->SetTextFont(43);
    latex->SetTextColor(1);
    latex->SetTextSize(26);
    
    TMultiGraph *mgr_eff[2];
    mgr_eff[0]= new TMultiGraph("mgr_eff_tiltX","; Tilt Angle (degrees); Efficiency");
    mgr_eff[1]= new TMultiGraph("mgr_eff_tiltY","; Tilt Angle (degrees); Efficiency");
    //mgr_eff[0]= new TMultiGraph("mgr_eff_tiltX","Tilt around X; Tilt angle (degrees); Efficiency");
    //mgr_eff[1]= new TMultiGraph("mgr_eff_tiltY","Tilt around Y; Tilt angle (degrees); Efficiency");
    
    TMultiGraph *mgr_res[2];
    mgr_res[0]= new TMultiGraph("mgr_res_tiltX","; Tilt angle (degrees); Spatial Resolution (mm)");
    mgr_res[1]= new TMultiGraph("mgr_res_tiltY","; Tilt angle (degrees); Spatial Resolution (mm)");
    //mgr_res[0]= new TMultiGraph("mgr_res_tiltX","Tilt around X; Tilt angle (degrees); Spatial Resolution (mm)");
    //mgr_res[1]= new TMultiGraph("mgr_res_tiltY","Tilt around Y; Tilt angle (degrees); Spatial Resolution (mm)");
    
    TMultiGraph *mgr_angres= new TMultiGraph("mgr_angres","; Tilt Angle (degrees); Angular Resolution (degrees)");
    //TMultiGraph *mgr_angres= new TMultiGraph("mgr_angres","B Layer; Tilt angle (degrees); Angular Resolution (degrees)");
    TString histname;
    
    TLegend *leg_eff[2];
    TLegend *leg_res[2];
    TLegend *leg_angres;
    leg_angres = new TLegend(0.5, 0.9, 0.7, 0.75);
    leg_angres->SetFillColor(0);
    leg_angres->SetLineColor(0);
    leg_angres->SetTextSize(0.05);
    leg_angres->SetFillStyle(0);
    
    for (int i=0; i<2; i++) {
        leg_eff[i] = new TLegend(0.5, 0.9, 0.8, 0.7);
        leg_eff[i]->SetFillColor(0);
        leg_eff[i]->SetLineColor(0);
        leg_eff[i]->SetNColumns(2);
        leg_eff[i]->SetTextSize(0.05);
        leg_eff[i]->SetFillStyle(0);
        leg_res[i] = new TLegend(0.5, 0.9, 0.8, 0.7);
        leg_res[i]->SetFillColor(0);
        leg_res[i]->SetLineColor(0);
        leg_res[i]->SetNColumns(2);
        leg_res[i]->SetTextSize(0.05);
        leg_res[i]->SetFillStyle(0);
    }
    
    
    TString legname;
    int icolor[5]={0};
    for (map<TString,TObject*>::iterator it=rootobjects.begin(); it != rootobjects.end(); it++) {
        
        histname = it->first;
        legname = histname(histname.Last('_')+1, histname.Length()-histname.Last('_'));
        
        if(legname == TString("B-L2") ) legname = TString("B - L2");
        else if(legname == TString("B-L1") ) legname = TString("B - L1");
        else if(legname == TString("A-L2") ) legname = TString("A - L2");
        else if(legname == TString("A-L1") ) legname = TString("A - L1");
        
        if ( histname.Index(eff_histoname) != -1 && histname.Index("R") == -1){
            //&& histname.Index("R") == -1 && histname.Index("B-L2") == -1 && histname.Index("A-L1") == -1) { // pick eff plots
            gre = (TGraphErrors*) it->second;
            gre->SetMarkerStyle(20);
            gre->SetMarkerSize(1.6);
            gre->SetLineWidth(2);
            gre->SetFillStyle(0);
            
            if (histname.Index("tiltX") != -1) {
                gre->SetMarkerColor(color[icolor[0]]);
                gre->SetLineColor(color[icolor[0]]);
                icolor[0]++;
                mgr_eff[0]->Add(gre);
                leg_eff[0]->AddEntry(gre, legname.Data(), "ep");
                
                
            }
            else if (histname.Index("tiltY") != -1){
                gre->SetMarkerColor(color[icolor[1]]);
                gre->SetLineColor(color[icolor[1]]);
                icolor[1]++;
                mgr_eff[1]->Add(gre);
                leg_eff[1]->AddEntry(gre, legname.Data(), "ep");
                
            }
            
        }
        else if (histname.Index(res_histoname) != -1 && histname.Index(angres_histoname) == -1){ //&& histname.Index("R") == -1){ // res plots
            gre = (TGraphErrors*) it->second;
            gre->SetMarkerStyle(20);
            gre->SetLineWidth(2);
            gre->SetFillStyle(0);
            
            if (histname.Index("tiltX") != -1) {
                gre->SetMarkerColor(color[icolor[2]]);
                gre->SetLineColor(color[icolor[2]]);
                icolor[2]++;
                mgr_res[0]->Add(gre);
                leg_res[0]->AddEntry(gre, legname.Data(), "ep");
            }
            else if (histname.Index("tiltY") != -1){
                gre->SetMarkerColor(color[icolor[3]]);
                gre->SetLineColor(color[icolor[3]]);
                icolor[3]++;
                mgr_res[1]->Add(gre);
                leg_res[1]->AddEntry(gre, legname.Data(), "ep");
                
            }
        }
        else if (histname.Index(angres_histoname) != -1 && histname.Index("R") == -1){
            gre = (TGraphErrors*) it->second;
            gre->SetMarkerStyle(20);
            gre->SetMarkerSize(1.6);
            gre->SetLineWidth(2);
            gre->SetFillStyle(0);
            gre->SetMarkerColor(color[icolor[4]]);
            gre->SetLineColor(color[icolor[4]]);
            icolor[4]++;
            mgr_angres->Add(gre);
	    if(histname.Index("tiltX") != -1) legname = TString("Tilt around x-axis");
	    if(histname.Index("tiltY") != -1) legname = TString("Tilt around y-axis");
            leg_angres->AddEntry(gre, legname.Data(), "ep");
        }
        
    }
    TString pdfname, Cfilename;
    TCanvas *cc = new TCanvas("cc","",800,600);
    formatCanvas1D(cc);
    for (int i=0; i<2; i++) {
        mgr_eff[i]->SetMinimum(0.7);
        mgr_eff[i]->SetMaximum(1.25);
        
        histname = mgr_eff[i]->GetName();
        pdfname = TString("./results/") + histname +TString(".pdf");
        Cfilename = TString("./results/") + histname +TString(".C");
        mgr_eff[i]->Draw("AP");
        leg_eff[i]->Draw();
        formatMultiGraph(mgr_eff[i]);
        latex->DrawLatex(0.17,0.85,"DESY Testbeam");
        latex->DrawLatex(0.17,0.77,"4.4 GeV   2016");
        
        cc->SaveAs(pdfname.Data());
        cc->SaveAs(Cfilename.Data());
        rootobjects.insert(pair<TString,TObject*>(histname,mgr_eff[i]));
        
        mgr_res[i]->SetMinimum(0.);
        if (i==0) mgr_res[i]->SetMaximum(0.35);
        if (i==1) mgr_res[i]->SetMaximum(0.25);
        histname = mgr_res[i]->GetName();
        pdfname = TString("./results/") + histname +TString(".pdf");
        Cfilename = TString("./results/") + histname +TString(".C");
        mgr_res[i]->Draw("AP");
        leg_res[i]->Draw();
        mgr_res[i]->GetXaxis()->SetRangeUser(-1,11);
        formatMultiGraph(mgr_res[i]);
        latex->DrawLatex(0.17,0.85,"DESY Testbeam");
        latex->DrawLatex(0.17,0.77,"4.4 GeV   2016");
        
        cc->SaveAs(pdfname.Data());
        cc->SaveAs(Cfilename.Data());
        
        rootobjects.insert(pair<TString,TObject*>(histname,mgr_res[i]));
        
    }
    
    histname = mgr_angres->GetName();
    mgr_angres->SetMinimum(0);
    mgr_angres->SetMaximum(0.65);
    pdfname = TString("./results/") + histname +TString(".pdf");
    Cfilename = TString("./results/") + histname +TString(".C");
    mgr_angres->Draw("AP");
    leg_angres->Draw();
    mgr_angres->GetXaxis()->SetRangeUser(-1,9);
    formatMultiGraph(mgr_angres);
    
    leg_angres = new TLegend(0.5, 0.9, 0.7, 0.75);
    latex->DrawLatex(0.17,0.85,"DESY Testbeam");
    latex->DrawLatex(0.17,0.77,"4.4 GeV   2016");
    cc->SaveAs(pdfname.Data());
    cc->SaveAs(Cfilename.Data());
    
    
}


void createOutputFile(TFile *fout){
    
    vector<TString> created_folders;
    fout->cd();
    for (map<TString,TObject*>::iterator it=rootobjects.begin(); it != rootobjects.end(); it++) {
        fout->cd();
        
        Ssiz_t char_num = it->first.First("_");
        if (char_num != -1) {
            TString foldername(it->first(0,char_num));
            if (find(created_folders.begin(), created_folders.end(), foldername) == created_folders.end()){
                fout->mkdir(foldername.Data());
                created_folders.push_back(foldername);
            }
            fout->cd(foldername.Data());
        }
        else{
            fout->cd();
        }
        it->second->Write();
        
    }
    
}

void runHitMaker(){
    
    TString runthis;
    
    for (unsigned int irun=0; irun<runlist.size(); irun++) {
        runlist[irun].print();
        runthis = TString("./hitmaker ") + TString::Itoa(runlist[irun].DataRun,10) + TString(" ") + runlist[irun].config;
        cout<<runthis<<endl;
        system(runthis.Data());
    }
}
void runAlignment(){
    TString runthis;
    for (unsigned int irun=0; irun<runlist.size(); irun++) {
        runlist[irun].print();
        runthis = TString("./alignment ") + TString::Itoa(runlist[irun].DataRun,10);
        cout<<runthis<<endl;
        system(runthis.Data());
    }
    
}

void runAnalyze(){
    
    TString runthis;
    TString runnum_str;
    runthis = TString("rm ")+ output_efficiency_txtfile + TString(" ")+ output_resolution_txtfile;
    system(runthis.Data());
    
    std::ofstream out_eff;
    out_eff.open(output_efficiency_txtfile.Data(), std::ofstream::out | std::ofstream::app);
    /*
     std::map<unsigned short,TString> IDlayermap={
     {0,TString("B-L2")},{1,TString("B-L1")},
     {2,TString("A-L2")},{3,TString("A-L1")},
     {4,TString("R-L2")},{5,TString("R-L1")},
     };
     */
    out_eff<<"RunNum;";
    for (map<unsigned short,TString>::iterator it=IDlayermap.begin(); it != IDlayermap.end(); it++) {
        out_eff<<it->second<<";"<<it->second<<"_err;";
    }

    out_eff<<endl;
    out_eff.close();
    
    
    std::ofstream out_res;
    out_res.open(output_resolution_txtfile.Data(), std::ofstream::out | std::ofstream::app);
    out_res<<"RunNum;";
    for (map<unsigned short,TString>::iterator it=IDlayermap.begin(); it != IDlayermap.end(); it++) {
            out_res<<it->second<<";"<<it->second<<"_err;";
    }
    out_res<<"Angular;Angular_err;"<<endl;
    out_res.close();
    
    for (unsigned int irun=0; irun<runlist.size(); irun++) {
        runlist[irun].print();
        runnum_str = TString::Itoa(runlist[irun].DataRun,10);
        runthis = TString("./analyze ") + runnum_str + TString(" 134");
        //runthis = TString("./analyze ") + runnum_str + TString(" ")+ runnum_str;
        cout<<runthis<<endl;
        system(runthis.Data());
    }
    
    
}


void readMeasRes(){
    cout<< "Reading "<<output_resolution_txtfile<<endl;
    ifstream fileinfo_file(output_resolution_txtfile.Data());
    string line;
    getline(fileinfo_file,line);
    
    Meas* ameas;
    int runnum;
  //  for (map<unsigned short,TString>::iterator it=IDlayermap.begin(); it != IDlayermap.end(); it++) {

    while(getline(fileinfo_file,line)){
        vector<string> elems = splitstring(line,';');
        
        
        if(elems.size()!=15){
            cout<<"line in "<< output_resolution_txtfile<<" is not valid"<<endl;
            continue;
        }
        runnum = atoi(elems[0].c_str());
        if (runnum == 0) continue;
        
        ameas = dynamic_cast<Meas*> (measurements[runnum]);
        
        ameas->set_res_RL2(atof(elems[1].c_str()));
        ameas->set_res_RL2_err(atof(elems[2].c_str()));
        
        ameas->set_res_RL1(atof(elems[3].c_str()));
        ameas->set_res_RL1_err(atof(elems[4].c_str()));
        
        ameas->set_res_AL1(atof(elems[5].c_str()));
        ameas->set_res_AL1_err(atof(elems[6].c_str()));
        
        ameas->set_res_AL2(atof(elems[7].c_str()));
        ameas->set_res_AL2_err(atof(elems[8].c_str()));
        
        ameas->set_res_BL2(atof(elems[9].c_str()));
        ameas->set_res_BL2_err(atof(elems[10].c_str()));
        
        ameas->set_res_BL1(atof(elems[11].c_str()));
        ameas->set_res_BL1_err(atof(elems[12].c_str()));
        cout<<runnum<<" BL1 res err: "<<atof(elems[12].c_str())<<endl;
        
        ameas->set_ang_res(atof(elems[13].c_str()));
        ameas->set_ang_res_err(atof(elems[14].c_str()));
        
        measurements[runnum] = ameas;
        
    }
}


void readMeasEff(){
    cout<< "Reading "<<output_efficiency_txtfile<<endl;
    
    ifstream fileinfo_file(output_efficiency_txtfile.Data());
    string line;
    getline(fileinfo_file,line);
    
    Meas* ameas;
    int runnum;
    while(getline(fileinfo_file,line)){
        vector<string> elems = splitstring(line,';');
        
        if(elems.size()!=13){
            cout<<"line in "<< output_efficiency_txtfile<<" is not valid"<<endl;
            continue;
        }
        runnum = atoi(elems[0].c_str());
        if (runnum == 0) continue;

        ameas = dynamic_cast<Meas*> (measurements[runnum]);
        
        ameas->set_eff_RL2(atof(elems[1].c_str()));
        ameas->set_eff_RL2_err(atof(elems[2].c_str()));
        
        ameas->set_eff_RL1(atof(elems[3].c_str()));
        ameas->set_eff_RL1_err(atof(elems[4].c_str()));
        
        ameas->set_eff_AL1(atof(elems[5].c_str()));
        ameas->set_eff_AL1_err(atof(elems[6].c_str()));
        
        ameas->set_eff_AL2(atof(elems[7].c_str()));
        ameas->set_eff_AL2_err(atof(elems[8].c_str()));
        
        ameas->set_eff_BL2(atof(elems[9].c_str()));
        ameas->set_eff_BL2_err(atof(elems[10].c_str()));
        
        ameas->set_eff_BL1(atof(elems[11].c_str()));
        ameas->set_eff_BL1_err(atof(elems[12].c_str()));
        
        measurements[runnum] = ameas;
        
    }
    
}

void readRunList(TString filename){
    cout<< "Reading runlist "<<filename<<endl;
    
    vector<RunInfo> fileinfo;
    ifstream fileinfo_file(filename.Data());
    string firstline;
    getline(fileinfo_file,firstline);
    
    string s;
    while(!fileinfo_file.eof()){
        
        RunInfo info;
        
        getline(fileinfo_file,s,';');
        info.PedRun = atoi(s.c_str());
        
        getline(fileinfo_file,s,';');
        info.DataRun = atoi(s.c_str());
        
        getline(fileinfo_file,s,';');
        if (s == string("")) info.BeamEnergy = NOTSET;
        else info.BeamEnergy = atof(s.c_str());
        
        getline(fileinfo_file,s,';');
        if (s == string("")) info.Distance1 = 0;
        else info.Distance1 = atof(s.c_str());
        
        getline(fileinfo_file,s,';');
        if (s == string("")) info.Angle1 = 0;
        else info.Angle1 = atof(s.c_str());
        
        getline(fileinfo_file,s,';');
        if (s == string("")) info.Distance2 = 0;
        else info.Distance2 = atof(s.c_str());
        
        getline(fileinfo_file,s,';');
        if (s == string("")) info.Angle2 = 0;
        else info.Angle2 = atof(s.c_str());
        
        getline(fileinfo_file,info.config,';');
        getline(fileinfo_file,info.map,';');
        getline(fileinfo_file,info.Comment,';');
        
        if (info.PedRun != 0 && info.DataRun !=0 && info.BeamEnergy !=0) {
            runlist.push_back(info);
            
            Meas* ameas = new Meas();
            
            ameas->set_RunNum(info.DataRun);
            ameas->set_Angle1(info.Angle1);
            ameas->set_Angle2(info.Angle2);
            ameas->set_BeamEnergy(info.BeamEnergy);
            
            
            measurements.insert(make_pair(info.DataRun, ameas));
            //info.print();
        }
    }
    return;
}
