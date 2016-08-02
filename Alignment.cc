/*
 * Created by Eda Yildirim
 *  (2016 JGU)
 *
 *  email:eda.yildirim@cern.ch
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <stdlib.h>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TObject.h"
#include "TROOT.h"

#include "Settings.h"
#include "Alignment.h"


TString input_tree_name = TString("hitsTree");

TString hitpos_histname = TString("hitposition_");
TString hitamp_histname = TString("hitamplitude_");
TString dx_histname = TString("dx_");
TString dy_histname = TString("dy_");

TString aligned_hitpos_histname = TString("aligned_hitposition_");
TString aligned_hitamp_histname = TString("aligned_hitamplitude_");
TString aligned_corr_histname = TString("aligned_corr_");
TString aligned_dx_histname = TString("aligned_dx_");
TString aligned_dy_histname = TString("aligned_dy_");

TString output_alignmenttxt = TString("alignment_");
TString output_alignmentroot = TString("alignment_");

map<TString,TObject*> rootobjects;

unsigned int event_num;
vector<unsigned short> *layerID;
vector<unsigned short> *cluster_size;
vector<float> *hit_amplitude;
vector<float> *hit_position;

int main(int argc, char *argv[]){
    
    TString runnum = TString(argv[1]);
    
    TString fin_name = outputPath + TString("hitmaker_") + runnum + TString(".root");
    
    TFile * fin = TFile::Open(fin_name.Data());
    TTree* input_tree = (TTree*)fin->Get(input_tree_name.Data());
    setBranches(input_tree);
    
    output_alignmenttxt = outputPath + output_alignmenttxt + runnum + TString(".txt");
    
    output_alignmentroot = outputPath + output_alignmentroot + runnum + TString(".root");
    TFile *fout = new TFile(output_alignmentroot.Data(),"RECREATE");
    
    
    // create histograms
    createHistos();
    
    
    Long64_t nentries = input_tree->GetEntries();
    
    for (Long64_t ientry=0; ientry<nentries;ientry++) {
        
        input_tree->GetEntry(ientry);
        
        if(isGoodEvent()){
            fillHistos();
        }
    }
    
    
    alignLayers();
    
    if (PlotHistosAfterAlignment) {
        
        double alignmentpar[NLayer];
        readAlignmentParameters(alignmentpar);
        
        for (Long64_t ientry=0; ientry<nentries;ientry++) {
            
            input_tree->GetEntry(ientry);
            
            if(isGoodEvent()){
                fillHistosAfterAlignment(alignmentpar);
            }
        }
    }
    createOutputFile(fout);
    
    return 0;
}
void createOutputFile(TFile *fout){
    fout->cd();
    
    for (map<TString,TObject*>::iterator it=rootobjects.begin(); it != rootobjects.end(); it++) {
        it->second->Write();
    }
    
}


void readAlignmentParameters(double alignmentpar[]){
    ifstream alignmentFile (output_alignmenttxt.Data());
    if( ! alignmentFile.is_open() ){
        cout<< "Unable to open file: " << output_alignmenttxt << endl;
        return;
    }
    
    string alignmentline, a_layername;
    int a_layerID;
    double a_alignmentpar;
    getline(alignmentFile,alignmentline);
    while ( getline(alignmentFile,alignmentline) ){
        
        size_t charpos, prev_charpos;
        // layer ID
        charpos = alignmentline.find(";");
        a_layerID = stoi( alignmentline.substr(0,charpos) );
        // layer name
        prev_charpos = charpos;
        charpos = alignmentline.find(";",prev_charpos+1);
        a_layername = alignmentline.substr(prev_charpos+1,charpos-prev_charpos-1);
        // alignment parameter
        prev_charpos = charpos;
        charpos = alignmentline.find(";",prev_charpos+1);
        
        a_alignmentpar = stod( alignmentline.substr(prev_charpos+1,charpos-prev_charpos-1));
        alignmentpar[a_layerID] = a_alignmentpar;
    }
    
}

void alignLayers(){
    cout<< "Alignment Parameters"<<endl;
    
    TString layername, histname, fitname;
    TF1* fitfunc;
    TH1D* h1D;
    
    double mean_pos[2][3];
    for (unsigned int i=0; i<3; i++) {
        if (i==0) {
            mean_pos[0][i]=0;
            mean_pos[1][i]=0;
        }
        else{
            layername = IDlayermap[xlayers[i]];
            histname =  dx_histname + layername;
            h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
            fitname = TString("fit_") + histname;
            fitfunc = new TF1(fitname.Data(),"gaus");
            h1D->Fit(fitfunc,"Q");
            mean_pos[0][i]=fitfunc->GetParameter(1);
            
            layername = IDlayermap[ylayers[i]];
            histname =  dy_histname + layername;
            h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
            fitname = TString("fit_") + histname;
            fitfunc = new TF1(fitname.Data(),"gaus");
            h1D->Fit(fitfunc,"Q");
            mean_pos[1][i]=fitfunc->GetParameter(1);
            
        }
    }
    
    // Write alignment parameters
    ofstream outtxt (output_alignmenttxt.Data());
    if( ! outtxt.is_open() ){
        cout<< "Unable to open file: " << output_alignmenttxt << endl;
        return;
    }
    cout<<"LayerID; LayerName; AlignmentShift;"<<endl;
    outtxt<<"LayerID;LayerName;AlignmentShift;"<<endl;
    
    double alignmentshift=0;
    for (unsigned int i=0; i<3; i++) {
        // xlayers
        alignmentshift = mean_pos[0][i];
        cout<<xlayers[i]<<"; "<<IDlayermap[xlayers[i]]<<"; "<< alignmentshift<<";"<<endl;
        outtxt<<xlayers[i]<<";"<<IDlayermap[xlayers[i]]<<";"<< alignmentshift<<";"<<endl;
        
        // ylayers
        alignmentshift = mean_pos[1][i];
        cout<<ylayers[i]<<"; "<<IDlayermap[ylayers[i]]<<"; "<< alignmentshift<<";"<<endl;
        outtxt<<ylayers[i]<<";"<<IDlayermap[ylayers[i]]<<";"<< alignmentshift<<";"<<endl;
    }
    
    outtxt.close();
}

void fillHistosAfterAlignment(double alignmentpar[]){
    TString layername, histname;
    TH1D* h1D;
    TH2D* h2D;
    
    double newpos;
    
    double px[3];
    double py[3];
    
    for (unsigned int ihit=0; ihit<layerID->size(); ihit++) {
        
        layername = IDlayermap[layerID->at(ihit)];
        
        histname = aligned_hitpos_histname + layername;
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        newpos = hit_position->at(ihit) - alignmentpar[layerID->at(ihit)];
        h1D->Fill(newpos);
        
        for (unsigned int i=0; i<3; i++) {
            if (layerID->at(ihit) == xlayers[i]) {
                px[i] = newpos;
            }
            else if (layerID->at(ihit) == ylayers[i])
                py[i] = newpos;
        }
        
    }
    
    for (unsigned int i=1; i<3; i++) {
        // correlation
        
        histname = aligned_corr_histname + IDlayermap[xlayers[i]];
        h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
        h2D->Fill(px[0], px[i]);
        
        histname = aligned_corr_histname + IDlayermap[ylayers[i]];
        h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
        h2D->Fill(py[0], py[i]);
        
        histname = aligned_dx_histname + IDlayermap[xlayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(px[i]-px[0]);
        
        histname = aligned_dy_histname + IDlayermap[ylayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(py[i]-py[0]);
        
    }
    
}


void fillHistos(){
    TString layername, histname;
    TH1D* h1D;
    TH2D* h2D;
    
    double newpos;
    
    double px[3];
    double py[3];
    
    for (unsigned int ihit=0; ihit<layerID->size(); ihit++) {
        
        layername = IDlayermap[layerID->at(ihit)];
        
        histname = hitpos_histname + layername;
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(hit_position->at(ihit));
        
        histname = hitamp_histname + layername;
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(hit_amplitude->at(ihit));
        
        for (unsigned int i=0; i<3; i++) {
            if (layerID->at(ihit) == xlayers[i]) {
                px[i] = hit_position->at(ihit);
            }
            else if (layerID->at(ihit) == ylayers[i])
                py[i] = hit_position->at(ihit);
        }
        
    }
    
    for (unsigned int i=1; i<3; i++) {
        // correlation
        
        histname = dx_histname + IDlayermap[xlayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(px[i]-px[0]);
        
        histname = dy_histname + IDlayermap[ylayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(py[i]-py[0]);
        
    }
    
    
    
    
}


bool isGoodEvent(){
    
    short nhits[NLayer]={0};
    bool onlyonehit = true;
    
    // also check if every layer has only 1 cluster
    for (unsigned int i=0; i<layerID->size(); i++) {
        nhits[layerID->at(i)]++;
    }
    // well now all values should be 1
    for (unsigned int i=0; i<NLayer; i++) {
        if (nhits[i]!=1) onlyonehit = false;
    }
    
    return onlyonehit;
}


void createHistos(){
    TString histname, title;
    TH1D* h;
    TH1D* h1D;
    TH2D* h2D;
    
    for (map<unsigned short,TString>::iterator it=IDlayermap.begin(); it != IDlayermap.end(); it++) {
        
        
        histname = hitpos_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Hit position of all hits in ") + it->second + TString(";Hit Position (mm);Number of entries");
            h = new TH1D(histname.Data(), title.Data(), 1000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h));
        }
        histname = hitamp_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Hit amplitude of all clusters in ") + it->second + TString(";Hit Amplitude;Number of entries");
            h = new TH1D(histname.Data(), title.Data(), 1000,0, 10000);
            rootobjects.insert(pair<TString,TObject*>(histname,h));
        }
        
        if (PlotHistosAfterAlignment) {
            histname = aligned_hitpos_histname + it->second;
            if (rootobjects.find(histname) == rootobjects.end()) {
                title = TString("Hit position (after alignment) of all hits in ")+ it->second + TString(";Hit Position (mm);Number of entries");
                h = new TH1D(histname.Data(), title.Data(), 1000,0, 100);
                rootobjects.insert(pair<TString,TObject*>(histname,h));
            }
            
        }
        
    }
    
    for (unsigned int i=1; i<3; i++) {
        // correlation
        histname = dx_histname + IDlayermap[xlayers[i]];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("dx;") + IDlayermap[xlayers[i]] +TString("-")+ IDlayermap[xlayers[0]] + TString("; Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 1000,-50, 50);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        histname = dy_histname + IDlayermap[ylayers[i]];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("dy;") + IDlayermap[ylayers[i]] +TString("-")+ IDlayermap[ylayers[0]] + TString("; Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 1000,-50, 50);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        if (PlotHistosAfterAlignment) {
            // correlation
            histname = aligned_corr_histname + IDlayermap[xlayers[i]];
            if (rootobjects.find(histname) == rootobjects.end()) {
                title = TString("Aligned Correlation;") + IDlayermap[xlayers[0]] + TString(";") + IDlayermap[xlayers[i]];
                h2D = new TH2D(histname.Data(), title.Data(), 1000,0, 100, 1000,0, 100);
                rootobjects.insert(pair<TString,TObject*>(histname,h2D));
            }
            histname = aligned_corr_histname + IDlayermap[ylayers[i]];
            if (rootobjects.find(histname) == rootobjects.end()) {
                title = TString("Aligned Correlation;")+IDlayermap[ylayers[0]]+TString(";") + IDlayermap[ylayers[i]];
                h2D = new TH2D(histname.Data(), title.Data(), 1000,0, 100, 1000,0, 100);
                rootobjects.insert(pair<TString,TObject*>(histname,h2D));
            }
            histname = aligned_dx_histname + IDlayermap[xlayers[i]];
            if (rootobjects.find(histname) == rootobjects.end()) {
                title = TString("dx;") + IDlayermap[xlayers[i]] +TString("-")+ IDlayermap[xlayers[0]] + TString("; Number of entries");
                h1D = new TH1D(histname.Data(), title.Data(), 100,-20, 20);
                rootobjects.insert(pair<TString,TObject*>(histname,h1D));
            }
            histname = aligned_dy_histname + IDlayermap[ylayers[i]];
            if (rootobjects.find(histname) == rootobjects.end()) {
                title = TString("dy;") + IDlayermap[ylayers[i]] +TString("-")+ IDlayermap[ylayers[0]] + TString("; Number of entries");
                h1D = new TH1D(histname.Data(), title.Data(), 100,-20, 20);
                rootobjects.insert(pair<TString,TObject*>(histname,h1D));
            }
            
        }
        
    }
    
    
    
}

void setBranches(TTree* input_tree){
    
    input_tree->SetBranchAddress("event_num",&event_num);
    input_tree->SetBranchAddress("layerID",&layerID);
    input_tree->SetBranchAddress("cluster_size",&cluster_size);
    input_tree->SetBranchAddress("hit_amplitude",&hit_amplitude);
    input_tree->SetBranchAddress("hit_position",&hit_position);
}




