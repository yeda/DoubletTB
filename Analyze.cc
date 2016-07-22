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

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TObject.h"
#include "TROOT.h"
#include "TMath.h"

#include "Settings.h"
#include "Analyze.h"

TString input_tree_name = TString("hitsTree");

TString hitpos_histname = TString("hitposition_");
TString hitmap_histname = TString("hitmap_");
TString hitamp_histname = TString("hitamplitude_");
TString corr_histname = TString("corr_");
TString spatialRes_histname = TString("spatialRes_");
TString angularRes_histname = TString("angularRes_");

map<TString,TObject*> rootobjects;

unsigned int event_num;
vector<unsigned short> *layerID;
vector<unsigned short> *cluster_size;
vector<float> *hit_amplitude;
vector<float> *hit_position;

double alignmentpar[NLayer];

// for efficieancy calculation
unsigned int layercount[NLayer];
unsigned int expectedcount[NLayer];

int main(int argc, char *argv[]){

    // get tree from input file
    TString runnum = TString(argv[1]);
    TString fin_name = outputPath + TString("hitmaker_") + runnum + TString(".root");
    TFile * fin = TFile::Open(fin_name.Data());
    TTree* input_tree = (TTree*)fin->Get(input_tree_name.Data());
    setBranches(input_tree);

    TString fout_name = outputPath + TString("analyze_") + runnum + TString(".root");
    TFile *fout = new TFile(fout_name.Data(),"RECREATE");
    
    // get alignment parameters
    TString alignmentrun = TString(argv[2]);
    TString alignmenttxt = outputPath + TString("alignment_") + alignmentrun + TString(".txt");
    readAlignmentParameters(alignmenttxt, alignmentpar);

    
    // create histograms
    createHistos();

    // run through events
    Long64_t nentries = input_tree->GetEntries();

    for (Long64_t ientry=0; ientry<nentries;ientry++) {
        
        input_tree->GetEntry(ientry);
        
        calculateEfficiency();
        if(isGoodEvent()){

            fillHistos();
            
        }
    }
 
    printResolution(runnum);
    printEfficiency(runnum);
    createOutputFile(fout);
    
    return 0;
}


void findRefLayers(int layerid, int* arr){
    TString ilayerName = IDlayermap[layerid];
    
    if (ilayerName.First("X") != -1) {
        int reflayercount=0;
        for (unsigned int i=0; i<3; i++) {
            if (xlayers[i]!=layerid) {
                arr[reflayercount] = xlayers[i];
                reflayercount++;
            }
        }
    }
    else if (ilayerName.First("Y") != -1) {
        int reflayercount=0;
        for (unsigned int i=0; i<3; i++) {
            if (ylayers[i]!=layerid) {
                arr[reflayercount] = ylayers[i];
                reflayercount++;
            }
        }
    }
    else {
        cout<<"Reference layers of layer "<<ilayerName<<" cannot be found!"<<endl;
    }
    
    
}


void calculateEfficiency(){
    int reflayers[2];
    vector<unsigned short> reflayershit[2];
    vector<unsigned short> layerhit;
    double  x1, x2, y1, y2, x, y;
    
    for (unsigned short ilayer=0; ilayer<NLayer; ilayer++) {
        findRefLayers(ilayer,reflayers);
        //cout<<"Layer: "<< ilayer<<" Ref layers "<< reflayers[0]<< " "<<reflayers[1]<<endl;
        
        // count number of hits in each referance layer
        reflayershit[0].clear(); reflayershit[1].clear();
        layerhit.clear();
        for (unsigned short ihit=0; ihit<(layerID->size()); ihit++) {
            if (layerID->at(ihit) == ilayer) {
                layerhit.push_back(ihit);
            }
            else if (layerID->at(ihit) == reflayers[0]) {
                reflayershit[0].push_back(ihit);
            }
            else if (layerID->at(ihit) == reflayers[1]) {
                reflayershit[1].push_back(ihit);
            }
        } // end of loop over hits in each referance layer
        
        
        if (reflayershit[0].size() == 1 && reflayershit[1].size() == 1) {
            // calculate expected hit position on the layer
            
            // line eq -> ax+by+c=0
            // (y1-y2) * x + (x2-x1) * y + (x1y2-x2y1) = 0
            unsigned int ihit1, ihit2;
            ihit1 = reflayershit[0].at(0);
            ihit2 = reflayershit[1].at(0);

            x1 = hit_position->at(ihit1) - alignmentpar[layerID->at(ihit1)];
            x2 = hit_position->at(ihit2) - alignmentpar[layerID->at(ihit2)];
            y1 = layerZposition[layerID->at(ihit1)];
            y2 = layerZposition[layerID->at(ihit2)];

            y = layerZposition[ilayer];
            x = getExpectedHit(x1,y1,x2,y2,y);
        
            unsigned int n_pairs = 0;
            // check if there is a hit there
            
            for (unsigned int jj=0; jj<layerhit.size(); jj++) {
                double dpos = x - (hit_position->at(layerhit[jj])- alignmentpar[layerID->at(layerhit[jj])]);
                if ( fabs(dpos) < layerResolution[ilayer]){
                    n_pairs++;
                }
            }
            
            expectedcount[ilayer]++;
            if (n_pairs>0) layercount[ilayer]++;
        }
        
        
    }
}

double getExpectedHit(double x1, double y1, double x2, double y2, double at_y){
    double eq_a, eq_b, eq_c, x;
    // calculate expected hit position on the layer
    
    // line eq -> ax+by+c=0
    // (y1-y2) * x + (x2-x1) * y + (x1y2-x2y1) = 0
    
    eq_a = y1-y2;
    eq_b = x2-x1;
    eq_c = x1*y2 - x2*y1;
    
    x = -1*(eq_c + at_y * eq_b) / eq_a;
    
    return x;
}

void fillHistos(){
    
    TString histname,layername;
    TH1D *h1D;
    TH2D* h2D;

    double newpos;
    // [pointnumber][x,y,z]
    double p[3][3];
    
    for (unsigned int ihit=0; ihit<layerID->size(); ihit++) {
        
        layername = IDlayermap[layerID->at(ihit)];
        
        histname = hitpos_histname + layername;
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        newpos = hit_position->at(ihit) - alignmentpar[layerID->at(ihit)];
        h1D->Fill(newpos);

        histname = hitamp_histname + layername;
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(hit_amplitude->at(ihit));

        for (unsigned int i=0; i<3; i++) {
            if (layerID->at(ihit) == xlayers[i]) {
                p[i][0] = newpos;
                p[i][2] = layerZposition[layerID->at(ihit)];
            }
            else if (layerID->at(ihit) == ylayers[i])
                p[i][1] = newpos;
        }

    }
    
    for (unsigned int i=1; i<3; i++) {
        // correlation
        histname = corr_histname + IDlayermap[xlayers[i]];
        h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
        h2D->Fill(p[0][0], p[i][0]);

        histname = corr_histname + IDlayermap[ylayers[i]];
        h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
        h2D->Fill(p[0][1], p[i][1]);

        // spatial resolution
        histname = spatialRes_histname + IDlayermap[xlayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(p[0][0] - p[i][0]);

        histname = spatialRes_histname + IDlayermap[ylayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(p[0][1] - p[i][1]);

    }
    
    for (unsigned int i=0; i<3; i++) {
        layername = IDlayermap[xlayers[i]];
        TString tempname = layername(0, layername.Length() -1);
        histname = hitmap_histname + tempname;
        h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
        h2D->Fill(p[i][0], p[i][1]);
    }

    // Angular resolution
    
    // match pointnumber with layers
    // it is ordered as xlayers[i]
    //      so p[0] is ref, p[1] is Up, p[2] is Down
    double p_exp[3];
    p_exp[2] = p[2][2];
    p_exp[0] = getExpectedHit(p[0][0],p[0][2],p[1][0],p[1][2],p_exp[2]);
    p_exp[1] = getExpectedHit(p[0][1],p[0][2],p[1][1],p[1][2],p_exp[2]);
    
    
    histname = angularRes_histname + TString("Down");
    h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
    // angle between down_meas - up_meas - down_exp
    h1D->Fill(getAngleABC(p[2],p[1],p_exp));

    
}



double getAngleABC(double A[] , double B[], double C[]){
    // measuring angle on side B
    
    //In pseudo-code, the vector BA (call it v1) is:
    double v1[3]; // it is actually a vector
    v1[0] = A[0] - B[0];
    v1[1] = A[1] - B[1];
    v1[2] = A[2] - B[2];
    //Similarly the vector BC (call it v2) is:
    double v2[3];
    v2[0] = C[0] - B[0];
    v2[1] = C[1] - B[1];
    v2[2] = C[2] - B[2];
    //The dot product of v1 and v2 is a function of the cosine of the angle between them (it's scaled by the product of their magnitudes). So first normalize v1 and v2:
    
    double v1mag = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
    double v1norm[3];
    v1norm[0] = v1[0] / v1mag;
    v1norm[1] = v1[1] / v1mag;
    v1norm[2] = v1[2] / v1mag;
    
    double v2mag = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
    double v2norm[3];
    v2norm[0] = v2[0] / v2mag;
    v2norm[1] = v2[1] / v2mag;
    v2norm[2] = v2[2] / v2mag;
    
    //Then calculate the dot product:
    
    double res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2];
    //And finally, recover the angle:
    
    double angle = acos(res);
    // convert to degrees
    
    return (angle*180.0)/TMath::Pi();
    
}

void printEfficiency(TString runnum){
    std::ofstream outfile;
    
    outfile.open(output_efficiency_txtfile.Data(), std::ofstream::out | std::ofstream::app);
    
    outfile << runnum << ";";
    for (unsigned short i_layer=0; i_layer<NLayer; i_layer++){
        double eff = double (layercount[i_layer])/double (expectedcount[i_layer]);
        cout<<"Layer: "<<IDlayermap[i_layer]<<" efficiency: "<<eff<<"       found hits: "<<layercount[i_layer]<<" expected hits: "<<expectedcount[i_layer]<<endl;
        outfile<< eff << ";";
    }
    outfile<<endl;
    outfile.close();
}




void createOutputFile(TFile *fout){
    fout->cd();
    
    for (map<TString,TObject*>::iterator it=rootobjects.begin(); it != rootobjects.end(); it++) {
        it->second->Write();
    }
    
}

void readAlignmentParameters(TString alignmenttxt, double alignmentpar[]){
    ifstream alignmentFile (alignmenttxt.Data());
    if( ! alignmentFile.is_open() ){
        cout<< "Unable to open file: " << alignmenttxt << endl;
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

void printResolution(TString runnum){
    TString layername, histname, fitname;
    TF1* fitfunc;
    TH1D* h1D;
    double max_inxaxis;

    std::ofstream outfile;
    
    outfile.open(output_resolution_txtfile.Data(), std::ofstream::out | std::ofstream::app);
    outfile << runnum << ";";

    for (unsigned int i=1; i<3; i++) {
        // point resolution
        for (unsigned int j=0; j<2; j++) {
            if (j==0) layername =IDlayermap[xlayers[i]];
            if (j==1) layername =IDlayermap[ylayers[i]];
            
            histname =  spatialRes_histname + layername;
            h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
            max_inxaxis =  h1D->GetXaxis()->GetBinCenter(h1D->GetMaximumBin());
            fitname = TString("fit_") + histname;
            fitfunc = new TF1(fitname.Data(),"gaus", max_inxaxis-0.5,max_inxaxis+0.5);
            h1D->Fit(fitfunc,"QR");
            
            double spat_res =fitfunc->GetParameter(2) / sqrt(2);
            double spat_res_err =fitfunc->GetParError(2) / sqrt(2);
            cout<< "spatial resolution of "<< layername<< " is "<< spat_res <<" +/- "<<spat_res_err<< " (mm)"<< endl;
            outfile<<spat_res<<";"<<spat_res_err<<";";
        }
    }
    
    histname = angularRes_histname + TString("Down");
    h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
    max_inxaxis =  h1D->GetXaxis()->GetBinCenter(h1D->GetMaximumBin());
    fitname = TString("fit_") + histname;
    fitfunc = new TF1(fitname.Data(),"gaus", max_inxaxis-0.2,max_inxaxis+0.2);
    h1D->Fit(fitfunc,"QR");
    
    double ang_res = fitfunc->GetParameter(2) / sqrt(2);
    double ang_res_err = fitfunc->GetParError(2) / sqrt(2);
    cout<< "angular resolution of Down is "<< ang_res <<" +/- "<< ang_res_err<<" degrees"<<endl;
    outfile<<ang_res<<";"<<ang_res_err<<";"<<endl;
    outfile.close();

}

void createHistos(){
    TString histname,title;
    TH1D* h1D;
    TH2D* h2D;
    
    for (map<unsigned short,TString>::iterator it=IDlayermap.begin(); it != IDlayermap.end(); it++) {
        
        // hit position
        histname = hitpos_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Hit position of all hits in ") + it->second + TString(";Hit Position (mm);Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 1000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        // hit amplitude
        histname = hitamp_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Hit amplitude of all clusters in ") + it->second + TString(";Hit Amplitude;Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 1000,0, 10000);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
    }
    
    for (unsigned int i=1; i<3; i++) {
        // point resolution
        histname = spatialRes_histname + IDlayermap[xlayers[i]];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Spatial Resolution;")+IDlayermap[xlayers[0]]+TString(" - ")+IDlayermap[xlayers[i]]+TString(" (mm);Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 10000,-100, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        histname = spatialRes_histname + IDlayermap[ylayers[i]];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Spatial Resolution;")+IDlayermap[ylayers[0]]+TString(" - ")+IDlayermap[ylayers[i]]+TString(" (mm);Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 10000, -100, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        
        
        // correlation
        histname = corr_histname + IDlayermap[xlayers[i]];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Correlation;") + IDlayermap[xlayers[0]] + TString(";") + IDlayermap[xlayers[i]];
            h2D = new TH2D(histname.Data(), title.Data(), 1000,0, 100, 1000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h2D));
        }
        histname = corr_histname + IDlayermap[ylayers[i]];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Correlation;")+IDlayermap[ylayers[0]]+TString(";") + IDlayermap[ylayers[i]];
            h2D = new TH2D(histname.Data(), title.Data(), 1000,0, 100, 1000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h2D));
        }
        
    }
    // angular resolution
    histname = angularRes_histname + TString("Down");
    if (rootobjects.find(histname) == rootobjects.end()) {
        title = TString("Angular Resolution (Down); Angle difference (degrees);Number of entries");
        h1D = new TH1D(histname.Data(), title.Data(), 1000,-10, 10);
        rootobjects.insert(pair<TString,TObject*>(histname,h1D));
    }
    
    TString s[3] = {"Down","Up","Ref"};
    for (unsigned int i=0; i<3; i++) {
        // point resolution
        histname = hitmap_histname + s[i];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("HitMap ")+ s[i]+TString(";X (mm);Y (mm)");
            h2D = new TH2D(histname.Data(), title.Data(), 1000,0, 100,1000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h2D));
        }
    }
    
    
    /*
     histname = angularRes_histname + TString("Up");
     if (rootobjects.find(histname) == rootobjects.end()) {
     title = TString("Angular Resolution (Up); Angle difference (degrees);Number of entries");
     h1D = new TH1D(histname.Data(), title.Data(), 1000,-10, 10);
     rootobjects.insert(pair<TString,TObject*>(histname,h1D));
     }
     */
    
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


void setBranches(TTree* input_tree){
    
    input_tree->SetBranchAddress("event_num",&event_num);
    input_tree->SetBranchAddress("layerID",&layerID);
    input_tree->SetBranchAddress("cluster_size",&cluster_size);
    input_tree->SetBranchAddress("hit_amplitude",&hit_amplitude);
    input_tree->SetBranchAddress("hit_position",&hit_position);
}




