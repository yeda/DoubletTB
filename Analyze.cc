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
TString exphitmap_histname = TString("exphitmap_");
TString exphitpos_histname = TString("exphitpos_");
TString hitamp_histname = TString("hitamplitude_");
TString hitamp2D_histname = TString("hitamplitude2D_");
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
unsigned int layercount[NLayer]={0};
unsigned int expectedcount[NLayer]={0};

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
        fillHitMap2D();
        if(isGoodEvent()){
            fillHistos();
        }
    }
    //
    
    // modified res plots, mean shifted to zero
    for (Long64_t ientry=0; ientry<nentries;ientry++) {
        input_tree->GetEntry(ientry);
        if(isGoodEvent()){
            fillModifiedHistos();
        }
    }
    
    printResolution(runnum);
    printEfficiency(runnum);
    
    createOutputFile(fout);
    
    return 0;
}

void findRefLayers(int layerid, int* arr){
    
    // find if Layer is in x axis
    bool is_x_axis=false;
    for (unsigned int i=0; i<3; i++) {
        if (xlayers[i]==layerid) is_x_axis = true;
    }
    
    int reflayercount=0;
    
    if (is_x_axis) {
        for (unsigned int i=0; i<3; i++) {
            if (xlayers[i]!=layerid) {
                arr[reflayercount] = xlayers[i];
                reflayercount++;
            }
        }
    }
    else {
        for (unsigned int i=0; i<3; i++) {
            if (ylayers[i]!=layerid) {
                arr[reflayercount] = ylayers[i];
                reflayercount++;
            }
        }
    }
    
}



void calculateEfficiency(){
    int reflayersX[2];
    int reflayersY[2];
    
    map< int, vector<float> > hit;
    double  x1, x2,x3, y1, y2,y3, z1,z2,z3, x, y;
    TH1D *h1D;
    TH2D *h2D;
    TString histname;
    
    double newpos;
    for (int i=0; i<NLayer; i++) {
        vector<float> bla;
        bla.clear();
        hit[i] = bla;
    }
    
    for (unsigned int ihit=0; ihit<layerID->size(); ihit++) {
        newpos = hit_position->at(ihit) - alignmentpar[layerID->at(ihit)];
        hit[layerID->at(ihit)].push_back(newpos);
    }
    
    
    for (int i_layer=0; i_layer<3; i_layer++) {
        //current xlayer = xlayers[i_layer]
        //current ylayer = ylayers[i_layer]
        findRefLayers(xlayers[i_layer],reflayersX);
        findRefLayers(ylayers[i_layer],reflayersY);
        
        bool reflayershave1hit=true;
        for (int i_ref=0; i_ref<2; i_ref++) {
            if ( hit[reflayersX[i_ref]].size() !=1 || hit[reflayersY[i_ref]].size() !=1 ) reflayershave1hit = false;
        }
        
        if (reflayershave1hit) {
            expectedcount[ xlayers[i_layer] ]++;
            expectedcount[ ylayers[i_layer] ]++;
            
            // calculate expected hit position on the layer
            
            x1 = hit[ reflayersX[0] ][0];
            y1 = hit[ reflayersY[0] ][0];
            z1 = layerZposition[ reflayersX[0] ];
            
            x2 = hit[ reflayersX[1] ][0];
            y2 = hit[ reflayersY[1] ][0];
            z2 = layerZposition[ reflayersX[1] ];
            
            vector<float> xhits = hit[ xlayers[i_layer] ];
            vector<float> yhits = hit[ ylayers[i_layer] ];
            z3 = zpos[i_layer];
            unsigned int n_pairs_x = 0;
            unsigned int n_pairs_y = 0;
            
            // X layer
            for (int i_hit_x=0; i_hit_x<xhits.size(); i_hit_x++) {
                x3 = xhits[i_hit_x];
                x = getExpectedHit(x1,z1,x2,z2,z3);
                if ( fabs(x-x3)<layerResolution[ xlayers[i_layer] ] ) {
                    n_pairs_x++;
                }
            }
            if (n_pairs_x>0) layercount[ xlayers[i_layer] ]++;
            else{
                x = getExpectedHit(x1,z1,x2,z2,z3);
                y = getExpectedHit(y1,z1,y2,z2,z3);
                histname = exphitmap_histname + IDlayermap[xlayers[i_layer]];
                h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
                h2D->Fill(x,y);
                
            }
            
            // Y layer
            for (int i_hit_y=0; i_hit_y<yhits.size(); i_hit_y++) {
                y3 = yhits[i_hit_y];
                y = getExpectedHit(y1,z1,y2,z2,z3);
                if ( fabs(y-y3)<layerResolution[ ylayers[i_layer] ] ) {
                    n_pairs_y++;
                }
            }
            if (n_pairs_y>0) layercount[ ylayers[i_layer] ]++;
            else{
                x = getExpectedHit(x1,z1,x2,z2,z3);
                y = getExpectedHit(y1,z1,y2,z2,z3);
                histname = exphitmap_histname + IDlayermap[ylayers[i_layer]];
                h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
                h2D->Fill(x,y);
                
            }
            
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

void fillHitMap2D(){
    TString histname,layername,tempname;
    TH2D* h2D;
    float x,y;
    vector<float> vec_x;
    vector<float> vec_y;
    
    map< int, vector<float> > hit = getHits("hitamp");
    
    for (int i_layer=0; i_layer<3; i_layer++) {
        
        vec_x = hit[xlayers[i_layer]];
        vec_y = hit[ylayers[i_layer]];
        
        
        
        layername = IDlayermap[xlayers[i_layer]];
        tempname = layername(0, 1);
        histname = hitamp2D_histname + tempname;
        h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
        
        if (vec_x.size() != 0 && vec_y.size() != 0) {
            for (int i_x=0; i_x<vec_x.size(); i_x++) {
                x = vec_x[i_x];
                for (int i_y=0; i_y<vec_y.size(); i_y++) {
                    y = vec_y[i_y];
                    h2D->Fill(x,y);
                }
            }
        }
        else if(vec_x.size() == 0 && vec_y.size() != 0){
            x=0;
            for (int i_y=0; i_y<vec_y.size(); i_y++) {
                y = vec_y[i_y];
                h2D->Fill(x,y);
            }
        }
        else if(vec_x.size() != 0 && vec_y.size() == 0){
            y=0;
            for (int i_x=0; i_x<vec_x.size(); i_x++) {
                x = vec_x[i_x];
                h2D->Fill(x,y);
            }
        }
        else {
            x=0;
            y=0;
            h2D->Fill(x,y);
        }
    }
}

map< int, vector<float> > getHits(string str){
    map< int, vector<float> > hit;
    
    double newpos;
    for (int i=0; i<NLayer; i++) {
        vector<float> bla;
        bla.clear();
        hit[i] = bla;
    }
    
    for (unsigned int ihit=0; ihit<layerID->size(); ihit++) {
        if (str == string("hitpos"))
            newpos = hit_position->at(ihit) - alignmentpar[layerID->at(ihit)];
        else if (str == string("hitamp"))
            newpos = hit_amplitude->at(ihit);
        
        hit[layerID->at(ihit)].push_back(newpos);
    }
    
    return hit;
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
    
    // match pointnumber with layers
    // it is ordered as xlayers[i]
    //      so p[0] is ref, p[1] is A, p[2] is B
    
    double p_exp[3];
    p_exp[2] = p[2][2];
    p_exp[0] = getExpectedHit(p[0][0],p[0][2],p[1][0],p[1][2],p_exp[2]);
    p_exp[1] = getExpectedHit(p[0][1],p[0][2],p[1][1],p[1][2],p_exp[2]);
    
    
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
        h1D->Fill(p_exp[0] - p[i][0]);
        
        histname = spatialRes_histname + IDlayermap[ylayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(p_exp[1] - p[i][1]);
        
    }
    
    for (unsigned int i=0; i<3; i++) {
        layername = IDlayermap[xlayers[i]];
        TString tempname = layername(0, 1);
        histname = hitmap_histname + tempname;
        h2D = dynamic_cast<TH2D*> (rootobjects[histname]);
        h2D->Fill(p[i][0], p[i][1]);
    }
    
    // Angular resolution
    
    
    histname = angularRes_histname + TString("B");
    h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
    // angle between B_meas - A_meas - B_exp
    h1D->Fill(getAngleABC(p[2],p[1],p_exp));
    
    
}


void fillModifiedHistos(){
    
    TString histname,layername;
    TH1D *h1D;
    TH1D *h1D_mod;
    
    double newpos;
    // [pointnumber][x,y,z]
    double p[3][3];
    
    for (unsigned int ihit=0; ihit<layerID->size(); ihit++) {
        
        layername = IDlayermap[layerID->at(ihit)];
        newpos = hit_position->at(ihit) - alignmentpar[layerID->at(ihit)];
        
        for (unsigned int i=0; i<3; i++) {
            if (layerID->at(ihit) == xlayers[i]) {
                p[i][0] = newpos;
                p[i][2] = layerZposition[layerID->at(ihit)];
            }
            else if (layerID->at(ihit) == ylayers[i])
                p[i][1] = newpos;
        }
        
    }
    
    // match pointnumber with layers
    // it is ordered as xlayers[i]
    //      so p[0] is R, p[1] is A, p[2] is B
    
    double p_exp[3];
    p_exp[2] = p[2][2];
    p_exp[0] = getExpectedHit(p[0][0],p[0][2],p[1][0],p[1][2],p_exp[2]);
    p_exp[1] = getExpectedHit(p[0][1],p[0][2],p[1][1],p[1][2],p_exp[2]);
    
    
    for (unsigned int i=1; i<3; i++) {
        
        // spatial resolution
        histname = spatialRes_histname + IDlayermap[xlayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        histname = spatialRes_histname + IDlayermap[xlayers[i]]+ TString("_mod");
        h1D_mod = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D_mod->Fill(p_exp[0] - p[i][0]-h1D->GetMean());
        
        histname = spatialRes_histname + IDlayermap[ylayers[i]];
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        histname = spatialRes_histname + IDlayermap[ylayers[i]] + TString("_mod");
        h1D_mod = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D_mod->Fill(p_exp[1] - p[i][1]-h1D->GetMean());
        
    }
    
    // Angular resolution
    
    
    histname = angularRes_histname + TString("B");
    h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
    histname = angularRes_histname + TString("B_mod");
    h1D_mod = dynamic_cast<TH1D*> (rootobjects[histname]);
    // angle between down_meas - up_meas - down_exp
    h1D_mod->Fill(getAngleABC(p[2],p[1],p_exp)-h1D->GetMean());
    
    
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
    // here
    outfile.open(output_efficiency_txtfile.Data(), std::ofstream::out | std::ofstream::app);
    
    outfile << runnum << ";";
    for (map<unsigned short,TString>::iterator it=IDlayermap.begin(); it != IDlayermap.end(); it++) {

        double eff = double (layercount[it->first])/double (expectedcount[it->first]);
        double eff_err = sqrt(eff*(1-eff)/double (expectedcount[it->first]));
        cout<<"Layer: "<<IDlayermap[it->first]<<" efficiency: "<<eff<<" +/- "<< eff_err<<"       found hits: "<<layercount[it->first]<<" expected hits: "<<expectedcount[it->first]<<endl;
        outfile<< eff << ";"<<eff_err<<";";
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
    //LayerID; LayerName; AlignmentShift;
    
    while ( getline(alignmentFile,alignmentline) ){
        vector<string> elems = splitstring(alignmentline,';');
        // layer ID
        a_layerID = stoi( elems[0].c_str() );
        // layer name
        a_layername = string (elems[1]);
        // alignment parameter
        a_alignmentpar = atof(elems[2].c_str());
        alignmentpar[a_layerID] = a_alignmentpar;
    }
    
}

void printResolution(TString runnum){
    TString layername, histname, fitname;
    TF1* fitfunc1;
    TF1* fitfunc2;
    TH1D* h1D;
    double max_inxaxis;
    double mean,rms;
    
    std::ofstream outfile;
    
    outfile.open(output_resolution_txtfile.Data(), std::ofstream::out | std::ofstream::app);
    outfile << runnum << ";";
    
    for (map<unsigned short,TString>::iterator it=IDlayermap.begin(); it != IDlayermap.end(); it++) {
        
        // point resolution
        
        layername =it->second;
        
        histname =  spatialRes_histname + layername + TString("_mod");
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        mean = h1D->GetMean();
        rms = h1D->GetRMS();
        
        max_inxaxis =  h1D->GetXaxis()->GetBinCenter(h1D->GetMaximumBin());
        
        fitname = TString("fitgaus_") + histname;
        fitfunc1 = new TF1(fitname.Data(),"gaus", mean-3*rms,mean+3*rms);
        h1D->Fit(fitfunc1,"QR");
        
        fitname = TString("fitgauspol_") + histname;
        fitfunc2 = new TF1(fitname.Data(),"gaus+[3]", max_inxaxis-0.4,max_inxaxis+0.4);
        fitfunc2->SetLineColor(kGreen);
        fitfunc2->SetParameters(fitfunc1->GetParameter(0),fitfunc1->GetParameter(1),fitfunc1->GetParameter(2),0);
        fitfunc2->SetParLimits(3,0,10000);
        
        h1D->Fit(fitfunc2,"+QR");
        double spat_res =fitfunc2->GetParameter(2);
        double stat_err =fitfunc2->GetParError(2);
        double sys_err = fabs(stat_err - fitfunc1->GetParError(2));
        double spat_res_err = stat_err+sys_err;
        
        cout<< "spatial resolution of "<< layername<< " is "<< spat_res <<" +/- "<<spat_res_err<< " (mm)      "<<stat_err<<"(stat) + "<<sys_err<<"(sys)" << endl;
        outfile<<spat_res<<";"<<spat_res_err<<";";
    }
    
    
    histname = angularRes_histname + TString("B_mod");
    h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
    mean = h1D->GetMean();
    rms = h1D->GetRMS();
    max_inxaxis =  h1D->GetXaxis()->GetBinCenter(h1D->GetMaximumBin());
    
    fitname = TString("fitgaus_") + histname;
    fitfunc1 = new TF1(fitname.Data(),"gaus", mean-3*rms,mean+3*rms);
    h1D->Fit(fitfunc1,"R");
    
    fitname = TString("fitgauspol_") + histname;
    fitfunc2 = new TF1(fitname.Data(),"gaus+[3]", max_inxaxis-0.4,max_inxaxis+0.4);
    fitfunc2->SetLineColor(kGreen);
    fitfunc2->SetParameters(fitfunc1->GetParameter(0),fitfunc1->GetParameter(1),fitfunc1->GetParameter(2),0);
    fitfunc2->SetParLimits(3,0,10000);
    h1D->Fit(fitfunc2,"+R");
    
    double ang_res =fitfunc2->GetParameter(2) ;
    double stat_err =fitfunc2->GetParError(2) ;
    double sys_err = fabs(stat_err - (fitfunc1->GetParError(2)));
    double ang_res_err = stat_err+sys_err;
    
    
    cout<< "angular resolution of B is "<< ang_res <<" +/- "<< ang_res_err<<" degrees      "<<stat_err<<"(stat) + "<<sys_err<<"(sys)" <<endl;
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
        histname = exphitpos_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Expected HitPos ")+ it->second + TString(";Hit Position (mm);Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 10000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        // hit amplitude
        histname = hitamp_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Hit amplitude of all clusters in ") + it->second + TString(";Hit Amplitude;Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 1000,0, 10000);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        
        histname = exphitmap_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("expHitMap ")+ it->second + TString(";X (mm);Y (mm)");
            h2D = new TH2D(histname.Data(), title.Data(), 1000,0, 100,1000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h2D));
        }
        
        
    }
    
    for (unsigned int i=1; i<3; i++) {
        // spatial resolution
        histname = spatialRes_histname + IDlayermap[xlayers[i]];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Spatial Resolution;")+IDlayermap[xlayers[i]]+TString("_{expected} - ")+IDlayermap[xlayers[i]]+TString("_{measured} (mm);Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 10000,-100, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        histname = spatialRes_histname + IDlayermap[ylayers[i]];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Spatial Resolution;")+IDlayermap[ylayers[i]]+TString("_{expected} - ")+IDlayermap[ylayers[i]]+TString("_{measured} (mm);Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 10000, -100, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        // modified spatial resolution
        histname = spatialRes_histname + IDlayermap[xlayers[i]]+TString("_mod");
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Spatial Resolution;")+IDlayermap[xlayers[i]]+TString("_{expected} - ")+IDlayermap[xlayers[i]]+TString("_{measured} (mm);Number of entries");
            h1D = new TH1D(histname.Data(), title.Data(), 10000,-100, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h1D));
        }
        
        histname = spatialRes_histname + IDlayermap[ylayers[i]]+TString("_mod");
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Spatial Resolution;")+IDlayermap[ylayers[i]]+TString("_{expected} - ")+IDlayermap[ylayers[i]]+TString("_{measured} (mm);Number of entries");
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
    histname = angularRes_histname + TString("B");
    if (rootobjects.find(histname) == rootobjects.end()) {
        title = TString("Angular Resolution (B); Angle difference (degrees);Number of entries");
        h1D = new TH1D(histname.Data(), title.Data(), 1000,-30, 30);
        rootobjects.insert(pair<TString,TObject*>(histname,h1D));
    }
    // angular resolution
    histname = angularRes_histname + TString("B_mod");
    if (rootobjects.find(histname) == rootobjects.end()) {
        title = TString("Angular Resolution (B); Angle difference (degrees);Number of entries");
        h1D = new TH1D(histname.Data(), title.Data(), 1000,-30, 30);
        rootobjects.insert(pair<TString,TObject*>(histname,h1D));
    }
    
    TString s[3] = {"B","A","R"};
    for (unsigned int i=0; i<3; i++) {
        histname = hitmap_histname + s[i];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("HitMap ")+ s[i]+TString(";X (mm);Y (mm)");
            h2D = new TH2D(histname.Data(), title.Data(), 1000,0, 100,1000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h2D));
        }
        histname = hitamp2D_histname + s[i];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Hit Amplitude 2D ")+ s[i] + TString(";X (mm);Y (mm)");
            h2D = new TH2D(histname.Data(), title.Data(), 1000,0, 10000,1000,0, 10000);
            rootobjects.insert(pair<TString,TObject*>(histname,h2D));
        }
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


void setBranches(TTree* input_tree){
    
    input_tree->SetBranchAddress("event_num",&event_num);
    input_tree->SetBranchAddress("layerID",&layerID);
    input_tree->SetBranchAddress("cluster_size",&cluster_size);
    input_tree->SetBranchAddress("hit_amplitude",&hit_amplitude);
    input_tree->SetBranchAddress("hit_position",&hit_position);
}




