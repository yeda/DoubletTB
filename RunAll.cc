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
#include "TObject.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "RunAll.h"


using namespace std;

const double NOTSET = 100000;

TString output_resolution_txtfile= TString("./resolution.txt");
TString output_efficiency_txtfile= TString("./efficiency.txt");

TString output_root_filename = TString("measurementResults.root");

TString eff_histoname= TString("eff_");
TString res_histoname= TString("res_");
TString angres_histoname= TString("angres");

vector<RunInfo> runlist;
map<int, Meas* > measurements;

map<TString,TObject*> rootobjects;

int main(int argc, char *argv[]){

    readRunList(runlistfile);

  createTxtFiles();
    readMeasEff();
    readMeasRes();
    
    TFile *fout = new TFile(output_root_filename.Data(),"RECREATE");

    makePlots();
    
    createOutputFile(fout);
    
    return 0;
}

void createHistos(){
    TString histname, title;
    TGraphErrors* gr;
    
    vector<TString> eff_histlist = {"DownX", "DownY", "UpX", "UpY", "RefX", "RefY"};
    vector<TString> res_histlist = {"UpX", "UpY", "DownX", "DownY", "Angular"};
    
    // eff histos
    for (unsigned int k=0; k<eff_histlist.size(); k++) {
        histname = eff_histoname + eff_histlist[k];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Efficiency ") + eff_histlist[k] + TString(";Efficiency;Tilt angle (degrees)");
            gr = new TGraphErrors(measurements.size());
            gr->SetName(histname.Data());
            gr->SetTitle(title.Data());
            rootobjects.insert(pair<TString,TObject*>(histname,gr));
        }
    }
    
    // res histos
    for (unsigned int k=0; k<res_histlist.size(); k++) {
        histname = res_histoname + res_histlist[k];
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Resolution ") + res_histlist[k] + TString(";Resolution (mm);Tilt angle (degrees)");
            gr = new TGraphErrors(measurements.size());
            gr->SetName(histname.Data());
            gr->SetTitle(title.Data());
            rootobjects.insert(pair<TString,TObject*>(histname,gr));
        }
    }

}

void makePlots(){
    TString histname, title;
    TGraphErrors* gr;
    
    createHistos();
    cout << "Here "<<endl;
    
    int i_point = 0;
    for (map< int, Meas* >::iterator i_meas=measurements.begin(); i_meas != measurements.end(); i_meas++) {
        Meas* ameas = i_meas->second;

        ///////////// Efficiency
        histname = eff_histoname + TString("DownX");
        gr = dynamic_cast<TGraphErrors*> (rootobjects[histname]);
        gr->SetPoint(i_point, ameas->get_Angle1(), ameas->get_eff_DownX());
        //gr->SetPointError(i_point ,0, 0);

        histname = eff_histoname + TString("DownY");
        gr = dynamic_cast<TGraphErrors*> (rootobjects[histname]);
        gr->SetPoint(i_point, ameas->get_Angle1(), ameas->get_eff_DownY());

        histname = eff_histoname + TString("UpX");
        gr = dynamic_cast<TGraphErrors*> (rootobjects[histname]);
        gr->SetPoint(i_point, ameas->get_Angle1(), ameas->get_eff_UpX());

        histname = eff_histoname + TString("UpY");
        gr = dynamic_cast<TGraphErrors*> (rootobjects[histname]);
        gr->SetPoint(i_point, ameas->get_Angle1(), ameas->get_eff_UpY());

        histname = eff_histoname + TString("RefX");
        gr = dynamic_cast<TGraphErrors*> (rootobjects[histname]);
        gr->SetPoint(i_point, ameas->get_Angle1(), ameas->get_eff_RefX());

        histname = eff_histoname + TString("RefY");
        gr = dynamic_cast<TGraphErrors*> (rootobjects[histname]);
        gr->SetPoint(i_point, ameas->get_Angle1(), ameas->get_eff_RefY());

        i_point++;
        
        
        
        
    }

    
    
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

void createTxtFiles(){
    
    TString runthis;
    
    for (unsigned int irun=0; irun<runlist.size(); irun++) {
        runlist[irun].print();
        runthis = TString("./hitmaker ") + TString::Itoa(runlist[irun].DataRun,10) + TString(" ") + runlist[irun].config;
        cout<<runthis<<endl;
        system(runthis.Data());
    }
    
    system("./alignment 134");
   
    runthis = TString("rm ")+ output_efficiency_txtfile + TString(" ")+ output_resolution_txtfile;
    system(runthis.Data());
    
    std::ofstream out_eff;
    out_eff.open(output_efficiency_txtfile.Data(), std::ofstream::out | std::ofstream::app);
    /*
     std::map<unsigned short,TString> IDlayermap={
     {0,TString("DownX")},{1,TString("DownY")},
     {2,TString("UpX")},{3,TString("UpY")},
     {4,TString("RefX")},{5,TString("RefY")},
     };
     */
    out_eff<<"RunNum;DownX;DownY;UpX;UpY;RefX;RefY;"<<endl;
    out_eff.close();
    
    
    std::ofstream out_res;
    out_res.open(output_resolution_txtfile.Data(), std::ofstream::out | std::ofstream::app);
    out_res<<"RunNum;UpX;UpXerr;UpY;UpYerr;DownX;DownXerr;DownY;DownYerr;Angular;Angularerr;"<<endl;
    out_res.close();
    
    for (unsigned int irun=0; irun<runlist.size(); irun++) {
        runlist[irun].print();
        runthis = TString("./analyze ") + TString::Itoa(runlist[irun].DataRun,10) + TString(" 134");
        cout<<runthis<<endl;
        system(runthis.Data());
    }
    

}


void readMeasRes(){
    cout<< "Reading "<<output_resolution_txtfile<<endl;
    ifstream fileinfo_file(output_resolution_txtfile.Data());
    string firstline;
    getline(fileinfo_file,firstline);
    
    Meas* ameas;
    int runnum;
    
    string s;
    while(!fileinfo_file.eof()){
        getline(fileinfo_file,s,';');
        runnum = atoi(s.c_str());
        
        if (runnum == 0) {
            getline(fileinfo_file,s,';');
            //ameas->set_res_UpX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            //ameas->set_res_UpX_err(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            //ameas->set_res_UpY(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            //ameas->set_res_UpY_err(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            //ameas->set_res_DownX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            //ameas->set_res_DownX_err(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            //ameas->set_res_DownY(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            //ameas->set_res_DownY_err(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            //ameas->set_ang_res(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            //ameas->set_ang_res_err(atof(s.c_str()));
        }
        else {
            ameas = measurements[runnum];
            
            getline(fileinfo_file,s,';');
            ameas->set_res_UpX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            ameas->set_res_UpX_err(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            ameas->set_res_UpY(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            ameas->set_res_UpY_err(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            ameas->set_res_DownX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            ameas->set_res_DownX_err(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            ameas->set_res_DownY(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            ameas->set_res_DownY_err(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            ameas->set_ang_res(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            ameas->set_ang_res_err(atof(s.c_str()));
        }
    }
    
}

void readMeasEff(){
    cout<< "Reading "<<output_efficiency_txtfile<<endl;
    
    ifstream fileinfo_file(output_efficiency_txtfile.Data());
    string firstline;
    getline(fileinfo_file,firstline);
    
    Meas* ameas;
    int runnum;
    double temp_d[2];
    string s;
    while(!fileinfo_file.eof()){
        getline(fileinfo_file,s,';');
        runnum = atoi(s.c_str());
        
        if (runnum == 0) {
            getline(fileinfo_file,s,';');
            //ameas->set_eff_DownX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            //ameas->set_eff_DownY(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            //ameas->set_eff_UpX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            //ameas->set_eff_UpY(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            //ameas->set_eff_RefX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            //ameas->set_eff_RefY(atof(s.c_str()));
        }
        else {
            ameas = dynamic_cast<Meas*> (measurements[runnum]);
            
            getline(fileinfo_file,s,';');
            ameas->set_eff_DownX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            ameas->set_eff_DownY(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            ameas->set_eff_UpX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            ameas->set_eff_UpY(atof(s.c_str()));
            
            getline(fileinfo_file,s,';');
            ameas->set_eff_RefX(atof(s.c_str()));
            getline(fileinfo_file,s,';');
            ameas->set_eff_RefY(atof(s.c_str()));
            
        }
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

