/*
 * Created by Eda Yildirim
 *  (2016 JGU)
 *
 *  email:eda.yildirim@cern.ch
 */

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <fstream>
#include <stdlib.h>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObject.h"
#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"

#include "Settings.h"
#include "HitMaker.h"


map<TString,TObject*> rootobjects;

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


TString hitmap_histname = TString("hitmap_");
TString signal_histname = TString("signal_");
TString output_tree_name = TString("hitsTree");
TString NclusterPerEvent = TString("nClusterPerEvent_");
TString NhitPerEvent = TString("nHitPerEvent_");
TString hitamp_histname = TString("hitamplitude_");
TString seedhitamp_histname = TString("seedchanamplitude_");
TString hitamp2D_histname = TString("hitamplitude2D_");
TString clusize_histname = TString("clustersize_");



std::map<int,TString> apvIDmap;
map<TString, float> HitAmplitudeCut;

int main(int argc, char *argv[]){
    
    TString runnum = TString(argv[1]);
    
    if (TString(argv[2]) == TString("desyTBUntil11_05_2016_15_30.config")) {
        apvIDmap = apvIDmap1;
    }
    else if (TString(argv[2]) == TString("desyTB.config")) {
        apvIDmap = apvIDmap2;
    }
    else {
        cout<<"APV ID Map cannot be found!!!"<<endl;
    }
    
    readHitAmplitudeCuts(runnum);
    
    TString fin_name = inputPath + TString("run")+ runnum + TString(".root");
    TFile * fin = TFile::Open(fin_name.Data());
    
    TString fout_name = outputPath + TString("hitmaker_")+runnum+TString(".root");
    TFile *fout = new TFile(fout_name.Data(),"RECREATE");
    
    unsigned int tevent_num;
    vector<unsigned short> tlayerID;
    vector<unsigned short> tcluster_size;
    vector<float> thit_amplitude;
    vector<float> thit_position;
    
    
    //create output tree
    TTree *output_tree = new TTree(output_tree_name.Data(),"tree of hits");
    output_tree->Branch("event_num",&tevent_num);
    output_tree->Branch("layerID",&tlayerID);
    output_tree->Branch("cluster_size",&tcluster_size);
    output_tree->Branch("hit_amplitude",&thit_amplitude);
    output_tree->Branch("hit_position",&thit_position);
    
    rootobjects.insert(pair<TString,TObject*>(output_tree_name,output_tree));
    
    
    createHitMapHistos();
    createOtherHistos();
    
    TTree* rawTree = (TTree*)fin->Get(inputRawTreeName.Data());
    TTree* dataTree = (TTree*)fin->Get(inputDataTreeName.Data());
    setBranches(rawTree,dataTree);
    TH1D *h,*h1D;
    
    TString histname;
    Long64_t nentries = rawTree->GetEntries();
    
    if (rawTree->GetEntries() != dataTree->GetEntries()) {
        cout << "raw and data trees don't have same number of events"<< endl;
    }
    int iplotgoodevent=0;
    int iplotbadevent=0;
    int iplotnothitevent=0;
    
    // store all clusters in an event to a vector
    unsigned int count=0;
    for (Long64_t ientry=0; ientry<nentries;ientry++) {
        rawTree->GetEntry(ientry);
        dataTree->GetEntry(ientry);
        
        tevent_num = apv_evt;
        tlayerID.clear();
        tcluster_size.clear();
        thit_amplitude.clear();
        thit_position.clear();
        
        if (apv_evt > MaxEventNum) break;
        
        if (apv_evt < MAXEVENTTOPLOT) plotEvent("event");
        vector<Cluster> clusters;
        clusters.clear();
        // the mm_strip and its apv_qmax paired and stored in a map for each layer... and mapped with layer name
        map < TString, map<unsigned int, short> > layers;
        for (int i = 0; i < apv_id->size(); i++) {
            if (apv_qmax->at(i) < SignalCut) continue;
            
            TString layername = apvIDmap[apv_id->at(i)];
            
            // check if this layer is stored in layers map
            if ( layers.find(layername) == layers.end() ) {
                //cout<< "Creating Layer "<<layername<< endl;
                map<unsigned int, short> emptyMap;
                
                layers.insert(make_pair(layername, emptyMap) );
            }
            // well now the layer is in the layers map
            layers[layername].insert( make_pair(mm_strip->at(i), apv_qmax->at(i)) );
        }
        // now clustering
        for (map< TString, map<unsigned int, short> >::iterator i_layer=layers.begin(); i_layer != layers.end(); i_layer++) {
            
            // since maps are great, and they are already sorted according to the key which is mm_strip
            double tot_sig = 0;
            double weight = 0;
            unsigned int clusize = 0;
            unsigned int prev_chan = 0;
            double seedchanAmp = 0;
            unsigned int clu_count=0;
            
            i_layer->second.insert( make_pair(1000000, 0));
            
            for (map<unsigned int, short>::iterator i_chan=i_layer->second.begin(); i_chan != i_layer->second.end(); i_chan++) {
                
                // start new cluster
                if (clusize==0) {
                    tot_sig = double(i_chan->second);
                    weight = double(i_chan->first) * double(i_chan->second);
                    clusize = 1;
                    prev_chan = i_chan->first;
                    seedchanAmp =double(i_chan->second);
                }
                else{
                    if (i_chan->first == prev_chan+1 || i_chan->first == prev_chan+2) {
                        tot_sig += double (i_chan->second);
                        weight += double (i_chan->first * i_chan->second);
                        clusize++;
                        if (double (i_chan->second) > seedchanAmp) {
                            seedchanAmp =double (i_chan->second);
                        }
                        prev_chan = i_chan->first;
                        
                    } else {
                        TString layerdirection = i_layer->first;
                        layerdirection = layerdirection(layerdirection.Length() -1, layerdirection.Length());
                        // process and save the cluster
                        if (  (layerdirection == TString("X") && clusize >CluSizeCutX) || (layerdirection == TString("Y") && clusize >CluSizeCutY)  )
                            if (tot_sig > HitAmplitudeCut[i_layer->first]){
                                Cluster acluster;
                                if (layerIDmap[i_layer->first]==2) {
                                    acluster._layerID = 3;
                                    tlayerID.push_back(3);
                                }
                                else if(layerIDmap[i_layer->first]==3){
                                    acluster._layerID = 2;
                                    tlayerID.push_back(2);
                                }
                                else {
                                    acluster._layerID = layerIDmap[i_layer->first];
                                    tlayerID.push_back(layerIDmap[i_layer->first]);
                                }
                                acluster._cluster_size = clusize;
                                acluster._hit_amplitude = tot_sig;
                                acluster._hit_position = (weight/tot_sig) * pitchsize;
                                
                                
                                histname = hitamp_histname + i_layer->first;
                                h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
                                h1D->Fill(acluster._hit_amplitude);
                                
                                histname = seedhitamp_histname + i_layer->first;
                                h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
                                h1D->Fill(seedchanAmp);

                                histname = clusize_histname + i_layer->first;
                                h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
                                h1D->Fill(clusize);
                                
                                tcluster_size.push_back(clusize) ;
                                thit_amplitude.push_back(tot_sig);
                                thit_position.push_back((weight/tot_sig) * pitchsize);
                                
                                clusters.push_back(acluster);
                                clu_count++;
                                
                            }
                        
                        
                        tot_sig=0;
                        weight=0;
                        clusize=0;
                        
                        tot_sig = double(i_chan->second);
                        weight = double(i_chan->first) * double(i_chan->second);
                        clusize = 1;
                        prev_chan = i_chan->first;
                    seedchanAmp =double(i_chan->second);
                    }
                    
                    
                }
                //                cout<<endl;
                
                
            }// end of iteration over channels
            
            histname = NclusterPerEvent + i_layer->first;
            h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
            h1D->Fill(clu_count);
            //            cout<<"layer "<<i_layer->first<< " ncluster "<<clu_layer.size()<<endl;
            
        } // end of iteration over layers
        
        output_tree->Fill();
        // check if the event is good
        // event has to have only 1 hit in each layer!
        
        histname = NhitPerEvent;
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(clusters.size());
        
        
        
        if (isGoodEvent(clusters) == true) {
            count++;
            processEvent(clusters);
            if(iplotgoodevent<MAXEVENTTOPLOT && apv_evt > MAXEVENTTOPLOT ){
                plotEvent("good");
                iplotgoodevent++;
            }
        }
        else if (isNoHitOnLayer(clusters) == true){
            //processEvent(clusters);
            if(iplotnothitevent<MAXEVENTTOPLOT && apv_evt > MAXEVENTTOPLOT ){
                plotEvent("nohit");
                iplotnothitevent++;
            }

        }
        else{
            if(iplotbadevent<MAXEVENTTOPLOT && apv_evt > MAXEVENTTOPLOT ){
                plotEvent("bad");
                iplotbadevent++;
            }
        }
        
        
    } // end of loop over events
    
    
    
    createOutputFile(fout);
    
    
    cout<<"Number of good events: "<<count <<endl;
    cout<<"Total number of events: "<<nentries<<endl;
    return 0;
}






bool isGoodEvent(vector<Cluster> clus){
    
    short nhits[NLayer]={0};
    bool onlyonehit = true;
    
    // also check if every layer has only 1 cluster
    for (unsigned int i=0; i<clus.size(); i++) {
        nhits[clus[i]._layerID]++;
    }
    // well now all values should be 1
    for (unsigned int i=0; i<NLayer; i++) {
        if (nhits[i]!=1) onlyonehit = false;
    }
    
    return onlyonehit;
}

bool isNoHitOnLayer(vector<Cluster> clus){
    
    short nhits[NLayer]={0};
    bool nohit = false;
    
    // also check if every layer has only 1 cluster
    for (unsigned int i=0; i<clus.size(); i++) {
        nhits[clus[i]._layerID]++;
    }
    // well now all values should be 1
    for (unsigned int i=0; i<NLayer; i++) {
        if (nhits[i]==0) nohit = true;
    }
    
    return nohit;
}


void processEvent(vector<Cluster> clusters){
    TString histname,layername;
    TH1D *h1D;
    TH2D* h2D;
    
    Point p[3];
    
    for (unsigned int icluster=0; icluster<clusters.size(); icluster++) {
        
        for (unsigned int i=0; i<3; i++) {
            if (clusters[icluster]._layerID == xlayers[i]) {
                p[i].x = clusters[icluster]._hit_position;
                p[i].z = zpos[i];
            }
            else if (clusters[icluster]._layerID == ylayers[i])
                p[i].y = clusters[icluster]._hit_position;
        }
    }
    
    for (unsigned int icluster=0; icluster<clusters.size(); icluster++) {
        layername = IDlayermap[clusters[icluster]._layerID];
        
        histname = signal_histname + layername;
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(clusters[icluster]._hit_amplitude);
        
        
        histname = hitmap_histname + layername;
        h1D = dynamic_cast<TH1D*> (rootobjects[histname]);
        h1D->Fill(clusters[icluster]._hit_position);
        
    }
    
}

void createOtherHistos(){
    
    TString histname,title;
    TH1D* h;
    TH2D* h2D;
    
    histname = NhitPerEvent;
    if (rootobjects.find(histname) == rootobjects.end()) {
        title = TString("Number of hits in an event;nHits;Number of entries");
        h = new TH1D(histname.Data(), title.Data(), 50,0, 50);
        rootobjects.insert(pair<TString,TObject*>(histname,h));
    }
    
    for (map< TString, unsigned short>::iterator i_layer=layerIDmap.begin(); i_layer != layerIDmap.end(); i_layer++) {
        
        histname = NclusterPerEvent + i_layer->first;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Number of clusters per event in ")+i_layer->first+TString(";Ncluster/event;Number of entries");
            h = new TH1D(histname.Data(), title.Data(), 20,0, 20);
            rootobjects.insert(pair<TString,TObject*>(histname,h));
        }
        
        histname = hitamp_histname + i_layer->first;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Hit amplitude of all clusters in ")+i_layer->first+TString(";Hit Amplitude;Number of entries");
            h = new TH1D(histname.Data(), title.Data(), 1000,0, 10000);
            rootobjects.insert(pair<TString,TObject*>(histname,h));
        }
        histname = seedhitamp_histname + i_layer->first;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("seed chan amplitude of all clusters in ")+i_layer->first+TString(";Hit Amplitude;Number of entries");
            h = new TH1D(histname.Data(), title.Data(), 1000,0, 10000);
            rootobjects.insert(pair<TString,TObject*>(histname,h));
        }

        histname = clusize_histname + i_layer->first;
        if (rootobjects.find(histname) == rootobjects.end()) {
            title = TString("Cluster Size ")+i_layer->first+TString(";Cluster Size;Number of entries");
            h = new TH1D(histname.Data(), title.Data(), 100,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h));
        }
        
    }
    
    
    
}


void plotEvent(TString prefix){
    TH1D* hmaxsignalVSchannel;
    TH1D* hhighChan_signalVStime;
    TH2D* htimeVSchannel;
    
    TString event_name = TString("event_XXX");
    event_name = prefix + event_name;
    event_name.ReplaceAll("XXX",TString::Itoa(apv_evt,10));
    TString histname,title;
    
    
    // create hists
    for (map<int,TString>::iterator it=apvIDmap.begin(); it != apvIDmap.end(); it++) {
        
        // time vs strip number
        histname = event_name + it->second + TString("_timeVSchannel");
        title = histname+TString(";strip number;time(ns)");
        
        if (rootobjects.find(histname) == rootobjects.end()) {
            htimeVSchannel = new TH2D(histname.Data(), title.Data(), 360, 0, 360, 27, 0, 675);
            rootobjects.insert(pair<TString,TObject*>(histname,htimeVSchannel));
        }
        
        // max signal vs strip number
        histname = event_name + it->second + TString("_maxsignalVSchannel");
        title = histname+TString(";max signal;strip number");
        if (rootobjects.find(histname) == rootobjects.end()) {
            hmaxsignalVSchannel = new TH1D(histname.Data(), title.Data(), 360, 0, 360);
            rootobjects.insert(pair<TString,TObject*>(histname,hmaxsignalVSchannel));
        }
        
        // signal vs time for highest signal strip
  /*
        histname = event_name + it->second + TString("_highChan_signalVStime");
        title = histname+TString(";signal;time (ns)");
        if (rootobjects.find(histname) == rootobjects.end()) {
            
            hhighChan_signalVStime = new TH1D(histname.Data(), title.Data(), 27, 0, 675);
            rootobjects.insert(pair<TString,TObject*>(histname,hhighChan_signalVStime));
        }
  */
        
    }
    
    
    short time_ns;
    vector<short> timedist;
    short max_signal;
    short max_max_signal=0;
    unsigned int highest_chan_entry;
    TString layername;
    
    //vector<vector<short> >& vecRef =apv_q;
    // now fill them
    for (unsigned int istrip=0; istrip<mm_strip->size(); istrip++) {
        layername = apvIDmap[apv_id->at(istrip)];
        
        timedist.clear();
        timedist = apv_q->at(istrip);
        max_signal = apv_qmax->at(istrip);
        
        if (max_signal> max_max_signal) {
            max_max_signal = max_signal;
            highest_chan_entry =istrip;
        }
        
        histname = event_name + layername + TString("_maxsignalVSchannel");
        hmaxsignalVSchannel = dynamic_cast<TH1D*> (rootobjects[histname]);
        hmaxsignalVSchannel->Fill(mm_strip->at(istrip),max_signal);
        
        for (unsigned int itime=0; itime<timedist.size(); itime++) {
            time_ns = itime*25;
            histname = event_name + layername + TString("_timeVSchannel");
            htimeVSchannel = dynamic_cast<TH2D*> (rootobjects[histname]);
            htimeVSchannel->Fill(mm_strip->at(istrip),time_ns,timedist[itime]);
        }
    }
    
    // for highest channel
    timedist.clear();
    layername = apvIDmap[apv_id->at(highest_chan_entry)];
    //    &timedist = apv_q->at(highest_chan_entry);
    /*

    for (unsigned int itime=0; itime<timedist.size(); itime++) {
        time_ns = itime*25;
        histname = event_name + layername + TString("_highChan_signalVStime");
        hhighChan_signalVStime = dynamic_cast<TH1D*> (rootobjects[histname]);
        // hhighChan_signalVStime->Fill(time_ns,vecRef[highest_chan_entry][itime]);
         
        hhighChan_signalVStime->Fill(time_ns,(apv_q->at(highest_chan_entry)).at(itime));
        // hhighChan_signalVStime->Fill(time_ns,vecRef[highest_chan_entry][itime]);
         
    }
    */
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


void setBranches(TTree* rawtree, TTree* datatree){
    
    rawtree->SetBranchAddress("apv_evt", &apv_evt);
    rawtree->SetBranchAddress("time_s", &time_s);
    rawtree->SetBranchAddress("time_us", &time_us);
    rawtree->SetBranchAddress("apv_fecNo", &apv_fecNo);
    rawtree->SetBranchAddress("apv_id", &apv_id);
    rawtree->SetBranchAddress("apv_ch", &apv_ch);
    rawtree->SetBranchAddress("mm_id", &mm_id);
    rawtree->SetBranchAddress("mm_readout", &mm_readout);
    rawtree->SetBranchAddress("mm_strip", &mm_strip);
    rawtree->SetBranchAddress("apv_q", &apv_q);
    rawtree->SetBranchAddress("apv_presamples", &apv_presamples);
    
    datatree->SetBranchAddress("apv_qmax", &apv_qmax);
    datatree->SetBranchAddress("apv_tbqmax", &apv_tbqmax);
}


void createHitMapHistos(){
    TString histname;
    TH1D* h;
    for (map<int,TString>::iterator it=apvIDmap.begin(); it != apvIDmap.end(); it++) {
        
        histname = signal_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            h = new TH1D(histname.Data(), histname.Data(), 1000,0, 10000);
            rootobjects.insert(pair<TString,TObject*>(histname,h));
        }
        
        histname = hitmap_histname + it->second;
        if (rootobjects.find(histname) == rootobjects.end()) {
            h = new TH1D(histname.Data(), histname.Data(), 1000,0, 100);
            rootobjects.insert(pair<TString,TObject*>(histname,h));
        }
    }
    
}


void readHitAmplitudeCuts(TString runnum){
    cout<< "Reading HitAmplitudeCuts from "<<hitAmplitudeCut_file<<endl;
    
    
    string line;
    string tmpstr;
    bool runnum_found = false;

    ifstream fileinfo_file(hitAmplitudeCut_file.Data());
    //      0           1          2    3       4   5       6   7   8   9   10      11  12
    // MMDataRunNum;BeamEnergy;Angle1;Angle2;DownX;DownY;RefX;RefY;UpX;UpY;config;map;Comment;
    // first line
    getline(fileinfo_file,line);
    vector<string> titles = splitstring(line,';');
    while(getline(fileinfo_file,line)){
        vector<string> elems = splitstring(line,';');
        
        if (elems[0] == runnum) {
            for (int i=4; i<10; i++) {
                if (elems.size() < 12){
                    cout << "There is something wrong with entry of run "<<runnum<<" in file "<<hitAmplitudeCut_file<<endl;
                    return;
                }
                
                float value = atof(elems[i].c_str());
                HitAmplitudeCut.insert(make_pair( titles[i],value ));
                runnum_found = true;

            }
        }
    }
    
    if (runnum_found == false) {
        cout<<" No entry found for runnum "<<runnum<<" in file "<<hitAmplitudeCut_file<<endl;
    }

}
/*

vector<string> splitstring(string s, char delim) {
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

*/
