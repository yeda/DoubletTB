


{
#include <sstream>
    TGaxis::SetMaxDigits(3);
    
    TCanvas *cc = new TCanvas("cc","",800,600);
    cc->SetRightMargin(0.125);
    cc->SetLeftMargin(0.125);
    cc->SetBottomMargin(0.13);
    cc->SetTopMargin(0.07);
    
    TFile *f = TFile::Open("output/hitmaker_134.root");
    TString tmp_str;
    TString ss;
    TH2D *htime;
    
    // time vs channel
   /* 
    for(int i=0; i<2; i++){
        if (i==0){
            htime = (TH2D*)f->Get("goodevent/goodevent_13B-L2_timeVSchannel");
            htime->SetAxisRange(150,250,"X");
        }
        if (i==1){
            htime = (TH2D*)f->Get("goodevent/goodevent_13B-L1_timeVSchannel");
            htime->SetAxisRange(230,330,"X");
        }
        htime->GetXaxis()->SetLabelSize(0.045);
        htime->GetXaxis()->SetTitleSize(0.05);
        htime->GetXaxis()->SetTitleOffset(1.2);
        
        htime->GetYaxis()->SetLabelSize(0.045);
        htime->GetYaxis()->SetTitleSize(0.05);
        htime->GetYaxis()->SetTitleOffset(1.2);
        
        htime->SetTitle(";Strip Number;Time (ns)");
        htime->SetStats(0);
        
        cc->cd();
        htime->Draw("colz");
        ss = TString("results/") + TString(htime->GetName()) +  TString(".pdf");
        cc->SaveAs(ss.Data());
        ss = TString("results/") + TString(htime->GetName()) +  TString(".C");
        cc->SaveAs(ss.Data());
        
    }
    */
    
    // residual plots
    TLatex *latex =new TLatex();
    latex->SetNDC();
    latex->SetTextFont(43);
    latex->SetTextColor(1);
    latex->SetTextSize(26);
    
    
    
    TFile *_file0 = TFile::Open("output/selfaligned/analyze_110.root");
    TFile *_file1 = TFile::Open("output/analyze_110.root");
    TH1D * h1D;
    TF1* fgaus;
    TF1* fpol;
    ostringstream sss;
    
    for (int i=0; i<3; i++){
        if (i==0){
            h1D = (TH1D*)_file0->Get("spatialRes_A-L2_mod");
            h1D->SetAxisRange(-1.3,1.7,"X");
            fgaus = (TF1*)h1D->GetFunction("fitgaus_spatialRes_A-L2_mod");
            fpol = (TF1*)h1D->GetFunction("fitgauspol_spatialRes_A-L2_mod");
            h1D->SetTitle(";Position_{measured}-Position_{expected} (mm);Number of Entries");
        }
        else if (i==1){
            h1D = (TH1D*)_file0->Get("spatialRes_A-L1_mod");
            h1D->SetAxisRange(-1.3,1.7,"X");
            fgaus = (TF1*)h1D->GetFunction("fitgaus_spatialRes_A-L1_mod");
            fpol = (TF1*)h1D->GetFunction("fitgauspol_spatialRes_A-L1_mod");
            h1D->SetTitle(";Position_{measured}-Position_{expected} (mm);Number of Entries");
        }
        else{
            h1D = (TH1D*)_file1->Get("angularRes_B_mod");
            h1D->SetAxisRange(-2.8,3.2,"X");
            fgaus = (TF1*)h1D->GetFunction("fitgaus_angularRes_B_mod");
            fpol = (TF1*)h1D->GetFunction("fitgauspol_angularRes_B_mod");
            h1D->SetTitle(";#theta_{measured}-#theta_{expected} (degrees);Number of Entries");
        }
        h1D->GetXaxis()->SetLabelSize(0.045);
        h1D->GetXaxis()->SetTitleSize(0.05);
        h1D->GetXaxis()->SetTitleOffset(1.2);
        
        h1D->GetYaxis()->SetLabelSize(0.045);
        h1D->GetYaxis()->SetTitleSize(0.05);
        h1D->GetYaxis()->SetTitleOffset(1.2);
        
        h1D->SetLineWidth(1);
        h1D->SetLineColor(kBlack);
        h1D->SetStats(0);
        h1D->SetFillColor(kOrange);
        
        fpol->SetLineWidth(3);
        fpol->SetLineColor(kBlue);
        fgaus->SetLineWidth(3);
        fgaus->SetLineColor(kRed);
        
        double res =fpol->GetParameter(2);
        double stat_err =fpol->GetParError(2);
        double sys_err = fabs(stat_err - fgaus->GetParError(2));
        double res_err = stat_err+sys_err;
        cout << res <<" "<< stat_err <<" "<< sys_err <<" "<< res_err<<endl;
        sss.str("");
        sss<< string("#bf{#sigma = ")<< fixed<< setprecision(2)<<scientific<< res << string(" #pm ")<<fixed<<setprecision(2)<<scientific<<res_err<<string("}");
        tmp_str = TString(sss.str());
        
        TLegend *leg = new TLegend(0.55,0.65,0.9,0.92);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetLineColor(0);
        leg->SetTextSize(0.05);
        
        leg->AddEntry(h1D,"Data","f");
        leg->AddEntry(fgaus,"Gaussian","l");
        leg->AddEntry(fpol,"#splitline{Gaussian}{   + Constant}","l");
        
        cc->cd();
        h1D->Draw();
        leg->Draw();
        latex->SetTextColor(kOrange+7);
        latex->DrawLatex(0.54,0.5,tmp_str.Data());
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.17,0.8,"DESY Testbeam");
        latex->DrawLatex(0.17,0.72,"4.4 GeV   2016");
        
        ss = TString("results/") + TString(h1D->GetName()) +  TString(".pdf");
        cc->SaveAs(ss.Data());
        ss = TString("results/") + TString(h1D->GetName()) +  TString(".C");
        cc->SaveAs(ss.Data());
        
    }


 ////////////////   Hit Amp
 
TCanvas *ccc = new TCanvas("ccc","",800,600);
    ccc->SetRightMargin(0.125);
    ccc->SetLeftMargin(0.125);
    ccc->SetBottomMargin(0.13);
    ccc->SetTopMargin(0.07);
    


        TLegend *leg2 = new TLegend(0.55,0.65,0.9,0.92);
        leg2->SetFillColor(0);
        leg2->SetFillStyle(0);
        leg2->SetLineColor(0);
        leg2->SetTextSize(0.05);



  _file0 = TFile::Open("output/analyze_102.root");

TString legname;
Color_t col[4] = {kBlack, kRed, kBlue, kGreen+2};
TString histname;

for (int i=0; i<4; i++){
histname = TString("hitamplitude_");
if(i==0) {
	histname = histname + TString("B-L2");
	legname = TString("B - L2");

}
if(i==1){ 
	histname = histname + TString("B-L1");
	legname = TString("B - L1");
}
if(i==2){ 
	histname = histname + TString("A-L2");
	legname = TString("A - L1");
}
if(i==3){ 
	histname = histname + TString("A-L1");
	legname = TString("A - L2");
}
h1D = (TH1D*)_file0->Get(histname.Data());
h1D->SetTitle(";Total Cluster Charge (ADC);Number of Entries");
h1D->Rebin(10);
h1D->SetAxisRange(0,6000);

        h1D->GetXaxis()->SetLabelSize(0.045);
        h1D->GetXaxis()->SetTitleSize(0.05);
        h1D->GetXaxis()->SetTitleOffset(1.2);
        
        h1D->GetYaxis()->SetLabelSize(0.045);
        h1D->GetYaxis()->SetTitleSize(0.05);
        h1D->GetYaxis()->SetTitleOffset(1.2);
        
        h1D->SetLineWidth(2);
        h1D->SetStats(0);
	h1D->SetLineColor(col[i]);

	ccc->cd();
        h1D->Draw("same");
       
        leg2->AddEntry(h1D,legname.Data(),"l");
       
}	

  
        ccc->cd();
        leg2->Draw();
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.23,0.8,"DESY Testbeam");
        latex->DrawLatex(0.23,0.72,"4.4 GeV   2016");
 

        ss = TString("results/") + TString("HitAmpAll.pdf");
        ccc->SaveAs(ss.Data());
        ss = TString("results/") + TString("HitAmpAll.C");
        ccc->SaveAs(ss.Data());



  ////////////////   Clu size
 
TCanvas *cccc = new TCanvas("cccc","",800,600);
    cccc->SetRightMargin(0.125);
    cccc->SetLeftMargin(0.125);
    cccc->SetBottomMargin(0.13);
    cccc->SetTopMargin(0.07);
    


        TLegend *leg3 = new TLegend(0.6,0.65,0.95,0.92);
        leg3->SetFillColor(0);
        leg3->SetFillStyle(0);
        leg3->SetLineColor(0);
        leg3->SetTextSize(0.05);



  _file0 = TFile::Open("output/hitmaker_102.root");


for (int i=0; i<4; i++){
histname = TString("clustersize/clustersize_");
if(i==0) {
	histname = histname + TString("B-L2");
	legname = TString("B - L2");

}
if(i==1){ 
	histname = histname + TString("B-L1");
	legname = TString("B - L1");
}
if(i==2){ 
	histname = histname + TString("A-L1");
	legname = TString("A - L1");
}
if(i==3){ 
	histname = histname + TString("A-L2");
	legname = TString("A - L2");
}
h1D = (TH1D*)_file0->Get(histname.Data());
h1D->SetTitle(";Cluster Size;Number of Entries");
h1D->SetAxisRange(0,42000,"Y");
h1D->SetAxisRange(0,30,"X");

        h1D->GetXaxis()->SetLabelSize(0.045);
        h1D->GetXaxis()->SetTitleSize(0.05);
        h1D->GetXaxis()->SetTitleOffset(1.2);
        
        h1D->GetYaxis()->SetLabelSize(0.045);
        h1D->GetYaxis()->SetTitleSize(0.05);
        h1D->GetYaxis()->SetTitleOffset(1.2);
        
        h1D->SetLineWidth(2);
        h1D->SetStats(0);
	h1D->SetLineColor(col[i]);

	cccc->cd();
        h1D->Draw("same");
       
        leg3->AddEntry(h1D,legname.Data(),"l");
       
}	

  
        cccc->cd();
        leg3->Draw();
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.31,0.8,"DESY Testbeam");
        latex->DrawLatex(0.31,0.72,"4.4 GeV   2016");
 

        ss = TString("results/") + TString("CluSizeAll.pdf");
        cccc->SaveAs(ss.Data());
        ss = TString("results/") + TString("CluSizeAll.C");
        cccc->SaveAs(ss.Data());
 



}
