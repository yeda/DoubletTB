


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
    
    for(int i=0; i<2; i++){
        if (i==0){
            htime = (TH2D*)f->Get("goodevent/goodevent_13DownX_timeVSchannel");
            htime->SetAxisRange(150,250,"X");
        }
        if (i==1){
            htime = (TH2D*)f->Get("goodevent/goodevent_13DownY_timeVSchannel");
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
    
    
    // residual plots
    TLatex *latex =new TLatex();
    latex->SetNDC();
    latex->SetTextFont(43);
    latex->SetTextColor(1);
    latex->SetTextSize(26);
    
    
    
    TFile *_file0 = TFile::Open("output/analyze_84.root");
    TH1D * h1D;
    TF1* fgaus;
    TF1* fpol;
    ostringstream sss;
    
    for (int i=0; i<3; i++){
        if (i==0){
            h1D = (TH1D*)_file0->Get("spatialRes_DownX_mod");
            h1D->SetAxisRange(-1.3,1.7,"X");
            fgaus = (TF1*)h1D->GetFunction("fitgaus_spatialRes_DownX_mod");
            fpol = (TF1*)h1D->GetFunction("fitgauspol_spatialRes_DownX_mod");
            h1D->SetTitle(";X_{measured}-X_{expected} (mm);Number of Entries");
        }
        else if (i==1){
            h1D = (TH1D*)_file0->Get("spatialRes_DownY_mod");
            h1D->SetAxisRange(-1.3,1.7,"X");
            fgaus = (TF1*)h1D->GetFunction("fitgaus_spatialRes_DownY_mod");
            fpol = (TF1*)h1D->GetFunction("fitgauspol_spatialRes_DownY_mod");
            h1D->SetTitle(";Y_{measured}-Y_{expected} (mm);Number of Entries");
        }
        else{
            h1D = (TH1D*)_file0->Get("angularRes_Down_mod");
            h1D->SetAxisRange(-2.8,3.2,"X");
            fgaus = (TF1*)h1D->GetFunction("fitgaus_angularRes_Down_mod");
            fpol = (TF1*)h1D->GetFunction("fitgauspol_angularRes_Down_mod");
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
        
        double res =fpol->GetParameter(2) / sqrt(2);
        double stat_err =fpol->GetParError(2) / sqrt(2);
        double sys_err = fabs(stat_err - (fgaus->GetParError(2) / sqrt(2)));
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
    
}
