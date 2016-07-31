 


{
TCanvas *cc = new TCanvas("cc","",800,600);
    cc->SetRightMargin(0.125);
    cc->SetLeftMargin(0.125);
    cc->SetBottomMargin(0.13);
    cc->SetTopMargin(0.07);

TFile *f = TFile::Open("output/hitmaker_134.root");

TString ss;
TH2D *htime;
for(int i=0; i<2; i++){
if (i==0){
	htime = (TH2D*)f->Get("event/event_13DownX_timeVSchannel");
	htime->SetAxisRange(150,250,"X");
}
if (i==1){
	htime = (TH2D*)f->Get("event/event_13DownY_timeVSchannel");
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

}







}
