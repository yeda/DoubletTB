{

TFile *_file71 = TFile::Open("output/analyze_71.root");
TH2D* h71_AL2 = (TH2D*)_file71->Get("exphitmap_A-L2");
TCanvas* cc71_AL2 = new TCanvas();
cc71_AL2->SetLogz(1);
h71_AL2->SetTitle("tilt 0 A-L2");
h71_AL2->SetAxisRange(5,25,"Y");
h71_AL2->SetAxisRange(55,75,"X");
h71_AL2->SetAxisRange(0,10,"Z");
h71_AL2->Draw("colz");

TH2D* h71_AL1 = (TH2D*)_file71->Get("exphitmap_A-L1");
TCanvas* cc71_AL1 = new TCanvas();
cc71_AL1->SetLogz(1);
h71_AL1->SetTitle("tilt 0 A-L1");
h71_AL1->SetAxisRange(5,25,"Y");
h71_AL1->SetAxisRange(55,75,"X");
h71_AL1->SetAxisRange(0,10,"Z");
h71_AL1->Draw("colz");

TH2D* h71_BL2 = (TH2D*)_file71->Get("exphitmap_B-L2");
TCanvas* cc71_BL2 = new TCanvas();
cc71_BL2->SetLogz(1);
h71_BL2->SetTitle("tilt 0 B-L2");
h71_BL2->SetAxisRange(5,25,"Y");
h71_BL2->SetAxisRange(55,75,"X");
h71_BL2->SetAxisRange(0,10,"Z");
h71_BL2->Draw("colz");

TH2D* h71_BL1 = (TH2D*)_file71->Get("exphitmap_B-L1");
TCanvas* cc71_BL1 = new TCanvas();
cc71_BL1->SetLogz(1);
h71_BL1->SetTitle("tilt 0 B-L1");
h71_BL1->SetAxisRange(5,25,"Y");
h71_BL1->SetAxisRange(55,75,"X");
h71_BL1->SetAxisRange(0,10,"Z");
h71_BL1->Draw("colz");

cc71_AL1->Print("ineff_T0_R71.pdf(");
cc71_AL2->Print("ineff_T0_R71.pdf");
cc71_BL1->Print("ineff_T0_R71.pdf");
cc71_BL2->Print("ineff_T0_R71.pdf)");

}
