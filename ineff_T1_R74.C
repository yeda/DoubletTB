{

TFile *_file74 = TFile::Open("output/analyze_74.root");
TH2D* h74_AL2 = (TH2D*)_file74->Get("exphitmap_A-L2");
TCanvas* cc74_AL2 = new TCanvas();
cc74_AL2->SetLogz(1);
h74_AL2->SetTitle("tilt 1 A-L2");
h74_AL2->SetAxisRange(5,25,"Y");
h74_AL2->SetAxisRange(55,75,"X");
h74_AL2->SetAxisRange(0,10,"Z");
h74_AL2->Draw("colz");

TH2D* h74_AL1 = (TH2D*)_file74->Get("exphitmap_A-L1");
TCanvas* cc74_AL1 = new TCanvas();
cc74_AL1->SetLogz(1);
h74_AL1->SetTitle("tilt 1 A-L1");
h74_AL1->SetAxisRange(5,25,"Y");
h74_AL1->SetAxisRange(55,75,"X");
h74_AL1->SetAxisRange(0,10,"Z");
h74_AL1->Draw("colz");

TH2D* h74_BL2 = (TH2D*)_file74->Get("exphitmap_B-L2");
TCanvas* cc74_BL2 = new TCanvas();
cc74_BL2->SetLogz(1);
h74_BL2->SetTitle("tilt 1 B-L2");
h74_BL2->SetAxisRange(5,25,"Y");
h74_BL2->SetAxisRange(55,75,"X");
h74_BL2->SetAxisRange(0,10,"Z");
h74_BL2->Draw("colz");

TH2D* h74_BL1 = (TH2D*)_file74->Get("exphitmap_B-L1");
TCanvas* cc74_BL1 = new TCanvas();
cc74_BL1->SetLogz(1);
h74_BL1->SetTitle("tilt 1 B-L1");
h74_BL1->SetAxisRange(5,25,"Y");
h74_BL1->SetAxisRange(55,75,"X");
h74_BL1->SetAxisRange(0,10,"Z");
h74_BL1->Draw("colz");

cc74_AL1->Print("ineff_T1_R74.pdf(");
cc74_AL2->Print("ineff_T1_R74.pdf");
cc74_BL1->Print("ineff_T1_R74.pdf");
cc74_BL2->Print("ineff_T1_R74.pdf)");

}
