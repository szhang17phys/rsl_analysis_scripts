#include <iostream>
#include <TFile.h>
#include <TH2.h>


int singleOpCh_combine(){

    string rslxx ="rsl150";

    //For opch01=================================================
    TFile *file01_1 = new TFile("../results/tmp/opCh01_"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    if (!file01_1 || file01_1->IsZombie()) {
        std::cerr << "Error: Unable to open file1.root" << std::endl;
        return 1;
    }
    //Open the second root file---
    TFile *file01_2 = new TFile("../results/tmp/opCh01_"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file01_3 = new TFile("../results/tmp/opCh01_"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Get access to histos inside root files---
    TH2F *hist01_1 = dynamic_cast<TH2F*>(file01_1->Get("CRT_Opch")); 
    TH2F *hist01_2 = dynamic_cast<TH2F*>(file01_2->Get("CRT_Opch")); 
    TH2F *hist01_3 = dynamic_cast<TH2F*>(file01_3->Get("CRT_Opch")); 

    //Create new hist to store opch01-----------------------
    TH2F *hist_opch01 = new TH2F(*hist01_1);
    hist_opch01->Add(hist01_2);
    hist_opch01->Add(hist01_3);

    hist_opch01->SetName("opch01");
    hist_opch01->SetTitle("Response of OpCh01");
    hist_opch01->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch01->GetYaxis()->SetTitle("#photon / event");
    hist_opch01->SetMarkerStyle(21);
    hist_opch01->SetMarkerSize(1.0);




    //For opch03=================================================
    TFile *file03_1 = new TFile("../results/tmp/opCh03_"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    TFile *file03_2 = new TFile("../results/tmp/opCh03_"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file03_3 = new TFile("../results/tmp/opCh03_"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Get access to histos inside root files---
    TH2F *hist03_1 = dynamic_cast<TH2F*>(file03_1->Get("CRT_Opch")); 
    TH2F *hist03_2 = dynamic_cast<TH2F*>(file03_2->Get("CRT_Opch")); 
    TH2F *hist03_3 = dynamic_cast<TH2F*>(file03_3->Get("CRT_Opch")); 

    //Create new hist to store opch01---
    TH2F *hist_opch03 = new TH2F(*hist03_1);
    hist_opch03->Add(hist03_2);
    hist_opch03->Add(hist03_3);

    hist_opch03->SetName("opch03");
    hist_opch03->SetTitle("Response of OpCh03");
    hist_opch03->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch03->GetYaxis()->SetTitle("#photon / event");
    hist_opch03->SetMarkerStyle(21);
    hist_opch03->SetMarkerSize(1.0);



    //For opch16=================================================
    TFile *file16_1 = new TFile("../results/tmp/opCh16_"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    TFile *file16_2 = new TFile("../results/tmp/opCh16_"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file16_3 = new TFile("../results/tmp/opCh16_"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Get access to histos inside root files---
    TH2F *hist16_1 = dynamic_cast<TH2F*>(file16_1->Get("CRT_Opch")); 
    TH2F *hist16_2 = dynamic_cast<TH2F*>(file16_2->Get("CRT_Opch")); 
    TH2F *hist16_3 = dynamic_cast<TH2F*>(file16_3->Get("CRT_Opch")); 

    //Create new hist to store opch01---
    TH2F *hist_opch16 = new TH2F(*hist16_1);
    hist_opch16->Add(hist16_2);
    hist_opch16->Add(hist16_3);

    hist_opch16->SetName("opch16");
    hist_opch16->SetTitle("Response of OpCh16");
    hist_opch16->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch16->GetYaxis()->SetTitle("#photon / event");
    hist_opch16->SetMarkerStyle(21);
    hist_opch16->SetMarkerSize(1.0);






    //For opch22=================================================
    TFile *file22_1 = new TFile("../results/tmp/opCh22_"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    TFile *file22_2 = new TFile("../results/tmp/opCh22_"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file22_3 = new TFile("../results/tmp/opCh22_"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Get access to histos inside root files---
    TH2F *hist22_1 = dynamic_cast<TH2F*>(file22_1->Get("CRT_Opch")); 
    TH2F *hist22_2 = dynamic_cast<TH2F*>(file22_2->Get("CRT_Opch")); 
    TH2F *hist22_3 = dynamic_cast<TH2F*>(file22_3->Get("CRT_Opch")); 

    //Create new hist to store opch01---
    TH2F *hist_opch22 = new TH2F(*hist22_1);
    hist_opch22->Add(hist22_2);
    hist_opch22->Add(hist22_3);

    hist_opch22->SetName("opch22");
    hist_opch22->SetTitle("Response of OpCh22");
    hist_opch22->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch22->GetYaxis()->SetTitle("#photon / event");
    hist_opch22->SetMarkerStyle(21);
    hist_opch22->SetMarkerSize(1.0);




    
    //Output location:============================================
    string output_path = "/Users/shuaixiangzhang/Work/current/FNAL_Work2024/rsl_analyses/v4_analysis/results/fit_Develop_m2check/combine_3000results/";
    string output_name = rslxx + "_3000num_e67_crtCut.root";

    TFile *outputFile = new TFile(TString(output_path)+TString(output_name), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Unable to create combined_file.root" << std::endl;
        return 1;
    }


    // Write the combined histogram to the new ROOT file----------
    hist_opch01->Write();
    hist_opch03->Write();
    hist_opch16->Write();
    hist_opch22->Write();


    delete file01_1;
    delete file01_2;
    delete file01_3;
    delete file03_1;
    delete file03_2;
    delete file03_3;
    delete outputFile;

    return 0;
}
