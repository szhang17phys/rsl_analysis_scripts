#include <iostream>
#include <TFile.h>
#include <TH2.h>


int singleOpCh_combine(){
    //Input & output root files===========================

    string rslxx ="opCh01_rsl99";

    //Open the first ROOT file---
    TFile *file0 = new TFile("../results/tmp/"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    if (!file0 || file0->IsZombie()) {
        std::cerr << "Error: Unable to open file1.root" << std::endl;
        return 1;
    }
    //Open the second root file---
    TFile *file1 = new TFile("../results/tmp/"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file2 = new TFile("../results/tmp/"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Output location:------------
//    string output_path = "../results/combine_2000results/";

    string output_path = "/Users/shuaixiangzhang/Work/current/FNAL_Work2024/rsl_analyses/v4_analysis/results/fit_Develop_m2check/combine_3000results/";

    string output_name = rslxx + "_3000num_e67_crtCut.root";
    //===================================================


    //Get access to histos inside root files----------------------------
    TH2F *hist10 = dynamic_cast<TH2F*>(file0->Get("CRT_Opch")); 
    TH2F *hist11 = dynamic_cast<TH2F*>(file1->Get("CRT_Opch")); 
    TH2F *hist12 = dynamic_cast<TH2F*>(file2->Get("CRT_Opch")); 

    //Create new hist to store opch01-----------------------
    TH2F *hist_opch01 = new TH2F(*hist10);
    hist_opch01->Add(hist11);
    hist_opch01->Add(hist12);

    hist_opch01->SetName("opch01");
    hist_opch01->SetTitle("Response of OpCh01");
    hist_opch01->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch01->GetYaxis()->SetTitle("#photon / event");
    hist_opch01->SetMarkerStyle(21);
    hist_opch01->SetMarkerSize(1.0);


    


   //Create NEW root file--------------------------------------------------
    TFile *outputFile = new TFile(TString(output_path)+TString(output_name), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Unable to create combined_file.root" << std::endl;
        return 1;
    }

    // Write the combined histogram to the new ROOT file----------
    hist_opch01->Write();


    delete file0;
    delete file1;
    delete file2;   
    delete outputFile;

    return 0;
}
