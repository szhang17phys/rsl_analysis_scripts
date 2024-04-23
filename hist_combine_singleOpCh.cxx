#include <iostream>
#include <TFile.h>
#include <TH2.h>


int hist_combine_singleOpCh(){

    string rslxx ="rsl150";

    //For opch00=================================================
    TFile *file00_1 = new TFile("../results/tmp/opCh00_"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    if (!file00_1 || file00_1->IsZombie()) {
        std::cerr << "Error: Unable to open file1.root" << std::endl;
        return 1;
    }
    //Open the second root file---
    TFile *file00_2 = new TFile("../results/tmp/opCh00_"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file00_3 = new TFile("../results/tmp/opCh00_"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Get access to histos inside root files---
    TH2F *hist00_1 = dynamic_cast<TH2F*>(file00_1->Get("CRT_Opch")); 
    TH2F *hist00_2 = dynamic_cast<TH2F*>(file00_2->Get("CRT_Opch")); 
    TH2F *hist00_3 = dynamic_cast<TH2F*>(file00_3->Get("CRT_Opch")); 

    //Create new hist to store opch00-----------------------
    TH2F *hist_opch00 = new TH2F(*hist00_1);
    hist_opch00->Add(hist00_2);
    hist_opch00->Add(hist00_3);

    hist_opch00->SetName("opch00");
    hist_opch00->SetTitle("Response of OpCh00");
    hist_opch00->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch00->GetYaxis()->SetTitle("#photon / event");
    hist_opch00->SetMarkerStyle(21);
    hist_opch00->SetMarkerSize(1.0);




    //For opch02=================================================
    TFile *file02_1 = new TFile("../results/tmp/opCh02_"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    TFile *file02_2 = new TFile("../results/tmp/opCh02_"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file02_3 = new TFile("../results/tmp/opCh02_"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Get access to histos inside root files---
    TH2F *hist02_1 = dynamic_cast<TH2F*>(file02_1->Get("CRT_Opch")); 
    TH2F *hist02_2 = dynamic_cast<TH2F*>(file02_2->Get("CRT_Opch")); 
    TH2F *hist02_3 = dynamic_cast<TH2F*>(file02_3->Get("CRT_Opch")); 

    //Create new hist to store opch01---
    TH2F *hist_opch02 = new TH2F(*hist02_1);
    hist_opch02->Add(hist02_2);
    hist_opch02->Add(hist02_3);

    hist_opch02->SetName("opch02");
    hist_opch02->SetTitle("Response of OpCh02");
    hist_opch02->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch02->GetYaxis()->SetTitle("#photon / event");
    hist_opch02->SetMarkerStyle(21);
    hist_opch02->SetMarkerSize(1.0);



    //For opch17=================================================
    TFile *file17_1 = new TFile("../results/tmp/opCh17_"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    TFile *file17_2 = new TFile("../results/tmp/opCh17_"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file17_3 = new TFile("../results/tmp/opCh17_"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Get access to histos inside root files---
    TH2F *hist17_1 = dynamic_cast<TH2F*>(file17_1->Get("CRT_Opch")); 
    TH2F *hist17_2 = dynamic_cast<TH2F*>(file17_2->Get("CRT_Opch")); 
    TH2F *hist17_3 = dynamic_cast<TH2F*>(file17_3->Get("CRT_Opch")); 

    //Create new hist to store opch01---
    TH2F *hist_opch17 = new TH2F(*hist17_1);
    hist_opch17->Add(hist17_2);
    hist_opch17->Add(hist17_3);

    hist_opch17->SetName("opch17");
    hist_opch17->SetTitle("Response of OpCh17");
    hist_opch17->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch17->GetYaxis()->SetTitle("#photon / event");
    hist_opch17->SetMarkerStyle(21);
    hist_opch17->SetMarkerSize(1.0);






    //For opch23=================================================
    TFile *file23_1 = new TFile("../results/tmp/opCh23_"+TString(rslxx)+"_1000num_e67_hist.root", "READ");
    TFile *file23_2 = new TFile("../results/tmp/opCh23_"+TString(rslxx)+"_1000num2_e67_hist.root", "READ");
    TFile *file23_3 = new TFile("../results/tmp/opCh23_"+TString(rslxx)+"_1000num3_e67_hist.root", "READ");

    //Get access to histos inside root files---
    TH2F *hist23_1 = dynamic_cast<TH2F*>(file23_1->Get("CRT_Opch")); 
    TH2F *hist23_2 = dynamic_cast<TH2F*>(file23_2->Get("CRT_Opch")); 
    TH2F *hist23_3 = dynamic_cast<TH2F*>(file23_3->Get("CRT_Opch")); 

    //Create new hist to store opch01---
    TH2F *hist_opch23 = new TH2F(*hist23_1);
    hist_opch23->Add(hist23_2);
    hist_opch23->Add(hist23_3);

    hist_opch23->SetName("opch23");
    hist_opch23->SetTitle("Response of OpCh23");
    hist_opch23->GetXaxis()->SetTitle("Distance [cm]");
	hist_opch23->GetYaxis()->SetTitle("#photon / event");
    hist_opch23->SetMarkerStyle(21);
    hist_opch23->SetMarkerSize(1.0);




    
    //Output location:============================================
    string output_path = "/Users/shuaixiangzhang/Work/current/FNAL_Work2024/rsl_analyses/v4_analysis/results/fit_Develop/combine_3000results/";
    string output_name = "singleOpCh_" + rslxx + "_3000num_e67_crtCut.root";

    TFile *outputFile = new TFile(TString(output_path)+TString(output_name), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Unable to create combined_file.root" << std::endl;
        return 1;
    }


    // Write the combined histogram to the new ROOT file----------
    hist_opch00->Write();
    hist_opch02->Write();
    hist_opch17->Write();
    hist_opch23->Write();


    delete outputFile;

    return 0;
}