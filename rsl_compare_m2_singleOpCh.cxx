#include <iostream>

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>


void readFile(const std::string& fileName, double distances[], double mpvConv[], double sigConv[], int numHist[], int maxArraySize) {
    //Open the input file for reading
    std::ifstream inputFile(fileName);

    //Check if the file is open
    if(!inputFile.is_open()) {
        std::cerr << "Error: Cannot open input file '" << fileName << "' for reading" << std::endl;
        return;
    }

    //Temporary variables to store data from each line
    double distance, mpv, sig;
    int num;
    int index = 0;
    
    //Read data from each line and store them into arrays
    while (inputFile >> distance >> mpv >> sig >> num && index < maxArraySize) {
        distances[index] = distance;
        mpvConv[index] = mpv;
        sigConv[index] = sig;
        numHist[index] = num;
        index++;
    }

    //Close the input file
    inputFile.close();
}



double sigBias(double A, double sigA, double B, double sigB){
    double term1 = 4/((A+B)*(A+B));
    double term2 = sqrt(B*B*sigA*sigA + A*A*sigB*sigB);

    return term1*term2;
}



//======TMP FUNCTION=======================================
void rsl_compareTMP(const std::string& path){

    const int num = 60; //num of slices---

    std::string inputFile99_opch01 = path + "/output_rsl99_opch01.txt";
    std::string inputFile99_opch03 = path + "/output_rsl99_opch03.txt";
    std::string inputFile99_opch16 = path + "/output_rsl99_opch16.txt";
    std::string inputFile99_opch22 = path + "/output_rsl99_opch22.txt";


    double dis99_opch01[num];
    double mpv99_opch01[num];
    double sig99_opch01[num];
    int ent99_opch01[num];//entries---
    readFile(inputFile99_opch01, dis99_opch01, mpv99_opch01, sig99_opch01, ent99_opch01, num);
    double sigMPV99_opch01[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent99_opch01[i]>0){
            sigMPV99_opch01[i] = sig99_opch01[i]/sqrt(static_cast<double>(ent99_opch01[i]));
        }
        else sigMPV99_opch01[i] = 0.0;
    }


    double dis99_opch03[num];
    double mpv99_opch03[num];
    double sig99_opch03[num];
    int ent99_opch03[num];//entries---
    readFile(inputFile99_opch03, dis99_opch03, mpv99_opch03, sig99_opch03, ent99_opch03, num);
    double sigMPV99_opch03[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent99_opch03[i]>0){
            sigMPV99_opch03[i] = sig99_opch03[i]/sqrt(static_cast<double>(ent99_opch03[i]));
        }
        else sigMPV99_opch03[i] = 0.0;
    }


    double dis99_opch16[num];
    double mpv99_opch16[num];
    double sig99_opch16[num];
    int ent99_opch16[num];//entries---
    readFile(inputFile99_opch16, dis99_opch16, mpv99_opch16, sig99_opch16, ent99_opch16, num);
    double sigMPV99_opch16[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent99_opch16[i]>0){
            sigMPV99_opch16[i] = sig99_opch16[i]/sqrt(static_cast<double>(ent99_opch16[i]));
        }
        else sigMPV99_opch16[i] = 0.0;
    }



    double dis99_opch22[num];
    double mpv99_opch22[num];
    double sig99_opch22[num];
    int ent99_opch22[num];//entries---
    readFile(inputFile99_opch22, dis99_opch22, mpv99_opch22, sig99_opch22, ent99_opch22, num);
    double sigMPV99_opch22[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent99_opch22[i]>0){
            sigMPV99_opch22[i] = sig99_opch22[i]/sqrt(static_cast<double>(ent99_opch22[i]));
        }
        else sigMPV99_opch22[i] = 0.0;
    }


    //Drawing 1-------------------------------------------------------------------
    TCanvas canvas("fitCompare", "Different RSLs", 800, 600);
    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);

    //RSL99 opch01---
    TGraphErrors* scat99_opch01 = new TGraphErrors(num, dis99_opch01, mpv99_opch01, nullptr, sigMPV99_opch01);
    scat99_opch01->SetLineColor(kRed);
    scat99_opch01->SetMarkerStyle(21);//20: filled square
    scat99_opch01->SetMarkerSize(0.6);
    scat99_opch01->SetMarkerColor(kRed);
    scat99_opch01->Draw("AP");
    scat99_opch01->GetXaxis()->SetTitle("Distance [cm]");
    scat99_opch01->GetYaxis()->SetTitle("#photon / event");
    scat99_opch01->GetXaxis()->SetRangeUser(0, 600);
    scat99_opch01->GetYaxis()->SetRangeUser(0, 2000);
//    TGraph* lineGraph = new TGraph(num, dis99, mpv99);
//    lineGraph->SetLineColor(kRed);
//    lineGraph->SetLineWidth(1);
//    lineGraph->Draw("L");

    //RSL99 opch03---
//    TGraph* scatterGraph50 = new TGraph(num, dis50, mpv50);
    TGraphErrors* scat99_opch03 = new TGraphErrors(num, dis99_opch03, mpv99_opch03, nullptr, sigMPV99_opch03);
    scat99_opch03->SetLineColor(kBlue);
    scat99_opch03->SetMarkerStyle(20);//20: filled circle
    scat99_opch03->SetMarkerSize(0.6);
    scat99_opch03->SetMarkerColor(kBlue);
    scat99_opch03->Draw("P SAME");


    //RSL99 opch16---
    TGraphErrors* scat99_opch16 = new TGraphErrors(num, dis99_opch16, mpv99_opch16, nullptr, sigMPV99_opch16);
    scat99_opch16->SetLineColor(kGreen);
    scat99_opch16->SetMarkerStyle(33);//33: filled rhombus
    scat99_opch16->SetMarkerSize(0.6);
    scat99_opch16->SetMarkerColor(kGreen);
    scat99_opch16->Draw("P SAME");


    //RSL99 opch22---
    TGraphErrors* scat99_opch22 = new TGraphErrors(num, dis99_opch22, mpv99_opch22, nullptr, sigMPV99_opch22);
    scat99_opch22->SetLineColor(kViolet);
    scat99_opch22->SetMarkerStyle(44);//44: four pointed star
    scat99_opch22->SetMarkerSize(0.6);
    scat99_opch22->SetMarkerColor(kViolet);
    scat99_opch22->Draw("P SAME");

    
    legend->AddEntry(scat99_opch01, "RSL = 100cm, opch01", "pe");
    legend->AddEntry(scat99_opch03, "RSL = 100cm, opch03", "pe");
    legend->AddEntry(scat99_opch16, "RSL = 100cm, opch16", "pe");
    legend->AddEntry(scat99_opch22, "RSL = 100cm, opch22", "pe");

    legend->Draw();




//Create root file--------------------------------------------------------
    TFile* outputFile = new TFile(TString(path)+"/compare.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file for writing" << std::endl;
        return;
    }










    canvas.Write();

    outputFile->Close();


    //TEST------------------
/*    for(int i=0; i<60; ++i){
        std::cout<<"mpv: "<<mpv99[i]<<std::endl;
    }
*/
}

    
//==========Main Function===============================================
void rsl_compare_m2_singleOpCh(){

    rsl_compareTMP("/Users/shuaixiangzhang/Work/current/FNAL_Work2024/rsl_analyses/v4_analysis/results/fit_Develop_m2check/membrane2");
        

}

    

