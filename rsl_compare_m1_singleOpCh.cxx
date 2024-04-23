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

    std::string inputFile99_opch01 = path + "/output_rsl99_opch00.txt";
    std::string inputFile99_opch03 = path + "/output_rsl99_opch02.txt";
    std::string inputFile99_opch16 = path + "/output_rsl99_opch17.txt";
    std::string inputFile99_opch22 = path + "/output_rsl99_opch23.txt";

    std::string inputFile50_opch01 = path + "/output_rsl50_opch00.txt";
    std::string inputFile50_opch03 = path + "/output_rsl50_opch02.txt";
    std::string inputFile50_opch16 = path + "/output_rsl50_opch17.txt";
    std::string inputFile50_opch22 = path + "/output_rsl50_opch23.txt";

    std::string inputFile70_opch01 = path + "/output_rsl70_opch00.txt";
    std::string inputFile70_opch03 = path + "/output_rsl70_opch02.txt";
    std::string inputFile70_opch16 = path + "/output_rsl70_opch17.txt";
    std::string inputFile70_opch22 = path + "/output_rsl70_opch23.txt";

    std::string inputFile130_opch01 = path + "/output_rsl130_opch00.txt";
    std::string inputFile130_opch03 = path + "/output_rsl130_opch02.txt";
    std::string inputFile130_opch16 = path + "/output_rsl130_opch17.txt";
    std::string inputFile130_opch22 = path + "/output_rsl130_opch23.txt";

    std::string inputFile150_opch01 = path + "/output_rsl150_opch00.txt";
    std::string inputFile150_opch03 = path + "/output_rsl150_opch02.txt";
    std::string inputFile150_opch16 = path + "/output_rsl150_opch17.txt";
    std::string inputFile150_opch22 = path + "/output_rsl150_opch23.txt";



    //rsl99-------------------------------
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



    //rsl50-------------------------------
    double dis50_opch01[num];
    double mpv50_opch01[num];
    double sig50_opch01[num];
    int ent50_opch01[num];//entries---
    readFile(inputFile50_opch01, dis50_opch01, mpv50_opch01, sig50_opch01, ent50_opch01, num);
    double sigMPV50_opch01[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent50_opch01[i]>0){
            sigMPV50_opch01[i] = sig50_opch01[i]/sqrt(static_cast<double>(ent50_opch01[i]));
        }
        else sigMPV50_opch01[i] = 0.0;
    }

    double dis50_opch03[num];
    double mpv50_opch03[num];
    double sig50_opch03[num];
    int ent50_opch03[num];//entries---
    readFile(inputFile50_opch03, dis50_opch03, mpv50_opch03, sig50_opch03, ent50_opch03, num);
    double sigMPV50_opch03[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent50_opch03[i]>0){
            sigMPV50_opch03[i] = sig50_opch03[i]/sqrt(static_cast<double>(ent50_opch03[i]));
        }
        else sigMPV50_opch03[i] = 0.0;
    }

    double dis50_opch16[num];
    double mpv50_opch16[num];
    double sig50_opch16[num];
    int ent50_opch16[num];//entries---
    readFile(inputFile50_opch16, dis50_opch16, mpv50_opch16, sig50_opch16, ent50_opch16, num);
    double sigMPV50_opch16[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent50_opch16[i]>0){
            sigMPV50_opch16[i] = sig50_opch16[i]/sqrt(static_cast<double>(ent50_opch16[i]));
        }
        else sigMPV50_opch16[i] = 0.0;
    }

    double dis50_opch22[num];
    double mpv50_opch22[num];
    double sig50_opch22[num];
    int ent50_opch22[num];//entries---
    readFile(inputFile50_opch22, dis50_opch22, mpv50_opch22, sig50_opch22, ent50_opch22, num);
    double sigMPV50_opch22[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent50_opch22[i]>0){
            sigMPV50_opch22[i] = sig50_opch22[i]/sqrt(static_cast<double>(ent50_opch22[i]));
        }
        else sigMPV50_opch22[i] = 0.0;
    }


    //rsl70-------------------------------
    double dis70_opch01[num];
    double mpv70_opch01[num];
    double sig70_opch01[num];
    int ent70_opch01[num];//entries---
    readFile(inputFile70_opch01, dis70_opch01, mpv70_opch01, sig70_opch01, ent70_opch01, num);
    double sigMPV70_opch01[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent70_opch01[i]>0){
            sigMPV70_opch01[i] = sig70_opch01[i]/sqrt(static_cast<double>(ent70_opch01[i]));
        }
        else sigMPV70_opch01[i] = 0.0;
    }

    double dis70_opch03[num];
    double mpv70_opch03[num];
    double sig70_opch03[num];
    int ent70_opch03[num];//entries---
    readFile(inputFile70_opch03, dis70_opch03, mpv70_opch03, sig70_opch03, ent70_opch03, num);
    double sigMPV70_opch03[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent70_opch03[i]>0){
            sigMPV70_opch03[i] = sig70_opch03[i]/sqrt(static_cast<double>(ent70_opch03[i]));
        }
        else sigMPV70_opch03[i] = 0.0;
    }

    double dis70_opch16[num];
    double mpv70_opch16[num];
    double sig70_opch16[num];
    int ent70_opch16[num];//entries---
    readFile(inputFile70_opch16, dis70_opch16, mpv70_opch16, sig70_opch16, ent70_opch16, num);
    double sigMPV70_opch16[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent70_opch16[i]>0){
            sigMPV70_opch16[i] = sig70_opch16[i]/sqrt(static_cast<double>(ent70_opch16[i]));
        }
        else sigMPV70_opch16[i] = 0.0;
    }

    double dis70_opch22[num];
    double mpv70_opch22[num];
    double sig70_opch22[num];
    int ent70_opch22[num];//entries---
    readFile(inputFile70_opch22, dis70_opch22, mpv70_opch22, sig70_opch22, ent70_opch22, num);
    double sigMPV70_opch22[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent70_opch22[i]>0){
            sigMPV70_opch22[i] = sig70_opch22[i]/sqrt(static_cast<double>(ent70_opch22[i]));
        }
        else sigMPV70_opch22[i] = 0.0;
    }


    //rsl130-------------------------------
    double dis130_opch01[num];
    double mpv130_opch01[num];
    double sig130_opch01[num];
    int ent130_opch01[num];//entries---
    readFile(inputFile130_opch01, dis130_opch01, mpv130_opch01, sig130_opch01, ent130_opch01, num);
    double sigMPV130_opch01[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent130_opch01[i]>0){
            sigMPV130_opch01[i] = sig130_opch01[i]/sqrt(static_cast<double>(ent130_opch01[i]));
        }
        else sigMPV130_opch01[i] = 0.0;
    }

    double dis130_opch03[num];
    double mpv130_opch03[num];
    double sig130_opch03[num];
    int ent130_opch03[num];//entries---
    readFile(inputFile130_opch03, dis130_opch03, mpv130_opch03, sig130_opch03, ent130_opch03, num);
    double sigMPV130_opch03[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent130_opch03[i]>0){
            sigMPV130_opch03[i] = sig130_opch03[i]/sqrt(static_cast<double>(ent130_opch03[i]));
        }
        else sigMPV130_opch03[i] = 0.0;
    }

    double dis130_opch16[num];
    double mpv130_opch16[num];
    double sig130_opch16[num];
    int ent130_opch16[num];//entries---
    readFile(inputFile130_opch16, dis130_opch16, mpv130_opch16, sig130_opch16, ent130_opch16, num);
    double sigMPV130_opch16[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent130_opch16[i]>0){
            sigMPV130_opch16[i] = sig130_opch16[i]/sqrt(static_cast<double>(ent130_opch16[i]));
        }
        else sigMPV130_opch16[i] = 0.0;
    }

    double dis130_opch22[num];
    double mpv130_opch22[num];
    double sig130_opch22[num];
    int ent130_opch22[num];//entries---
    readFile(inputFile130_opch22, dis130_opch22, mpv130_opch22, sig130_opch22, ent130_opch22, num);
    double sigMPV130_opch22[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent130_opch22[i]>0){
            sigMPV130_opch22[i] = sig130_opch22[i]/sqrt(static_cast<double>(ent130_opch22[i]));
        }
        else sigMPV130_opch22[i] = 0.0;
    }


    //rsl150-------------------------------
    double dis150_opch01[num];
    double mpv150_opch01[num];
    double sig150_opch01[num];
    int ent150_opch01[num];//entries---
    readFile(inputFile150_opch01, dis150_opch01, mpv150_opch01, sig150_opch01, ent150_opch01, num);
    double sigMPV150_opch01[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent150_opch01[i]>0){
            sigMPV150_opch01[i] = sig150_opch01[i]/sqrt(static_cast<double>(ent150_opch01[i]));
        }
        else sigMPV150_opch01[i] = 0.0;
    }

    double dis150_opch03[num];
    double mpv150_opch03[num];
    double sig150_opch03[num];
    int ent150_opch03[num];//entries---
    readFile(inputFile150_opch03, dis150_opch03, mpv150_opch03, sig150_opch03, ent150_opch03, num);
    double sigMPV150_opch03[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent150_opch03[i]>0){
            sigMPV150_opch03[i] = sig150_opch03[i]/sqrt(static_cast<double>(ent150_opch03[i]));
        }
        else sigMPV150_opch03[i] = 0.0;
    }

    double dis150_opch16[num];
    double mpv150_opch16[num];
    double sig150_opch16[num];
    int ent150_opch16[num];//entries---
    readFile(inputFile150_opch16, dis150_opch16, mpv150_opch16, sig150_opch16, ent150_opch16, num);
    double sigMPV150_opch16[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent150_opch16[i]>0){
            sigMPV150_opch16[i] = sig150_opch16[i]/sqrt(static_cast<double>(ent150_opch16[i]));
        }
        else sigMPV150_opch16[i] = 0.0;
    }

    double dis150_opch22[num];
    double mpv150_opch22[num];
    double sig150_opch22[num];
    int ent150_opch22[num];//entries---
    readFile(inputFile150_opch22, dis150_opch22, mpv150_opch22, sig150_opch22, ent150_opch22, num);
    double sigMPV150_opch22[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent150_opch22[i]>0){
            sigMPV150_opch22[i] = sig150_opch22[i]/sqrt(static_cast<double>(ent150_opch22[i]));
        }
        else sigMPV150_opch22[i] = 0.0;
    }




    //Drawing 1======================================================================
    TCanvas canvas("fitCompare", "Different RSLs", 800, 600);
    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);

    //RSL99 opch01---------------------------------------------
    TGraphErrors* scat99_opch01 = new TGraphErrors(num, dis99_opch01, mpv99_opch01, nullptr, sigMPV99_opch01);
    scat99_opch01->SetLineColor(kRed);
    scat99_opch01->SetMarkerStyle(20);//21: filled square
    scat99_opch01->SetMarkerSize(0.6);
    scat99_opch01->SetMarkerColor(kRed);
//    scat99_opch01->Draw("AP");
    scat99_opch01->GetXaxis()->SetTitle("Distance [cm]");
    scat99_opch01->GetYaxis()->SetTitle("#photon / event");
    scat99_opch01->GetXaxis()->SetRangeUser(0, 300);
//    scat99_opch01->GetYaxis()->SetRangeUser(0, 2000);
//    TGraph* lineGraph = new TGraph(num, dis99, mpv99);
//    lineGraph->SetLineColor(kRed);
//    lineGraph->SetLineWidth(1);
//    lineGraph->Draw("L");

    //RSL99 opch03---
//    TGraph* scatterGraph50 = new TGraph(num, dis50, mpv50);
    TGraphErrors* scat99_opch03 = new TGraphErrors(num, dis99_opch03, mpv99_opch03, nullptr, sigMPV99_opch03);
    scat99_opch03->SetLineColor(kRed);
    scat99_opch03->SetMarkerStyle(21);//20: filled circle
    scat99_opch03->SetMarkerSize(0.6);
    scat99_opch03->SetMarkerColor(kRed);
//    scat99_opch03->Draw("P SAME");


    //RSL99 opch16---
    TGraphErrors* scat99_opch16 = new TGraphErrors(num, dis99_opch16, mpv99_opch16, nullptr, sigMPV99_opch16);
    scat99_opch16->SetLineColor(kRed);
    scat99_opch16->SetMarkerStyle(44);//22: normal triangle
    scat99_opch16->SetMarkerSize(0.7);
    scat99_opch16->SetMarkerColor(kRed);
//    scat99_opch16->Draw("P SAME");


    //RSL99 opch22---
    TGraphErrors* scat99_opch22 = new TGraphErrors(num, dis99_opch22, mpv99_opch22, nullptr, sigMPV99_opch22);
    scat99_opch22->SetLineColor(kRed);
    scat99_opch22->SetMarkerStyle(22);//44: four pointed star
    scat99_opch22->SetMarkerSize(0.6);
    scat99_opch22->SetMarkerColor(kRed);
//    scat99_opch22->Draw("P SAME");
    scat99_opch22->Draw("AP");

    legend->AddEntry(scat99_opch01, "RSL = 100cm", "pe");
//    legend->AddEntry(scat99_opch03, "RSL = 100cm, opch03", "pe");
//    legend->AddEntry(scat99_opch16, "RSL = 100cm, opch16", "pe");
//    legend->AddEntry(scat99_opch22, "RSL = 100cm, opch22", "pe");


    //RSL50 opch01---------------------------------------------
    TGraphErrors* scat50_opch01 = new TGraphErrors(num, dis50_opch01, mpv50_opch01, nullptr, sigMPV50_opch01);
    scat50_opch01->SetLineColor(kGreen);
    scat50_opch01->SetMarkerStyle(20);//21: filled square
    scat50_opch01->SetMarkerSize(0.6);
    scat50_opch01->SetMarkerColor(kGreen);
//    scat50_opch01->Draw("P SAME");

    //RSL50 opch03---
    TGraphErrors* scat50_opch03 = new TGraphErrors(num, dis50_opch03, mpv50_opch03, nullptr, sigMPV50_opch03);
    scat50_opch03->SetLineColor(kGreen);
    scat50_opch03->SetMarkerStyle(21);//20: filled circle
    scat50_opch03->SetMarkerSize(0.6);
    scat50_opch03->SetMarkerColor(kGreen);
//    scat50_opch03->Draw("P SAME");

    //RSL50 opch16---
    TGraphErrors* scat50_opch16 = new TGraphErrors(num, dis50_opch16, mpv50_opch16, nullptr, sigMPV50_opch16);
    scat50_opch16->SetLineColor(kGreen);
    scat50_opch16->SetMarkerStyle(44);//22: normal triangle
    scat50_opch16->SetMarkerSize(0.7);
    scat50_opch16->SetMarkerColor(kGreen);
//    scat50_opch16->Draw("P SAME");

    //RSL50 opch22---
    TGraphErrors* scat50_opch22 = new TGraphErrors(num, dis50_opch22, mpv50_opch22, nullptr, sigMPV50_opch22);
    scat50_opch22->SetLineColor(kGreen);
    scat50_opch22->SetMarkerStyle(22);//44: four pointed star
    scat50_opch22->SetMarkerSize(0.6);
    scat50_opch22->SetMarkerColor(kGreen);
    scat50_opch22->Draw("P SAME");

    legend->AddEntry(scat50_opch01, "RSL = 50cm", "pe");
//    legend->AddEntry(scat50_opch03, "RSL = 100cm, opch03", "pe");
//    legend->AddEntry(scat50_opch16, "RSL = 100cm, opch16", "pe");
//    legend->AddEntry(scat50_opch22, "RSL = 100cm, opch22", "pe");


    //RSL70 opch01---------------------------------------------
    TGraphErrors* scat70_opch01 = new TGraphErrors(num, dis70_opch01, mpv70_opch01, nullptr, sigMPV70_opch01);
    scat70_opch01->SetLineColor(kBlue);
    scat70_opch01->SetMarkerStyle(20);//21: filled square
    scat70_opch01->SetMarkerSize(0.6);
    scat70_opch01->SetMarkerColor(kBlue);
//    scat70_opch01->Draw("P SAME");

    //RS70 opch03---
    TGraphErrors* scat70_opch03 = new TGraphErrors(num, dis70_opch03, mpv70_opch03, nullptr, sigMPV70_opch03);
    scat70_opch03->SetLineColor(kBlue);
    scat70_opch03->SetMarkerStyle(21);//20: filled circle
    scat70_opch03->SetMarkerSize(0.6);
    scat70_opch03->SetMarkerColor(kBlue);
//    scat70_opch03->Draw("P SAME");

    //RSL70 opch16---
    TGraphErrors* scat70_opch16 = new TGraphErrors(num, dis70_opch16, mpv70_opch16, nullptr, sigMPV70_opch16);
    scat70_opch16->SetLineColor(kBlue);
    scat70_opch16->SetMarkerStyle(44);//22: normal triangle
    scat70_opch16->SetMarkerSize(0.6);
    scat70_opch16->SetMarkerColor(kBlue);
//    scat70_opch16->Draw("P SAME");

    //RSL70 opch22---
    TGraphErrors* scat70_opch22 = new TGraphErrors(num, dis70_opch22, mpv70_opch22, nullptr, sigMPV70_opch22);
    scat70_opch22->SetLineColor(kBlue);
    scat70_opch22->SetMarkerStyle(22);//44: four pointed star
    scat70_opch22->SetMarkerSize(0.6);
    scat70_opch22->SetMarkerColor(kBlue);
    scat70_opch22->Draw("P SAME");

    legend->AddEntry(scat70_opch01, "RSL = 70cm", "pe");


    //RSL130 opch01---------------------------------------------
    TGraphErrors* scat130_opch01 = new TGraphErrors(num, dis130_opch01, mpv130_opch01, nullptr, sigMPV130_opch01);
    scat130_opch01->SetLineColor(28);
    scat130_opch01->SetMarkerStyle(20);//21: filled square
    scat130_opch01->SetMarkerSize(0.6);
    scat130_opch01->SetMarkerColor(28);
//    scat130_opch01->Draw("P SAME");

    //RS130 opch03---
    TGraphErrors* scat130_opch03 = new TGraphErrors(num, dis130_opch03, mpv130_opch03, nullptr, sigMPV130_opch03);
    scat130_opch03->SetLineColor(28);
    scat130_opch03->SetMarkerStyle(21);//20: filled circle
    scat130_opch03->SetMarkerSize(0.6);
    scat130_opch03->SetMarkerColor(28);
//    scat130_opch03->Draw("P SAME");

    //RSL130 opch16---
    TGraphErrors* scat130_opch16 = new TGraphErrors(num, dis130_opch16, mpv130_opch16, nullptr, sigMPV130_opch16);
    scat130_opch16->SetLineColor(28);
    scat130_opch16->SetMarkerStyle(44);//22: normal triangle
    scat130_opch16->SetMarkerSize(0.6);
    scat130_opch16->SetMarkerColor(28);
//    scat130_opch16->Draw("P SAME");

    //RSL130 opch22---
    TGraphErrors* scat130_opch22 = new TGraphErrors(num, dis130_opch22, mpv130_opch22, nullptr, sigMPV130_opch22);
    scat130_opch22->SetLineColor(28);
    scat130_opch22->SetMarkerStyle(22);//44: four pointed star
    scat130_opch22->SetMarkerSize(0.6);
    scat130_opch22->SetMarkerColor(28);
    scat130_opch22->Draw("P SAME");

    legend->AddEntry(scat130_opch01, "RSL = 130cm", "pe");


    //RSL150 opch01---------------------------------------------
    TGraphErrors* scat150_opch01 = new TGraphErrors(num, dis150_opch01, mpv150_opch01, nullptr, sigMPV150_opch01);
    scat150_opch01->SetLineColor(kViolet);
    scat150_opch01->SetMarkerStyle(20);//21: filled square
    scat150_opch01->SetMarkerSize(0.6);
    scat150_opch01->SetMarkerColor(kViolet);
//    scat150_opch01->Draw("P SAME");

    //RS150 opch03---
    TGraphErrors* scat150_opch03 = new TGraphErrors(num, dis150_opch03, mpv150_opch03, nullptr, sigMPV150_opch03);
    scat150_opch03->SetLineColor(kViolet);
    scat150_opch03->SetMarkerStyle(21);//20: filled circle
    scat150_opch03->SetMarkerSize(0.6);
    scat150_opch03->SetMarkerColor(kViolet);
//    scat150_opch03->Draw("P SAME");

    //RSL150 opch16---
    TGraphErrors* scat150_opch16 = new TGraphErrors(num, dis150_opch16, mpv150_opch16, nullptr, sigMPV150_opch16);
    scat150_opch16->SetLineColor(kViolet);
    scat150_opch16->SetMarkerStyle(44);//22: normal triangle
    scat150_opch16->SetMarkerSize(0.7);
    scat150_opch16->SetMarkerColor(kViolet);
//    scat150_opch16->Draw("P SAME");

    //RSL150 opch22---
    TGraphErrors* scat150_opch22 = new TGraphErrors(num, dis150_opch22, mpv150_opch22, nullptr, sigMPV150_opch22);
    scat150_opch22->SetLineColor(kViolet);
    scat150_opch22->SetMarkerStyle(22);//44: four pointed star
    scat150_opch22->SetMarkerSize(0.6);
    scat150_opch22->SetMarkerColor(kViolet);
    scat150_opch22->Draw("P SAME");

    legend->AddEntry(scat150_opch01, "RSL = 150cm", "pe");

    legend->Draw();




//Drawing 2==============================================================
    TCanvas canvas2("fitDiff", "Different RSLs", 800, 600);
    TLegend* legend2 = new TLegend(0.6, 0.6, 0.9, 0.9);

    double diff50_opch01[num];//(mpv50 - mpv99) (opch01)---
    double diff50_opch03[num];
    double diff50_opch16[num];
    double diff50_opch22[num];

    double diff70_opch01[num];
    double diff70_opch03[num];
    double diff70_opch16[num];
    double diff70_opch22[num];

    double diff130_opch01[num];
    double diff130_opch03[num];
    double diff130_opch16[num];
    double diff130_opch22[num];

    double diff150_opch01[num];
    double diff150_opch03[num];
    double diff150_opch16[num];
    double diff150_opch22[num];

    double sigDiff50_opch01[num];//sigma of diff50---
    double sigDiff50_opch03[num];
    double sigDiff50_opch16[num];
    double sigDiff50_opch22[num];

    double sigDiff70_opch01[num];
    double sigDiff70_opch03[num];
    double sigDiff70_opch16[num];
    double sigDiff70_opch22[num];

    double sigDiff130_opch01[num];
    double sigDiff130_opch03[num];
    double sigDiff130_opch16[num];
    double sigDiff130_opch22[num];

    double sigDiff150_opch01[num];
    double sigDiff150_opch03[num];
    double sigDiff150_opch16[num];
    double sigDiff150_opch22[num];


    for(int i=0; i<num; ++i){
        diff50_opch01[i] = mpv50_opch01[i] - mpv99_opch01[i];
        diff50_opch03[i] = mpv50_opch03[i] - mpv99_opch03[i];
        diff50_opch16[i] = mpv50_opch16[i] - mpv99_opch16[i];
        diff50_opch22[i] = mpv50_opch22[i] - mpv99_opch22[i];

        diff70_opch01[i] = mpv70_opch01[i] - mpv99_opch01[i];
        diff70_opch03[i] = mpv70_opch03[i] - mpv99_opch03[i];
        diff70_opch16[i] = mpv70_opch16[i] - mpv99_opch16[i];
        diff70_opch22[i] = mpv70_opch22[i] - mpv99_opch22[i];

        diff130_opch01[i] = mpv130_opch01[i] - mpv99_opch01[i];
        diff130_opch03[i] = mpv130_opch03[i] - mpv99_opch03[i];
        diff130_opch16[i] = mpv130_opch16[i] - mpv99_opch16[i];
        diff130_opch22[i] = mpv130_opch22[i] - mpv99_opch22[i];

        diff150_opch01[i] = mpv150_opch01[i] - mpv99_opch01[i];
        diff150_opch03[i] = mpv150_opch03[i] - mpv99_opch03[i];
        diff150_opch16[i] = mpv150_opch16[i] - mpv99_opch16[i];
        diff150_opch22[i] = mpv150_opch22[i] - mpv99_opch22[i];

        sigDiff50_opch01[i] = sqrt(sigMPV50_opch01[i]*sigMPV50_opch01[i] + sigMPV99_opch01[i]*sigMPV99_opch01[i]);
        sigDiff50_opch03[i] = sqrt(sigMPV50_opch03[i]*sigMPV50_opch03[i] + sigMPV99_opch03[i]*sigMPV99_opch03[i]);        
        sigDiff50_opch16[i] = sqrt(sigMPV50_opch16[i]*sigMPV50_opch16[i] + sigMPV99_opch16[i]*sigMPV99_opch16[i]);
        sigDiff50_opch22[i] = sqrt(sigMPV50_opch22[i]*sigMPV50_opch22[i] + sigMPV99_opch22[i]*sigMPV99_opch22[i]);

        sigDiff70_opch01[i] = sqrt(sigMPV70_opch01[i]*sigMPV70_opch01[i] + sigMPV99_opch01[i]*sigMPV99_opch01[i]);
        sigDiff70_opch03[i] = sqrt(sigMPV70_opch03[i]*sigMPV70_opch03[i] + sigMPV99_opch03[i]*sigMPV99_opch03[i]);        
        sigDiff70_opch16[i] = sqrt(sigMPV70_opch16[i]*sigMPV70_opch16[i] + sigMPV99_opch16[i]*sigMPV99_opch16[i]);
        sigDiff70_opch22[i] = sqrt(sigMPV70_opch22[i]*sigMPV70_opch22[i] + sigMPV99_opch22[i]*sigMPV99_opch22[i]);

        sigDiff130_opch01[i] = sqrt(sigMPV130_opch01[i]*sigMPV130_opch01[i] + sigMPV99_opch01[i]*sigMPV99_opch01[i]);
        sigDiff130_opch03[i] = sqrt(sigMPV130_opch03[i]*sigMPV130_opch03[i] + sigMPV99_opch03[i]*sigMPV99_opch03[i]);        
        sigDiff130_opch16[i] = sqrt(sigMPV130_opch16[i]*sigMPV130_opch16[i] + sigMPV99_opch16[i]*sigMPV99_opch16[i]);
        sigDiff130_opch22[i] = sqrt(sigMPV130_opch22[i]*sigMPV130_opch22[i] + sigMPV99_opch22[i]*sigMPV99_opch22[i]);

        sigDiff150_opch01[i] = sqrt(sigMPV150_opch01[i]*sigMPV150_opch01[i] + sigMPV99_opch01[i]*sigMPV99_opch01[i]);
        sigDiff150_opch03[i] = sqrt(sigMPV150_opch03[i]*sigMPV150_opch03[i] + sigMPV99_opch03[i]*sigMPV99_opch03[i]);        
        sigDiff150_opch16[i] = sqrt(sigMPV150_opch16[i]*sigMPV150_opch16[i] + sigMPV99_opch16[i]*sigMPV99_opch16[i]);
        sigDiff150_opch22[i] = sqrt(sigMPV150_opch22[i]*sigMPV150_opch22[i] + sigMPV99_opch22[i]*sigMPV99_opch22[i]);

    }

    //diff50_opch01------------------------------------
    TGraphErrors* graphDiff50_opch01 = new TGraphErrors(num, dis50_opch01, diff50_opch01, nullptr, sigDiff50_opch01);
    graphDiff50_opch01->SetLineColor(kGreen);
    graphDiff50_opch01->SetMarkerStyle(20);
    graphDiff50_opch01->SetMarkerSize(0.6);
    graphDiff50_opch01->SetMarkerColor(kGreen);
//    graphDiff50_opch01->Draw("AP");
    graphDiff50_opch01->GetXaxis()->SetTitle("Distance [cm]");
    graphDiff50_opch01->GetYaxis()->SetTitle("Photon Num Diff");
    graphDiff50_opch01->GetXaxis()->SetRangeUser(0, 300);
//    graphDiff50_opch01->GetYaxis()->SetRangeUser(-500, 1000);

    //diff50_opch03
    TGraphErrors* graphDiff50_opch03 = new TGraphErrors(num, dis50_opch03, diff50_opch03, nullptr, sigDiff50_opch03);
    graphDiff50_opch03->SetLineColor(kGreen);
    graphDiff50_opch03->SetMarkerStyle(21);//20: filled circle
    graphDiff50_opch03->SetMarkerSize(0.6);
    graphDiff50_opch03->SetMarkerColor(kGreen);
//    graphDiff50_opch03->Draw("P SAME");

    graphDiff50_opch03->Draw("AP");
    graphDiff50_opch03->GetXaxis()->SetTitle("Distance [cm]");
    graphDiff50_opch03->GetYaxis()->SetTitle("Photon Num Diff");
    graphDiff50_opch03->GetXaxis()->SetRangeUser(0, 300);
//    graphDiff50_opch03->GetYaxis()->SetRangeUser(-500, 1000);

    //diff50_opch16
    TGraphErrors* graphDiff50_opch16 = new TGraphErrors(num, dis50_opch16, diff50_opch16, nullptr, sigDiff50_opch16);
    graphDiff50_opch16->SetLineColor(kGreen);
    graphDiff50_opch16->SetMarkerStyle(44);
    graphDiff50_opch16->SetMarkerSize(0.6);
    graphDiff50_opch16->SetMarkerColor(kGreen);
//    graphDiff50_opch16->Draw("P SAME");

    //diff50_opch22
    TGraphErrors* graphDiff50_opch22 = new TGraphErrors(num, dis50_opch22, diff50_opch22, nullptr, sigDiff50_opch22);
    graphDiff50_opch22->SetLineColor(kGreen);
    graphDiff50_opch22->SetMarkerStyle(22);
    graphDiff50_opch22->SetMarkerSize(0.6);
    graphDiff50_opch22->SetMarkerColor(kGreen);
    graphDiff50_opch22->Draw("P SAME");

    legend2->AddEntry(graphDiff50_opch01, "RSL50 - RSL100", "pe");



    //diff70_opch01------------------------------------
    TGraphErrors* graphDiff70_opch01 = new TGraphErrors(num, dis70_opch01, diff70_opch01, nullptr, sigDiff70_opch01);
    graphDiff70_opch01->SetLineColor(kBlue);
    graphDiff70_opch01->SetMarkerStyle(20);
    graphDiff70_opch01->SetMarkerSize(0.6);
    graphDiff70_opch01->SetMarkerColor(kBlue);
//    graphDiff70_opch01->Draw("P SAME");

    //diff70_opch03
    TGraphErrors* graphDiff70_opch03 = new TGraphErrors(num, dis70_opch03, diff70_opch03, nullptr, sigDiff70_opch03);
    graphDiff70_opch03->SetLineColor(kBlue);
    graphDiff70_opch03->SetMarkerStyle(21);//20: filled circle
    graphDiff70_opch03->SetMarkerSize(0.6);
    graphDiff70_opch03->SetMarkerColor(kBlue);
    graphDiff70_opch03->Draw("P SAME");

    //diff70_opch16
    TGraphErrors* graphDiff70_opch16 = new TGraphErrors(num, dis70_opch16, diff70_opch16, nullptr, sigDiff70_opch16);
    graphDiff70_opch16->SetLineColor(kBlue);
    graphDiff70_opch16->SetMarkerStyle(44);
    graphDiff70_opch16->SetMarkerSize(0.6);
    graphDiff70_opch16->SetMarkerColor(kBlue);
//    graphDiff70_opch16->Draw("P SAME");

    //diff70_opch22
    TGraphErrors* graphDiff70_opch22 = new TGraphErrors(num, dis70_opch22, diff70_opch22, nullptr, sigDiff70_opch22);
    graphDiff70_opch22->SetLineColor(kBlue);
    graphDiff70_opch22->SetMarkerStyle(22);
    graphDiff70_opch22->SetMarkerSize(0.6);
    graphDiff70_opch22->SetMarkerColor(kBlue);
    graphDiff70_opch22->Draw("P SAME");

    legend2->AddEntry(graphDiff70_opch01, "RSL70 - RSL100", "pe");



    //diff130_opch01------------------------------------
    TGraphErrors* graphDiff130_opch01 = new TGraphErrors(num, dis130_opch01, diff130_opch01, nullptr, sigDiff130_opch01);
    graphDiff130_opch01->SetLineColor(28);
    graphDiff130_opch01->SetMarkerStyle(20);
    graphDiff130_opch01->SetMarkerSize(0.6);
    graphDiff130_opch01->SetMarkerColor(28);
//    graphDiff130_opch01->Draw("P SAME");

    //diff130_opch03
    TGraphErrors* graphDiff130_opch03 = new TGraphErrors(num, dis130_opch03, diff130_opch03, nullptr, sigDiff130_opch03);
    graphDiff130_opch03->SetLineColor(28);
    graphDiff130_opch03->SetMarkerStyle(21);//20: filled circle
    graphDiff130_opch03->SetMarkerSize(0.6);
    graphDiff130_opch03->SetMarkerColor(28);
    graphDiff130_opch03->Draw("P SAME");

    //diff130_opch16
    TGraphErrors* graphDiff130_opch16 = new TGraphErrors(num, dis130_opch16, diff130_opch16, nullptr, sigDiff130_opch16);
    graphDiff130_opch16->SetLineColor(28);
    graphDiff130_opch16->SetMarkerStyle(44);
    graphDiff130_opch16->SetMarkerSize(0.6);
    graphDiff130_opch16->SetMarkerColor(28);
//    graphDiff130_opch16->Draw("P SAME");

    //diff130_opch22
    TGraphErrors* graphDiff130_opch22 = new TGraphErrors(num, dis130_opch22, diff130_opch22, nullptr, sigDiff130_opch22);
    graphDiff130_opch22->SetLineColor(28);
    graphDiff130_opch22->SetMarkerStyle(22);
    graphDiff130_opch22->SetMarkerSize(0.6);
    graphDiff130_opch22->SetMarkerColor(28);
    graphDiff130_opch22->Draw("P SAME");

    legend2->AddEntry(graphDiff130_opch01, "RSL130 - RSL100", "pe");


    //diff150_opch01------------------------------------
    TGraphErrors* graphDiff150_opch01 = new TGraphErrors(num, dis150_opch01, diff150_opch01, nullptr, sigDiff150_opch01);
    graphDiff150_opch01->SetLineColor(kViolet);
    graphDiff150_opch01->SetMarkerStyle(20);
    graphDiff150_opch01->SetMarkerSize(0.6);
    graphDiff150_opch01->SetMarkerColor(kViolet);
//    graphDiff150_opch01->Draw("P SAME");

    //diff150_opch03
    TGraphErrors* graphDiff150_opch03 = new TGraphErrors(num, dis150_opch03, diff150_opch03, nullptr, sigDiff150_opch03);
    graphDiff150_opch03->SetLineColor(kViolet);
    graphDiff150_opch03->SetMarkerStyle(21);//20: filled circle
    graphDiff150_opch03->SetMarkerSize(0.6);
    graphDiff150_opch03->SetMarkerColor(kViolet);
    graphDiff150_opch03->Draw("P SAME");

    //diff150_opch16
    TGraphErrors* graphDiff150_opch16 = new TGraphErrors(num, dis150_opch16, diff150_opch16, nullptr, sigDiff150_opch16);
    graphDiff150_opch16->SetLineColor(kViolet);
    graphDiff150_opch16->SetMarkerStyle(44);
    graphDiff150_opch16->SetMarkerSize(0.6);
    graphDiff150_opch16->SetMarkerColor(kViolet);
//    graphDiff150_opch16->Draw("P SAME");

    //diff150_opch22
    TGraphErrors* graphDiff150_opch22 = new TGraphErrors(num, dis150_opch22, diff150_opch22, nullptr, sigDiff150_opch22);
    graphDiff150_opch22->SetLineColor(kViolet);
    graphDiff150_opch22->SetMarkerStyle(22);
    graphDiff150_opch22->SetMarkerSize(0.6);
    graphDiff150_opch22->SetMarkerColor(kViolet);
    graphDiff150_opch22->Draw("P SAME");

    legend2->AddEntry(graphDiff150_opch01, "RSL150 - RSL100", "pe");


    legend2->Draw();





//Drawing 3================================================================
    TCanvas canvas3("fitBias", "fit Bias", 800, 600);
    TLegend* legend3 = new TLegend(0.6, 0.6, 0.9, 0.9);

    double bias50_opch01[num];//2*(mpv50-mpv99)/(mpv50+mpv99)---
    double bias50_opch03[num];
    double bias50_opch16[num];
    double bias50_opch22[num];

    double bias70_opch01[num];
    double bias70_opch03[num];
    double bias70_opch16[num];
    double bias70_opch22[num];

    double bias130_opch01[num];
    double bias130_opch03[num];
    double bias130_opch16[num];
    double bias130_opch22[num];

    double bias150_opch01[num];
    double bias150_opch03[num];
    double bias150_opch16[num];
    double bias150_opch22[num];

    double sigBias50_opch01[num];//sigma of bias50---
    double sigBias50_opch03[num];
    double sigBias50_opch16[num];
    double sigBias50_opch22[num];

    double sigBias70_opch01[num];
    double sigBias70_opch03[num];
    double sigBias70_opch16[num];
    double sigBias70_opch22[num];
 
    double sigBias130_opch01[num];
    double sigBias130_opch03[num];
    double sigBias130_opch16[num];
    double sigBias130_opch22[num];

    double sigBias150_opch01[num];
    double sigBias150_opch03[num];
    double sigBias150_opch16[num];
    double sigBias150_opch22[num];


    for(int i=0; i<num; ++i){
        bias50_opch01[i] = 2 * (mpv50_opch01[i] - mpv99_opch01[i]) / (mpv50_opch01[i] + mpv99_opch01[i]);
        bias50_opch03[i] = 2 * (mpv50_opch03[i] - mpv99_opch03[i]) / (mpv50_opch03[i] + mpv99_opch03[i]);
        bias50_opch16[i] = 2 * (mpv50_opch16[i] - mpv99_opch16[i]) / (mpv50_opch16[i] + mpv99_opch16[i]);
        bias50_opch22[i] = 2 * (mpv50_opch22[i] - mpv99_opch22[i]) / (mpv50_opch22[i] + mpv99_opch22[i]);

        bias70_opch01[i] = 2 * (mpv70_opch01[i] - mpv99_opch01[i]) / (mpv70_opch01[i] + mpv99_opch01[i]);
        bias70_opch03[i] = 2 * (mpv70_opch03[i] - mpv99_opch03[i]) / (mpv70_opch03[i] + mpv99_opch03[i]);
        bias70_opch16[i] = 2 * (mpv70_opch16[i] - mpv99_opch16[i]) / (mpv70_opch16[i] + mpv99_opch16[i]);
        bias70_opch22[i] = 2 * (mpv70_opch22[i] - mpv99_opch22[i]) / (mpv70_opch22[i] + mpv99_opch22[i]);

        bias130_opch01[i] = 2 * (mpv130_opch01[i] - mpv99_opch01[i]) / (mpv130_opch01[i] + mpv99_opch01[i]);
        bias130_opch03[i] = 2 * (mpv130_opch03[i] - mpv99_opch03[i]) / (mpv130_opch03[i] + mpv99_opch03[i]);
        bias130_opch16[i] = 2 * (mpv130_opch16[i] - mpv99_opch16[i]) / (mpv130_opch16[i] + mpv99_opch16[i]);
        bias130_opch22[i] = 2 * (mpv130_opch22[i] - mpv99_opch22[i]) / (mpv130_opch22[i] + mpv99_opch22[i]);

        bias150_opch01[i] = 2 * (mpv150_opch01[i] - mpv99_opch01[i]) / (mpv150_opch01[i] + mpv99_opch01[i]);
        bias150_opch03[i] = 2 * (mpv150_opch03[i] - mpv99_opch03[i]) / (mpv150_opch03[i] + mpv99_opch03[i]);
        bias150_opch16[i] = 2 * (mpv150_opch16[i] - mpv99_opch16[i]) / (mpv150_opch16[i] + mpv99_opch16[i]);
        bias150_opch22[i] = 2 * (mpv150_opch22[i] - mpv99_opch22[i]) / (mpv150_opch22[i] + mpv99_opch22[i]);

        sigBias50_opch01[i] = sigBias(mpv50_opch01[i], sigMPV50_opch01[i], mpv99_opch01[i], sigMPV99_opch01[i]);   
        sigBias50_opch03[i] = sigBias(mpv50_opch03[i], sigMPV50_opch03[i], mpv99_opch03[i], sigMPV99_opch03[i]); 
        sigBias50_opch16[i] = sigBias(mpv50_opch16[i], sigMPV50_opch16[i], mpv99_opch16[i], sigMPV99_opch16[i]); 
        sigBias50_opch22[i] = sigBias(mpv50_opch22[i], sigMPV50_opch22[i], mpv99_opch22[i], sigMPV99_opch22[i]); 

        sigBias70_opch01[i] = sigBias(mpv70_opch01[i], sigMPV70_opch01[i], mpv99_opch01[i], sigMPV99_opch01[i]);   
        sigBias70_opch03[i] = sigBias(mpv70_opch03[i], sigMPV70_opch03[i], mpv99_opch03[i], sigMPV99_opch03[i]); 
        sigBias70_opch16[i] = sigBias(mpv70_opch16[i], sigMPV70_opch16[i], mpv99_opch16[i], sigMPV99_opch16[i]); 
        sigBias70_opch22[i] = sigBias(mpv70_opch22[i], sigMPV70_opch22[i], mpv99_opch22[i], sigMPV99_opch22[i]); 

        sigBias130_opch01[i] = sigBias(mpv130_opch01[i], sigMPV130_opch01[i], mpv99_opch01[i], sigMPV99_opch01[i]);   
        sigBias130_opch03[i] = sigBias(mpv130_opch03[i], sigMPV130_opch03[i], mpv99_opch03[i], sigMPV99_opch03[i]); 
        sigBias130_opch16[i] = sigBias(mpv130_opch16[i], sigMPV130_opch16[i], mpv99_opch16[i], sigMPV99_opch16[i]); 
        sigBias130_opch22[i] = sigBias(mpv130_opch22[i], sigMPV130_opch22[i], mpv99_opch22[i], sigMPV99_opch22[i]); 

        sigBias150_opch01[i] = sigBias(mpv150_opch01[i], sigMPV150_opch01[i], mpv99_opch01[i], sigMPV99_opch01[i]);   
        sigBias150_opch03[i] = sigBias(mpv150_opch03[i], sigMPV150_opch03[i], mpv99_opch03[i], sigMPV99_opch03[i]); 
        sigBias150_opch16[i] = sigBias(mpv150_opch16[i], sigMPV150_opch16[i], mpv99_opch16[i], sigMPV99_opch16[i]); 
        sigBias150_opch22[i] = sigBias(mpv150_opch22[i], sigMPV150_opch22[i], mpv99_opch22[i], sigMPV99_opch22[i]); 

    }


    //bias50_opch01-----------------------------------------------------------------------
    TGraphErrors* graphBias50_opch01 = new TGraphErrors(num, dis50_opch01, bias50_opch01, nullptr, sigBias50_opch01);
    graphBias50_opch01->SetLineColor(kGreen);
    graphBias50_opch01->SetMarkerStyle(20);
    graphBias50_opch01->SetMarkerSize(0.6);
    graphBias50_opch01->SetMarkerColor(kGreen);
//    graphBias50_opch01->Draw("AP");
    graphBias50_opch01->GetXaxis()->SetTitle("Distance [cm]");
    graphBias50_opch01->GetYaxis()->SetTitle("2*(RSLxx - RSL100)/(RSLxx + RSL100)");
    graphBias50_opch01->GetXaxis()->SetRangeUser(0, 300);
//    graphBias50_opch01->GetYaxis()->SetRangeUser(-1.2, 1.2);

    //diff50_opch03---
    TGraphErrors* graphBias50_opch03 = new TGraphErrors(num, dis50_opch03, bias50_opch03, nullptr, sigBias50_opch03);
    graphBias50_opch03->SetLineColor(kGreen);
    graphBias50_opch03->SetMarkerStyle(21);
    graphBias50_opch03->SetMarkerSize(0.6);
    graphBias50_opch03->SetMarkerColor(kGreen);
//    graphBias50_opch03->Draw("P SAME");


/*    graphBias50_opch03->Draw("AP");
    graphBias50_opch03->GetXaxis()->SetTitle("Distance [cm]");
    graphBias50_opch03->GetYaxis()->SetTitle("2*(RSLxx - RSL100)/(RSLxx + RSL100)");
    graphBias50_opch03->GetXaxis()->SetRangeUser(0, 600);
    graphBias50_opch03->GetYaxis()->SetRangeUser(-1.2, 1.2);
*/
    //diff50_opch16---
    TGraphErrors* graphBias50_opch16 = new TGraphErrors(num, dis50_opch16, bias50_opch16, nullptr, sigBias50_opch16);
    graphBias50_opch16->SetLineColor(kGreen);
    graphBias50_opch16->SetMarkerStyle(44);
    graphBias50_opch16->SetMarkerSize(0.7);
    graphBias50_opch16->SetMarkerColor(kGreen);
//    graphBias50_opch16->Draw("P SAME");


    //diff50_opch22---
    TGraphErrors* graphBias50_opch22 = new TGraphErrors(num, dis50_opch22, bias50_opch22, nullptr, sigBias50_opch22);
    graphBias50_opch22->SetLineColor(kGreen);
    graphBias50_opch22->SetMarkerStyle(22);
    graphBias50_opch22->SetMarkerSize(0.6);
    graphBias50_opch22->SetMarkerColor(kGreen);
//    graphBias50_opch22->Draw("P SAME");
    graphBias50_opch22->Draw("AP");

    legend3->AddEntry(graphBias50_opch01, "RSL50", "pe");


    //bias70_opch01-----------------------------------------------------------------------
    TGraphErrors* graphBias70_opch01 = new TGraphErrors(num, dis70_opch01, bias70_opch01, nullptr, sigBias70_opch01);
    graphBias70_opch01->SetLineColor(kBlue);
    graphBias70_opch01->SetMarkerStyle(20);
    graphBias70_opch01->SetMarkerSize(0.6);
    graphBias70_opch01->SetMarkerColor(kBlue);
//    graphBias70_opch01->Draw("P SAME");

    //diff70_opch03---
    TGraphErrors* graphBias70_opch03 = new TGraphErrors(num, dis70_opch03, bias70_opch03, nullptr, sigBias70_opch03);
    graphBias70_opch03->SetLineColor(kBlue);
    graphBias70_opch03->SetMarkerStyle(21);
    graphBias70_opch03->SetMarkerSize(0.6);
    graphBias70_opch03->SetMarkerColor(kBlue);
//    graphBias70_opch03->Draw("P SAME");

    //diff70_opch16---
    TGraphErrors* graphBias70_opch16 = new TGraphErrors(num, dis70_opch16, bias70_opch16, nullptr, sigBias70_opch16);
    graphBias70_opch16->SetLineColor(kBlue);
    graphBias70_opch16->SetMarkerStyle(44);
    graphBias70_opch16->SetMarkerSize(0.7);
    graphBias70_opch16->SetMarkerColor(kBlue);
//    graphBias70_opch16->Draw("P SAME");

    //diff70_opch22---
    TGraphErrors* graphBias70_opch22 = new TGraphErrors(num, dis70_opch22, bias70_opch22, nullptr, sigBias70_opch22);
    graphBias70_opch22->SetLineColor(kBlue);
    graphBias70_opch22->SetMarkerStyle(22);
    graphBias70_opch22->SetMarkerSize(0.6);
    graphBias70_opch22->SetMarkerColor(kBlue);
    graphBias70_opch22->Draw("P SAME");

    legend3->AddEntry(graphBias70_opch01, "RSL70", "pe");


    //bias130_opch01-----------------------------------------------------------------------
    TGraphErrors* graphBias130_opch01 = new TGraphErrors(num, dis130_opch01, bias130_opch01, nullptr, sigBias130_opch01);
    graphBias130_opch01->SetLineColor(28);
    graphBias130_opch01->SetMarkerStyle(20);
    graphBias130_opch01->SetMarkerSize(0.6);
    graphBias130_opch01->SetMarkerColor(28);
//    graphBias130_opch01->Draw("P SAME");

    //diff130_opch03---
    TGraphErrors* graphBias130_opch03 = new TGraphErrors(num, dis130_opch03, bias130_opch03, nullptr, sigBias130_opch03);
    graphBias130_opch03->SetLineColor(28);
    graphBias130_opch03->SetMarkerStyle(21);
    graphBias130_opch03->SetMarkerSize(0.6);
    graphBias130_opch03->SetMarkerColor(28);
//    graphBias130_opch03->Draw("P SAME");

    //diff130_opch16---
    TGraphErrors* graphBias130_opch16 = new TGraphErrors(num, dis130_opch16, bias130_opch16, nullptr, sigBias130_opch16);
    graphBias130_opch16->SetLineColor(28);
    graphBias130_opch16->SetMarkerStyle(44);
    graphBias130_opch16->SetMarkerSize(0.7);
    graphBias130_opch16->SetMarkerColor(28);
//    graphBias130_opch16->Draw("P SAME");

    //diff130_opch22---
    TGraphErrors* graphBias130_opch22 = new TGraphErrors(num, dis130_opch22, bias130_opch22, nullptr, sigBias130_opch22);
    graphBias130_opch22->SetLineColor(28);
    graphBias130_opch22->SetMarkerStyle(22);
    graphBias130_opch22->SetMarkerSize(0.6);
    graphBias130_opch22->SetMarkerColor(28);
    graphBias130_opch22->Draw("P SAME");

    legend3->AddEntry(graphBias130_opch01, "RSL130", "pe");


    //bias150_opch01-----------------------------------------------------------------------
    TGraphErrors* graphBias150_opch01 = new TGraphErrors(num, dis150_opch01, bias150_opch01, nullptr, sigBias150_opch01);
    graphBias150_opch01->SetLineColor(kViolet);
    graphBias150_opch01->SetMarkerStyle(20);
    graphBias150_opch01->SetMarkerSize(0.6);
    graphBias150_opch01->SetMarkerColor(kViolet);
//    graphBias150_opch01->Draw("P SAME");

    //diff150_opch03---
    TGraphErrors* graphBias150_opch03 = new TGraphErrors(num, dis150_opch03, bias150_opch03, nullptr, sigBias150_opch03);
    graphBias150_opch03->SetLineColor(kViolet);
    graphBias150_opch03->SetMarkerStyle(21);
    graphBias150_opch03->SetMarkerSize(0.6);
    graphBias150_opch03->SetMarkerColor(kViolet);
//    graphBias150_opch03->Draw("P SAME");

    //diff150_opch16---
    TGraphErrors* graphBias150_opch16 = new TGraphErrors(num, dis150_opch16, bias150_opch16, nullptr, sigBias150_opch16);
    graphBias150_opch16->SetLineColor(kViolet);
    graphBias150_opch16->SetMarkerStyle(44);
    graphBias150_opch16->SetMarkerSize(0.7);
    graphBias150_opch16->SetMarkerColor(kViolet);
//    graphBias150_opch16->Draw("P SAME");

    //diff150_opch22---
    TGraphErrors* graphBias150_opch22 = new TGraphErrors(num, dis150_opch22, bias150_opch22, nullptr, sigBias150_opch22);
    graphBias150_opch22->SetLineColor(kViolet);
    graphBias150_opch22->SetMarkerStyle(22);
    graphBias150_opch22->SetMarkerSize(0.6);
    graphBias150_opch22->SetMarkerColor(kViolet);
    graphBias150_opch22->Draw("P SAME");

    legend3->AddEntry(graphBias150_opch01, "RSL150", "pe");




    legend3->Draw();


//Create root file==================================================================
    TFile* outputFile = new TFile(TString(path)+"/compare.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file for writing" << std::endl;
        return;
    }


    canvas.Write();
    canvas2.Write();
    canvas3.Write();

    outputFile->Close();

    //TEST------------------
/*    for(int i=0; i<60; ++i){
        std::cout<<"mpv: "<<mpv99[i]<<std::endl;
    }
*/
}

    
//==========Main Function===============================================
void rsl_compare_m1_singleOpCh(){

    rsl_compareTMP("/Users/shuaixiangzhang/Work/current/FNAL_Work2024/rsl_analyses/v4_analysis/results/fit_Develop/membrane1");
        

}