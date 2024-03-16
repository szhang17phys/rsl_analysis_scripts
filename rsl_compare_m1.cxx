#include <iostream>

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>


void readFile(const std::string& fileName, double distances[], double mpvConv[], double sigConv[], int maxArraySize) {
    //Open the input file for reading
    std::ifstream inputFile(fileName);

    //Check if the file is open
    if(!inputFile.is_open()) {
        std::cerr << "Error: Cannot open input file '" << fileName << "' for reading" << std::endl;
        return;
    }

    //Temporary variables to store data from each line
    double distance, mpv, sig;
    int index = 0;
    
    //Read data from each line and store them into arrays
    while (inputFile >> distance >> mpv >> sig && index < maxArraySize) {
        distances[index] = distance;
        mpvConv[index] = mpv;
        sigConv[index] = sig;
        index++;

        if(distance > 150){
            break;
        }
    }

    //Close the input file
    inputFile.close();
}






//======MAIN FUNCTION=======================================
void rsl_compare_m1(){

    const int num = 30; //num of slices---

    std::string inputFile99 = "../results/fit_Develop/membrane1/rsl99_fit.txt";
    double dis99[num];
    double mpv99[num];
    double sig99[num];
    readFile(inputFile99, dis99, mpv99, sig99, num);

    std::string inputFile70 = "../results/fit_Develop/membrane1/rsl70_fit.txt";
    double dis70[num];
    double mpv70[num];
    double sig70[num];
    readFile(inputFile70, dis70, mpv70, sig70, num);

    std::string inputFile50 = "../results/fit_Develop/membrane1/rsl50_fit.txt";
    double dis50[num];
    double mpv50[num];
    double sig50[num];
    readFile(inputFile50, dis50, mpv50, sig50, num);

    std::string inputFile150 = "../results/fit_Develop/membrane1/rsl150_fit.txt";
    double dis150[num];
    double mpv150[num];
    double sig150[num];
    readFile(inputFile150, dis150, mpv150, sig150, num);



    //Drawing--------------------------------
    TCanvas canvas("fitCompare", "Different RSLs", 800, 600);
    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);

    //RSL99
    TGraphErrors* scatterGraph = new TGraphErrors(num, dis99, mpv99, nullptr, sig99);
    scatterGraph->SetLineColor(kRed);
    scatterGraph->SetMarkerStyle(21);//20: filled circle
    scatterGraph->SetMarkerSize(0.6);
    scatterGraph->SetMarkerColor(kRed);
    scatterGraph->Draw("AP");
    scatterGraph->GetXaxis()->SetTitle("Distance [cm]");
    scatterGraph->GetYaxis()->SetTitle("#photon / event");
    scatterGraph->GetXaxis()->SetRangeUser(0, 300);
    scatterGraph->GetYaxis()->SetRangeUser(0, 200000);
    canvas.Update();
//    TGraph* lineGraph = new TGraph(num, dis99, mpv99);
//    lineGraph->SetLineColor(kRed);
//    lineGraph->SetLineWidth(1);
//    lineGraph->Draw("L");


    //RSL70
    TGraph* scatterGraph70 = new TGraph(num, dis70, mpv70);
    scatterGraph70->SetLineColor(kBlue);
    scatterGraph70->SetMarkerStyle(21);//20: filled square
    scatterGraph70->SetMarkerSize(0.6);
    scatterGraph70->SetMarkerColor(kBlue);
    scatterGraph70->Draw("P SAME");


    //RSL50
    TGraph* scatterGraph50 = new TGraph(num, dis50, mpv50);
    scatterGraph50->SetLineColor(kGreen);
    scatterGraph50->SetMarkerStyle(21);//20: filled square
    scatterGraph50->SetMarkerSize(0.6);
    scatterGraph50->SetMarkerColor(kGreen);
    scatterGraph50->Draw("P SAME");


    //RSL150
    TGraph* scatterGraph150 = new TGraph(num, dis150, mpv150);
    scatterGraph150->SetLineColor(kOrange);
    scatterGraph150->SetMarkerStyle(21);//20: filled square
    scatterGraph150->SetMarkerSize(0.6);
    scatterGraph150->SetMarkerColor(kOrange);
    scatterGraph150->Draw("P SAME");

    
    legend->AddEntry(scatterGraph150, "RSL = 150.0cm", "p");
    legend->AddEntry(scatterGraph, "RSL = 99.9cm", "pe");
    legend->AddEntry(scatterGraph70, "RSL = 70.0cm", "p");
    legend->AddEntry(scatterGraph50, "RSL = 50.0cm", "p");


    //Create root file---------------
    TFile* outputFile = new TFile("../results/fit_Develop/membrane1/fitResults_compare.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file for writing" << std::endl;
        return;
    }

    legend->Draw();
    canvas.Write();

    outputFile->Close();


    //Clean up memory----------
    delete scatterGraph;
//    delete lineGraph;
    delete scatterGraph70;
//    delete lineGraph70;
    delete scatterGraph50;


    //TEST------------------
/*    for(int i=0; i<60; ++i){
        std::cout<<"mpv: "<<mpv99[i]<<std::endl;
    }
*/
}

    

