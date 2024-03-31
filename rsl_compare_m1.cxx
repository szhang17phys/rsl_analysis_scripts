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

    std::string inputFile99 = path + "/rsl99_fitCLG1.txt";
    std::string inputFile50 = path + "/rsl50_fitCLG1.txt";
    std::string inputFile70 = path + "/rsl70_fitCLG1.txt";
    std::string inputFile100 = path + "/rsl100_fitCLG1.txt";
    std::string inputFile150 = path + "/rsl150_fitCLG1.txt";


    double dis99[num];
    double mpv99[num];
    double sig99[num];
    int ent99[num];//entries---
    readFile(inputFile99, dis99, mpv99, sig99, ent99, num);
    double sigMPV99[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent99[i]>0){
            sigMPV99[i] = sig99[i]/sqrt(static_cast<double>(ent99[i]));
        }
        else sigMPV99[i] = 0.0;
    }


    double dis50[num];
    double mpv50[num];
    double sig50[num];
    int ent50[num];
    readFile(inputFile50, dis50, mpv50, sig50, ent50, num);
    double sigMPV50[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent50[i]>0){
            sigMPV50[i] = sig50[i]/sqrt(static_cast<double>(ent50[i]));
        }
        else sigMPV50[i] = 0.0;
    }


    double dis70[num];
    double mpv70[num];
    double sig70[num];
    int ent70[num];
    readFile(inputFile70, dis70, mpv70, sig70, ent70, num);
    double sigMPV70[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent70[i]>0){
            sigMPV70[i] = sig70[i]/sqrt(static_cast<double>(ent70[i]));
        }
        else sigMPV70[i] = 0.0;
    }


    double dis100[num];
    double mpv100[num];
    double sig100[num];
    int ent100[num];//entries---
    readFile(inputFile100, dis100, mpv100, sig100, ent100, num);
    double sigMPV100[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent100[i]>0){
            sigMPV100[i] = sig100[i]/sqrt(static_cast<double>(ent100[i]));
        }
        else sigMPV100[i] = 0.0;
    }


    double dis150[num];
    double mpv150[num];
    double sig150[num];
    int ent150[num];
    readFile(inputFile150, dis150, mpv150, sig150, ent150, num);
    double sigMPV150[num];//sigma of MPV---
    for(int i=0; i<num; ++i){
        if(ent150[i]>0){
            sigMPV150[i] = sig150[i]/sqrt(static_cast<double>(ent150[i]));
        }
        else sigMPV150[i] = 0.0;
    }


/*    double minimum = 100;
    for(int i=0; i<60; ++i){
        std::cout<<"sigMPV: "<<sigMPV50[i]<<endl;
        std::cout<<"sigMPV: "<<sigMPV70[i]<<endl;
        std::cout<<"sigMPV: "<<sigMPV150[i]<<endl; 
        if (sigMPV50[i] <= minimum && sigMPV50[i]>0){
            minimum = sigMPV50[i];
        }  
        if (sigMPV70[i] <= minimum && sigMPV70[i]>0){
            minimum = sigMPV70[i];
        } 
        if (sigMPV150[i] <= minimum && sigMPV150[i]){
            minimum = sigMPV150[i];
        }                 

    }
    std::cout<<"\nMinimum sigMPV: "<<minimum<<std::endl;
*/


    //Drawing 1-------------------------------------------------------------------
    TCanvas canvas("fitCompare", "Different RSLs", 800, 600);
    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);

    //RSL99
    TGraphErrors* scatterGraph = new TGraphErrors(num, dis99, mpv99, nullptr, sigMPV99);
    scatterGraph->SetLineColor(kRed);
    scatterGraph->SetMarkerStyle(21);//20: filled circle
    scatterGraph->SetMarkerSize(0.6);
    scatterGraph->SetMarkerColor(kRed);
    scatterGraph->Draw("AP");
    scatterGraph->GetXaxis()->SetTitle("Distance [cm]");
    scatterGraph->GetYaxis()->SetTitle("#photon / event");
    scatterGraph->GetXaxis()->SetRangeUser(0, 300);
    scatterGraph->GetYaxis()->SetRangeUser(0, 200000);
//    TGraph* lineGraph = new TGraph(num, dis99, mpv99);
//    lineGraph->SetLineColor(kRed);
//    lineGraph->SetLineWidth(1);
//    lineGraph->Draw("L");

    //RSL50
//    TGraph* scatterGraph50 = new TGraph(num, dis50, mpv50);
    TGraphErrors* scatterGraph50 = new TGraphErrors(num, dis50, mpv50, nullptr, sigMPV50);
    scatterGraph50->SetLineColor(kGreen);
    scatterGraph50->SetMarkerStyle(21);//20: filled square
    scatterGraph50->SetMarkerSize(0.6);
    scatterGraph50->SetMarkerColor(kGreen);
    scatterGraph50->Draw("P SAME");

    //RSL70
//    TGraph* scatterGraph50 = new TGraph(num, dis50, mpv50);
    TGraphErrors* scatterGraph70 = new TGraphErrors(num, dis70, mpv70, nullptr, sigMPV70);
    scatterGraph70->SetLineColor(kBlue);
    scatterGraph70->SetMarkerStyle(21);//20: filled square
    scatterGraph70->SetMarkerSize(0.6);
    scatterGraph70->SetMarkerColor(kBlue);
    scatterGraph70->Draw("P SAME");

    //RSL100
//    TGraph* scatterGraph70 = new TGraph(num, dis70, mpv70);
    TGraphErrors* scatterGraph100 = new TGraphErrors(num, dis100, mpv100, nullptr, sigMPV100);
    scatterGraph100->SetLineColor(kViolet);
    scatterGraph100->SetMarkerStyle(21);//20: filled square
    scatterGraph100->SetMarkerSize(0.6);
    scatterGraph100->SetMarkerColor(kViolet);
    scatterGraph100->Draw("P SAME");


    //RSL150
//    TGraph* scatterGraph150 = new TGraph(num, dis150, mpv150);
    TGraphErrors* scatterGraph150 = new TGraphErrors(num, dis150, mpv150, nullptr, sigMPV150);
    scatterGraph150->SetLineColor(kOrange);
    scatterGraph150->SetMarkerStyle(21);//20: filled square
    scatterGraph150->SetMarkerSize(0.6);
    scatterGraph150->SetMarkerColor(kOrange);
    scatterGraph150->Draw("P SAME");
    
    legend->AddEntry(scatterGraph150, "RSL = 150.0cm", "pe");
    legend->AddEntry(scatterGraph100, "RSL = 100.0cm", "pe");
    legend->AddEntry(scatterGraph, "RSL = 99.9cm", "pe");
    legend->AddEntry(scatterGraph70, "RSL = 70.0cm", "pe");
    legend->AddEntry(scatterGraph50, "RSL = 50.0cm", "pe");
    legend->Draw();




//Create root file--------------------------------------------------------
    TFile* outputFile = new TFile(TString(path)+"/fitCLG1_compare.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file for writing" << std::endl;
        return;
    }




//Drawing 2-------------------------------------------------------------------
    TCanvas canvas2("fitDiff", "Different RSLs", 800, 600);
    TLegend* legend2 = new TLegend(0.6, 0.6, 0.9, 0.9);

    double diff50[num];//(mpv50 - mpv99)---
    double diff70[num];
    double diff100[num];
    double diff150[num];
    double sigDiff50[num];//sigma of diff50---
    double sigDiff70[num];
    double sigDiff100[num];    
    double sigDiff150[num];    

    for(int i=0; i<num; ++i){
        diff50[i] = mpv50[i] - mpv99[i];
        diff70[i] = mpv70[i] - mpv99[i];
        diff100[i] = mpv100[i] - mpv99[i];
        diff150[i] = mpv150[i] - mpv99[i];

        sigDiff50[i] = sqrt(sigMPV50[i]*sigMPV50[i] + sigMPV99[i]*sigMPV99[i]);
        sigDiff70[i] = sqrt(sigMPV70[i]*sigMPV70[i] + sigMPV99[i]*sigMPV99[i]);
        sigDiff100[i] = sqrt(sigMPV100[i]*sigMPV100[i] + sigMPV99[i]*sigMPV99[i]);       
        sigDiff150[i] = sqrt(sigMPV150[i]*sigMPV150[i] + sigMPV99[i]*sigMPV99[i]);  

    }

    //diff50
    TGraphErrors* graphDiff50 = new TGraphErrors(num, dis50, diff50, nullptr, sigDiff50);
    graphDiff50->SetLineColor(kGreen);
    graphDiff50->SetMarkerStyle(20);//20: filled circle
    graphDiff50->SetMarkerSize(0.6);
    graphDiff50->SetMarkerColor(kGreen);
    graphDiff50->Draw("AP");
    graphDiff50->GetXaxis()->SetTitle("Distance [cm]");
    graphDiff50->GetYaxis()->SetTitle("Photon Num Diff");
    graphDiff50->GetXaxis()->SetRangeUser(0, 300);
    graphDiff50->GetYaxis()->SetRangeUser(-10000, 40000);

    //diff70
    TGraphErrors* graphDiff70 = new TGraphErrors(num, dis70, diff70, nullptr, sigDiff70);
    graphDiff70->SetLineColor(kBlue);
    graphDiff70->SetMarkerStyle(20);//20: filled square
    graphDiff70->SetMarkerSize(0.6);
    graphDiff70->SetMarkerColor(kBlue);
    graphDiff70->Draw("P SAME");

    //diff100
    TGraphErrors* graphDiff100 = new TGraphErrors(num, dis100, diff100, nullptr, sigDiff100);
    graphDiff100->SetLineColor(kViolet);
    graphDiff100->SetMarkerStyle(20);//20: filled square
    graphDiff100->SetMarkerSize(0.6);
    graphDiff100->SetMarkerColor(kViolet);
    graphDiff100->Draw("P SAME");


    //diff150
    TGraphErrors* graphDiff150 = new TGraphErrors(num, dis150, diff150, nullptr, sigDiff150);
    graphDiff150->SetLineColor(kOrange);
    graphDiff150->SetMarkerStyle(20);//20: filled square
    graphDiff150->SetMarkerSize(0.6);
    graphDiff150->SetMarkerColor(kOrange);
    graphDiff150->Draw("P SAME");
    
    legend2->AddEntry(graphDiff50, "RSL50 - RSL100", "pe");
    legend2->AddEntry(graphDiff70, "RSL70 - RSL99", "pe");
    legend2->AddEntry(graphDiff100, "RSL100 - RSL99", "pe");
    legend2->AddEntry(graphDiff150, "RSL150 - RSL100", "pe");
    legend2->Draw();






//Drawing 3-------------------------------------------------------------------
    TCanvas canvas3("fitBias", "fit Bias", 800, 600);
    TLegend* legend3 = new TLegend(0.6, 0.6, 0.9, 0.9);

    double bias50[num];//2*(mpv50-mpv99)/(mpv50+mpv99)---
    double bias70[num];
    double bias100[num];
    double bias150[num];
    double sigBias50[num];//sigma of bias50---
    double sigBias70[num];
    double sigBias100[num];
    double sigBias150[num];  


    for(int i=0; i<num; ++i){
        bias50[i] = 2 * (mpv50[i] - mpv99[i]) / (mpv50[i] + mpv99[i]);
        bias70[i] = 2 * (mpv70[i] - mpv99[i]) / (mpv70[i] + mpv99[i]);
        bias100[i] = 2 * (mpv100[i] - mpv99[i]) / (mpv100[i] + mpv99[i]);
        bias150[i] = 2 * (mpv150[i] - mpv99[i]) / (mpv150[i] + mpv99[i]);

        sigBias50[i] = sigBias(mpv50[i], sigMPV50[i], mpv99[i], sigMPV99[i]);   
        sigBias70[i] = sigBias(mpv70[i], sigMPV70[i], mpv99[i], sigMPV99[i]);
        sigBias100[i] = sigBias(mpv100[i], sigMPV100[i], mpv99[i], sigMPV99[i]);
        sigBias150[i] = sigBias(mpv150[i], sigMPV150[i], mpv99[i], sigMPV99[i]);

    }

    //bias50
    TGraphErrors* graphBias50 = new TGraphErrors(num, dis50, bias50, nullptr, sigBias50);
    graphBias50->SetLineColor(kGreen);
    graphBias50->SetMarkerStyle(20);//20: filled circle
    graphBias50->SetMarkerSize(0.6);
    graphBias50->SetMarkerColor(kGreen);
//    graphBias50->Draw("AP");
    graphBias50->GetXaxis()->SetTitle("Distance [cm]");
    graphBias50->GetYaxis()->SetTitle("2*(RSLxx - RSL100)/(RSLxx + RSL100)");
    graphBias50->GetXaxis()->SetRangeUser(0, 300);
//    graphBias50->GetYaxis()->SetRangeUser(-10000, 40000);

    //diff70
    TGraphErrors* graphBias70 = new TGraphErrors(num, dis70, bias70, nullptr, sigBias70);
    graphBias70->SetLineColor(kBlue);
    graphBias70->SetMarkerStyle(20);//20: filled square
    graphBias70->SetMarkerSize(0.6);
    graphBias70->SetMarkerColor(kBlue);
//    graphBias70->Draw("P SAME");

    graphBias70->Draw("AP");
    graphBias70->GetXaxis()->SetTitle("Distance [cm]");
    graphBias70->GetYaxis()->SetTitle("2*(RSLxx - RSL100)/(RSLxx + RSL100)");
    graphBias70->GetXaxis()->SetRangeUser(0, 300);

    //diff100
    TGraphErrors* graphBias100 = new TGraphErrors(num, dis100, bias100, nullptr, sigBias100);
    graphBias100->SetLineColor(kViolet);
    graphBias100->SetMarkerStyle(20);//20: filled square
    graphBias100->SetMarkerSize(0.6);
    graphBias100->SetMarkerColor(kViolet);
    graphBias100->Draw("P SAME");

    //diff150
    TGraphErrors* graphBias150 = new TGraphErrors(num, dis150, bias150, nullptr, sigBias150);
    graphBias150->SetLineColor(kOrange);
    graphBias150->SetMarkerStyle(20);//20: filled square
    graphBias150->SetMarkerSize(0.6);
    graphBias150->SetMarkerColor(kOrange);
//    graphBias150->Draw("P SAME");
    
//    legend3->AddEntry(graphBias50, "RSL50", "pe");
    legend3->AddEntry(graphBias70, "RSL70", "pe");
    legend3->AddEntry(graphBias100, "RSL100", "pe");
//    legend3->AddEntry(graphBias150, "RSL150", "pe");
    legend3->Draw();





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
void rsl_compare_m1(){

    //For membrane1---
//    rsl_compareTMP("../results/fit_Develop/membrane1");
        

}

    

