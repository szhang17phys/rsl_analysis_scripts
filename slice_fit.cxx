#include <iostream>

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooLandau.h>
#include <RooFitResult.h>

//#include <RooNumIntOptions.h>

#include <RooGaussian.h>
#include "RooDataHist.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooLinkedListIter.h"

#include "TMinuit.h"



//Pay attention:
//Here is cathodeXA responses (change fitting initialization)---





//==========================================================
//Convoluted Landau + Gaussian fitting---
void CLG(TH1F* hist, TFile* outputFile, double fitRangeMin, double fitRangeMax) {
    //Create a RooRealVar for the variable you are fitting
    RooRealVar x("x", "Variable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    //Create a RooDataHist from the TH1F
    RooDataHist data("data", "Data Histogram", x, RooFit::Import(*hist));

    //Create RooRealVar for Landau parameters
    RooRealVar mpv("mpv", "Most Probable Value", 10, 1, 200);
    RooRealVar sigma("sigma", "Sigma", 5, 0.01, 100);
    //Create RooLandau PDF
    RooLandau landau("landau", "Landau PDF", x, mpv, sigma);

    //Create RooRealVar for Gaussian parameters
    RooRealVar mean("mean", "Mean", 10, 1, 200);
    RooRealVar width("width", "Width", 5, 0.01, 100);
    RooGaussian gaussian("gaussian", "Gaussian PDF", x, mean, width);

    //Create RooFFTConvPdf for the convoluted function
    RooFFTConvPdf convoluted("convoluted", "Convoluted Landau + Gaussian", x, landau, gaussian);

    //Create a RooCmdArg for the fit range
    RooCmdArg fitRangeArg = RooFit::Range(fitRangeMin, fitRangeMax);

    //Perform the fit, the default is maximum likelihood method--
    RooFitResult* fitResult = convoluted.fitTo(data, RooFit::Save(true), fitRangeArg );
//    RooFitResult* fitResult = convoluted.fitTo(data, RooFit::Save(true), RooFit::Minimizer("Minuit2", "Migrad"));

    //Access fit results---
    fitResult->Print("v");

    //To find the position of the peak of convoluted curve----------
    // To find the position of the peak of the convoluted curve
    double xMin = fitRangeMin;
    double xMax = fitRangeMax;
    double xStep = 0.1;  //Adjust the step size as needed
    double maxVal = -1.0;
    double xPeak = -1.0;
    for (double xPos = xMin; xPos <= xMax; xPos += xStep) {
        x.setVal(xPos);
        double val = convoluted.getVal(); 
        if (val > maxVal) {
            maxVal = val;
            xPeak = xPos;
        }
    }
    //--------------------------------------------------------------

    //clone the original histogram---
    TH1F* histClone = (TH1F*)hist->Clone("histClone");

    //Plot the result within the specified range
    RooPlot* frame = x.frame(RooFit::Range(fitRangeMin, fitRangeMax));
    data.plotOn(frame, RooFit::DataError(RooAbsData::SumW2));
    convoluted.plotOn(frame, RooFit::Range(fitRangeMin, fitRangeMax));

    //Create a canvas and draw the frame
    TCanvas canvas("canvas", "Convoluted Landau + Gaussian Fit");
    histClone->Draw();
    frame->Draw("SAME");//Draw fitting function---

    //evaluate fitting, chi square method--
    double chi2 = frame->chiSquare();//chi-squared value
    int ndf = convoluted.getParameters(data)->selectByAttrib("Constant", kFALSE)->getSize(); // number of degrees of freedom
    double chi2ndf = chi2 / ndf; // normalized chi-squared

    //legend---
    TLegend legend(0.6, 0.6, 0.9, 0.9);
    legend.AddEntry(histClone, "Original Histogram", "l");
    legend.AddEntry(frame->getObject(0), "Convoluted Landau + Gaussian Fit", "l");
    //Add fitting parameters to the legend
    RooLinkedListIter iter = convoluted.getParameters(data)->iterator();
    RooRealVar* var;
    while ((var = dynamic_cast<RooRealVar*>(iter.Next()))) {
        legend.AddEntry((TObject*)0, Form("%s = %.2f", var->GetName(), var->getVal()), "");
    }

    //Add the MPV of the convoluted function
    legend.AddEntry((TObject*)0, Form("MPV of Conv = %.2f", xPeak), "");

    //Add chi-squared information to the legend
    legend.AddEntry((TObject*)0, Form("#chi^2/ndf = %.2f", chi2ndf), "");

    legend.Draw("SAME");

    canvas.Write();//Save canvas into the ROOT file---
//    frame->Write();//save fitting function---

}
//----------------------------------------------------------









//==========================================================
//Perform Landau fit (with maximum likelihood estimation)---
void landauFit(TH1F* hist, TFile* outputFile){
    //Create a RooRealVar for the variable you are fitting---
    RooRealVar x("x", "Variable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    //Create a RooDataHist from the TH1F---
    RooDataHist data("data", "Data Histogram", x, RooFit::Import(*hist));

    //Create a Landau PDF---
    RooRealVar mpv("mpv", "Most Probable Value", 10, 1, 200);
    RooRealVar sigma("sigma", "Sigma", 5, 0.1, 100);
    RooLandau landau("landau", "Landau PDF", x, mpv, sigma);

    //Perform the fit using maximum likelihood---
    RooFitResult* fitResult = landau.fitTo(data, RooFit::Save(true));

    //Access fit results---
    fitResult->Print("v");

    //clone the original histogram---
    TH1F* histClone = (TH1F*)hist->Clone("histClone"); 

    //Plot the result---
    //SumW2: draw sqrt(content) error bars---
    RooPlot* frame = x.frame();
    data.plotOn(frame, RooFit::DataError(RooAbsData::SumW2));
    landau.plotOn(frame);

    TCanvas canvas("canvas", "Landau Fit with ML");
    histClone->Draw();
    frame->Draw("SAME");//Draw fitting function--- 

    TLegend legend(0.6, 0.6, 0.9, 0.9);
    legend.AddEntry(histClone, "Original Histogram", "l");
    legend.AddEntry(frame->getObject(0), "Landau Fit", "l");
    //Add fitting parameters to the legend
    RooLinkedListIter iter = landau.getParameters(data)->iterator();
    RooRealVar* var;
    while ((var = dynamic_cast<RooRealVar*>(iter.Next()))) {
        legend.AddEntry((TObject*)0, Form("%s = %.4f", var->GetName(), var->getVal()), "");
    }

    legend.Draw("SAME");

    canvas.Write();//Save canvas into the ROOT file---
    frame->Write();//save fitting function---
    //Save the fit parameters (if needed)
//    outputFile->cd();
//    fitResult->Write(("fitParameters_" + std::string(hist->GetName())).c_str());

}
//----------------------------------------------------------








//==========================================================
//Landau + Gaussian fit (with mle)---
void landauGaussianFit(TH1F* hist, TFile* outputFile) {
    //Create a RooRealVar for the variable you are fitting
    RooRealVar x("x", "Variable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    //Create a RooDataHist from the TH1F
    RooDataHist data("data", "Data Histogram", x, RooFit::Import(*hist));

    //Create RooRealVar for Landau parameters
    RooRealVar mpv("mpv", "Most Probable Value", 10, 1, 200);
    RooRealVar sigma("sigma", "Sigma", 5, 0.1, 100);
    //Create RooLandau PDF
    RooLandau landau("landau", "Landau PDF", x, mpv, sigma);

    //Create RooRealVar for Gaussian parameters
    RooRealVar mean("mean", "Mean", 10, 1, 200);
    RooRealVar width("width", "Width", 5, 0.1, 100);
    // Create RooGaussian PDF
    RooGaussian gaussian("gaussian", "Gaussian PDF", x, mean, width);

    //Create RooRealVar for the mixture coefficient
    RooRealVar landauFraction("landauFraction", "Landau Fraction", 0.001, 0, 0.001);

    //Create RooAddPdf to combine Landau and Gaussian
    RooAddPdf model("model", "Landau + Gaussian", RooArgList(landau, gaussian), RooArgList(landauFraction));

    //Perform the fit using maximum likelihood
    RooFitResult* fitResult = model.fitTo(data, RooFit::Save(true));

    //Perform the fit using least squares 
//    RooNumIntOptions::defaultMinimizerType(RooNumIntOptions::kLevenbergMarquardt);
//    RooFitResult* fitResult = model.fitTo(data, RooFit::Save(true), RooFit::Minimizer("Minuit2", "Migrad"));


    //clone the original histogram---
    TH1F* histClone = (TH1F*)hist->Clone("histClone");

    //Plot the result
    RooPlot* frame = x.frame();
    data.plotOn(frame, RooFit::DataError(RooAbsData::SumW2));
    model.plotOn(frame);

    //Create a canvas and draw the frame
    TCanvas canvas("canvas", "Landau + Gaussian Fit");
    histClone->Draw();
    frame->Draw("SAME");//Draw fitting function---

    //evaluate fitting, chi square method--
    double chi2 = frame->chiSquare();//chi-squared value
    int ndf = model.getParameters(data)->selectByAttrib("Constant", kFALSE)->getSize(); // number of degrees of freedom
    double chi2ndf = chi2 / ndf; // normalized chi-squared

    //legend---
    TLegend legend(0.6, 0.6, 0.9, 0.9);
    legend.AddEntry(histClone, "Original Histogram", "l");
    legend.AddEntry(frame->getObject(0), "Landau + Gaussian Fit", "l");
    //Add fitting parameters to the legend
    RooLinkedListIter iter = model.getParameters(data)->iterator();
    RooRealVar* var;
    while ((var = dynamic_cast<RooRealVar*>(iter.Next()))) {
        legend.AddEntry((TObject*)0, Form("%s = %.4f", var->GetName(), var->getVal()), "");
    }

    //Add chi-squared information to the legend
    legend.AddEntry((TObject*)0, Form("Chi^2/ndf = %.3f", chi2ndf), "");

    legend.Draw("SAME");

    canvas.Write();//Save canvas into the ROOT file---
//    frame->Write();//save fitting function---
}
//----------------------------------------------------------









//==========================================================
//The reason of resetting bin width is that the initial 
//value is too large,making it hard to do Landau fitting---
void adjustXAxisBinWidth(TH1F* hist, double scale) {
    //Clone the original histogram
    TH1F* histClone = dynamic_cast<TH1F*>(hist->Clone("histClone"));
    
    //Clear the original histogram
    hist->Reset();

    //Set the new binning to the original histogram
    int numBins = histClone->GetNbinsX();
    double xMin = histClone->GetXaxis()->GetXmin();
    double xMax = histClone->GetXaxis()->GetXmax();
    double binWidth = histClone->GetXaxis()->GetBinWidth(1);  
    double newWidth = static_cast<double>(binWidth / scale);

    hist->SetBins(numBins, xMin, xMin + numBins * newWidth);

    //Fill new hist with the content of the original hist
    for (int i = 1; i <= numBins; ++i) {
        double x = histClone->GetXaxis()->GetBinCenter(i);
        double content = histClone->GetBinContent(i);

        hist->SetBinContent(i, content);
    }

    delete histClone;
}
//----------------------------------------------------------






//==========================================================
//combine adjacent bins of histograms---
void combineBins(TH1F* originalHist) {
    //Create a temporary histogram with half the number of bins
    int newNumBins = originalHist->GetNbinsX() / 2;
    TH1F* combinedHist = new TH1F("combinedHist", "Combined Histogram", newNumBins, originalHist->GetXaxis()->GetXmin(), originalHist->GetXaxis()->GetXmax());

    //Fill new hist by summing two adjacent bins from original hist
    for (int i = 1; i <= newNumBins; ++i) {
        int binIndexOriginal = 2 * i - 1; //Consider odd bins of original hist
        double binContent = originalHist->GetBinContent(binIndexOriginal) + originalHist->GetBinContent(binIndexOriginal + 1);
        combinedHist->SetBinContent(i, binContent);
    }

    //Copy the content of combined hist back to original hist
    originalHist->Reset(); //Clear the original histogram
    originalHist->Add(combinedHist);

    //Delete the temporary histogram to avoid memory leaks
    delete combinedHist;
}
//----------------------------------------------------------










//==========================================================
//output the position of right border of the hist
double rightBorder(TH1F* hist){
    //Find the bin with the maximum entry
    int binMax = hist->GetMaximumBin();

    // Get the maximum entry
    double maxEntry = hist->GetBinContent(binMax);

    //Find the first bin with zero entry after the maximum
    int binWithZero = binMax + 1;//Start from the bin after the maximum
    //two adjacent bins have zero entries
    while(binWithZero < hist->GetNbinsX() && hist->GetBinContent(binWithZero) > 0){
        binWithZero++;
    }
    if(hist->GetBinContent(binWithZero+1) >=2){
        binWithZero++;//aovid empty bin by fluctuation
        while(binWithZero < hist->GetNbinsX() && hist->GetBinContent(binWithZero) >= 1){
            binWithZero++;
        }
    }

    //Check if the distance to the last bin is less than 5 bins
    int lastBin = hist->GetNbinsX();
    if ((lastBin - binWithZero) < 1) {
        //Output the position of the last bin
        double lastBinPosition = hist->GetBinCenter(lastBin);
        return lastBinPosition;
    }

    //Output the position of the final bin as a double
    double position = hist->GetBinCenter(binWithZero);

    return position;
}
//----------------------------------------------------------








//==========================================================
//output the position of left border of the hist
double leftBorder(TH1F* hist) {
    //Find the bin with the maximum entry
    int binMax = hist->GetMaximumBin();

    //Find the first bin with zero entry before the maximum
    int binWithZero = binMax - 1;  //Start from the bin before the maximum
    while(binWithZero > 1 && hist->GetBinContent(binWithZero) > 0){
        binWithZero--;
    }
    if(hist->GetBinContent(binWithZero-1) >=2){
        binWithZero--;//aovid empty bin by fluctuation
        while(binWithZero > 1 && hist->GetBinContent(binWithZero) >= 1){
            binWithZero--;
        }
    }

    //Check if the distance to the first bin is less than 5 bins
    if((binMax - binWithZero) < 1){
        // Output the position of the first bin
        double firstBinPosition = hist->GetBinCenter(1);
        return firstBinPosition;
    }

    //Output the position of the final bin as a double
    double position = hist->GetBinCenter(binWithZero);

    return position;
}
//----------------------------------------------------------







//==========================================================
//Main function!---
void slice_fit(){
    //change file name each time-----------------------
    string file_path = "../results/fit_Develop/";
    string file_suffix = "rsl99_2000num_e67_crtCut.root";
    string output_path = "../results/fit_Develop/";
    string output_name = "fit";

    //Choose the slice you want to look at!---
    //Define the X(distance) values where you want to extract data---
    Double_t distances[60];
    for (int i=0; i<60; ++i){
        distances[i] = 5*i + 2.5;
    }
    //-------------------------------------------

    //Open the input root file---
    TFile* inputFile = new TFile(TString(file_path)+TString(file_suffix), "READ");

    //Check if the file is open---
    if(!inputFile || inputFile->IsZombie()){
        std::cerr<<"Error: Cannot open input.root"<<std::endl;
        return;
    }

    //Access the TH2F from the file---
    TH2F* inputTH2F = (TH2F*)inputFile->Get("summedCathode8XA");

    if(!inputTH2F){
        std::cerr<<"Error: Cannot find TH2F in the input file"<< std::endl;
        inputFile->Close();
        return;
    }

    //Create an array of TH1F histograms to store extracted data---
    const int num = sizeof(distances)/sizeof(distances[0]);
    TH1F* hists[num];

    //Draw all TH1Fs in one canvas---
    TCanvas* canvas = new TCanvas("canvas", "All Histograms", 800, 600);
    //Add a legend(optional)---
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    //Extract values for each distance value and create TH1F histograms---
    for(int i=0; i<num; ++i){
        double x = distances[i];
        int binX = inputTH2F->GetXaxis()->FindBin(x);

        //Create a TH1F for each extracted data---
        hists[i] = new TH1F(Form("hist_%d", i), Form("Distance = %.1f cm", x), inputTH2F->GetNbinsY(), inputTH2F->GetYaxis()->GetXmin(), inputTH2F->GetYaxis()->GetXmax());

        hists[i]->GetXaxis()->SetTitle("# #gamma / 1000");
        hists[i]->GetYaxis()->SetTitle("Event Rate");

        //Fill the TH1F with data from the TH2F---
        for(int binY=1; binY<=inputTH2F->GetNbinsY(); ++binY){
            double y = inputTH2F->GetYaxis()->GetBinCenter(binY);
            double value = inputTH2F->GetBinContent(binX, binY);
            hists[i]->SetBinContent(binY, value);
        }

        //Set the line color for each TH1F(optional)---
        hists[i]->SetLineColor(i + 1);

        //Draw each TH1F on the canvas---
        if(i==0){
            hists[i]->Draw();
        } 
        else{
            hists[i]->Draw("SAME");
        }

        legend->AddEntry(hists[i], Form("Dis = %.1f cm", distances[i]), "l");

    }

    //Open the output ROOT file for writing---
    TFile* outputFile = new TFile(TString(output_path)+TString(output_name)+"_"+TString(file_suffix), "RECREATE");

    //plot containing several histos---
    legend->Draw("SAME");
    canvas->Write();
    legend->Write();

    double borderL = 0.0;
    double borderR = 0.0;

    //draw each histo in single canvas---
    for(int i=0; i<num; ++i){

        adjustXAxisBinWidth(hists[i], 1000.0);

        if(distances[i] <= 70){//hist has 100 bins
            combineBins(hists[i]);
        }
        if(distances[i] <= 40){//hist has 50 bins
            combineBins(hists[i]);
        }

//        landauFit(hists[i], outputFile);//apply Landau fitting---
//        landauGaussianFit(hists[i], outputFile);//Landau + Gaussian---
        borderR = rightBorder(hists[i]);
        borderL = leftBorder(hists[i]);        
        CLG(hists[i], outputFile, borderL, borderR);

        hists[i]->Write();
    }

    //Close the input and output files---
    inputFile->Close();
    outputFile->Close();

}
//----------------------------------------------------------




