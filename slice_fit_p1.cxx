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




//Used to store fitting variables----------
struct FitResults{
    double meanG; //mean value of Gaussian function
    double widthG; //sigma of Gaussain
    double mpvL; //mpv of Landau function
    double sigmaL; //width of Landau
    double mpvC; //mpv of convolution function
    double sigC; //sigma of convolution function
};



//==========================================================
//Convoluted Landau + Gaussian fitting---
//For Gaussain, mean is not zero!
FitResults CLG1(TH1F* hist, TFile* outputFile, double fitRangeMin, double fitRangeMax) {
    //Create a RooRealVar for the variable you are fitting
    RooRealVar x("x", "Variable", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    //Create a RooDataHist from the TH1F
    RooDataHist data("data", "Data Histogram", x, RooFit::Import(*hist));

    //Create RooRealVar for Landau parameters
    RooRealVar mpv("mpv", "Most Probable Value", 3, 0.1, 50);
    RooRealVar sigma("sigma", "Sigma", 1, 0.01, 20);
    //Create RooLandau PDF
    RooLandau landau("landau", "Landau PDF", x, mpv, sigma);

    //Create RooRealVar for Gaussian parameters
    RooRealVar mean("mean", "Mean", 3, 1, 50);
    RooRealVar width("width", "Width", 1, 0.01, 20);
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
    double xMin = fitRangeMin;
    double xMax = fitRangeMax;
    double xStep = 0.1;  //Adjust the step size as needed
    double maxVal = -1.0;
    double xPeak = -1.0; //mpv of conv function---
    for (double xPos = xMin; xPos <= xMax; xPos += xStep) {
        x.setVal(xPos);
        double val = convoluted.getVal(); 
        if (val > maxVal) {
            maxVal = val;
            xPeak = xPos;
        }
    }

    //To find FHWM------
    double halfL = 0.0;
    double halfR = 0.0;
    double fwhm = 0.0;
    for (double xPos = xMin; xPos <= xPeak; xPos += xStep) {
        x.setVal(xPos);
        double val = convoluted.getVal(); 
        if (val > 0.5*maxVal) {
            halfL = xPos;
            break;
        }
    }
    for (double xPos = xPeak; xPos <= xMax; xPos += xStep) {
        x.setVal(xPos);
        double val = convoluted.getVal(); 
        if (val < 0.5*maxVal) {
            halfR = xPos;
            break;
        }
    }
    fwhm = halfR - halfL;

    double sigConv = 0.0;
    sigConv = fwhm / 2.355;
    std::cout<<"\n\nPeak of conv: "<<xPeak<<";   FWHM: "<<fwhm<<";   sigConv: "<<sigConv<<"\n\n"<<std::endl;




    //store fitting results---
    FitResults results;

    //--------------------------------------------------------------

    //clone the original histogram---
    TH1F* histClone = (TH1F*)hist->Clone("histClone");

    //Plot the result within the specified range
    RooPlot* frame = x.frame(RooFit::Range(fitRangeMin, fitRangeMax));
    data.plotOn(frame, RooFit::DataError(RooAbsData::SumW2));
    convoluted.plotOn(frame, RooFit::Range(fitRangeMin, fitRangeMax));


/*
    //evaluate fitting, chi square method--
    double chi2 = frame->chiSquare();//chi-squared value
    int ndf = convoluted.getParameters(data)->selectByAttrib("Constant", kFALSE)->getSize(); // number of degrees of freedom
*/



    //Create a canvas and draw the frame
    TCanvas canvas("canvas", "Convoluted Landau + Gaussian Fit");
    histClone->Draw();
    frame->Draw("SAME");//Draw fitting function---



    //legend---------
    TLegend legend(0.5, 0.45, 0.9, 0.9);
    legend.SetTextSize(0.03);
    legend.AddEntry(histClone, "Original Histogram", "l");
    legend.AddEntry(frame->getObject(0), "Conv Landau + Gaussian", "l");
    //Add fitting parameters to the legend
    RooLinkedListIter iter = convoluted.getParameters(data)->iterator();
    RooRealVar* var;
    while ((var = dynamic_cast<RooRealVar*>(iter.Next()))) {
        legend.AddEntry((TObject*)0, Form("%s = %.2f", var->GetName(), var->getVal()), "");

        //store fitting results to struct---------
        const char* varName = var->GetName();
        double varValue = var->getVal();
        if (strcmp(varName, "mean") == 0) {
                results.meanG = varValue;
        } 
        else if (strcmp(varName, "width") == 0) {
                results.widthG = varValue;
        } 
        else if (strcmp(varName, "mpv") == 0) {
                results.mpvL = varValue;
        } 
        else if (strcmp(varName, "sigma") == 0) {
                results.sigmaL = varValue;
        }


    }

    //Add the MPV of the convoluted function
    legend.AddEntry((TObject*)0, Form("MPV of Conv = %.2f", xPeak), "");
    //Add the MPV of the convoluted function
    legend.AddEntry((TObject*)0, Form("#sigma of Conv = %.2f", sigConv), "");









    std::cout<<"\n\nMPV of convvv is: "<<xPeak<<std::endl;

    //Calculate chi2 and ndf manually----------------------------------
    double chi2 = 0.0;
    int ndf = 0;

    for (int bin = 1; bin <= histClone->GetNbinsX(); ++bin) {
        double xBinCenter = histClone->GetXaxis()->GetBinCenter(bin);
        if (xBinCenter < fitRangeMin || xBinCenter > fitRangeMax) {
            continue;  // Skip points outside the fitting range
        }

        double observed = histClone->GetBinContent(bin);
        RooArgSet argSet(x);// Create a RooArgSet for the variable 'x'
        x.setVal(xBinCenter);
    
//        double expected = convoluted.getValV(&argSet);
        double expected = convoluted.getValV();

        if (expected > 0.0) {
            double error = sqrt(histClone->GetBinContent(bin));
        //    double error = histClone->GetBinError(bin);
            if (error > 0.0) {
            double residual = (observed - expected) / error;
            chi2 += residual * residual;

            std::cout<<"\nobserve: "<<observed<<"   expect: "<<expected<<"   error: "<<error<<"\n"<<std::endl;

            ndf++;
            }
        }
    }

    // Retrieve the number of floating parameters from the RooFitResult
    int numFloatingParams = fitResult->floatParsFinal().getSize();

    //Minus the number of fitting parameters and add 2 empty bins--
    ndf = ndf - numFloatingParams + 2;



    //Add chi-squared information to the legend
//    legend.AddEntry((TObject*)0, Form("#chi^{2} / ndf = %.2f / %d", chi2, ndf), "");









    results.mpvC = xPeak; 
    results.sigC = sigConv;
    legend.Draw("SAME");

    canvas.Write();//Save canvas into the ROOT file---
//    frame->Write();//save fitting function---

    return results;

}
//----------------------------------------------------------








//==========================================================
//Convoluted Landau + Gaussian fitting---
//Set mean of Gaussian as zero---
FitResults CLG2(TH1F* hist, TFile* outputFile, double fitRangeMin, double fitRangeMax) {
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
    RooRealVar mean("mean", "Mean", 0.0, 0.0, 0.0);
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
    double xMin = fitRangeMin;
    double xMax = fitRangeMax;
    double xStep = 0.01;  //Adjust the step size as needed
    double maxVal = -1.0;
    double xPeak = -1.0; //mpv of conv function---
    for (double xPos = xMin; xPos <= xMax; xPos += xStep) {
        x.setVal(xPos);
        double val = convoluted.getVal(); 
        if (val > maxVal) {
            maxVal = val;
            xPeak = xPos;
        }
    }

    //To find FHWM------
    double halfL = 0.0;
    double halfR = 0.0;
    double fwhm = 0.0;
    for (double xPos = xMin; xPos <= xPeak; xPos += xStep) {
        x.setVal(xPos);
        double val = convoluted.getVal(); 
        if (val > 0.5*maxVal) {
            halfL = xPos;
            break;
        }
    }
    for (double xPos = xPeak; xPos <= xMax; xPos += xStep) {
        x.setVal(xPos);
        double val = convoluted.getVal(); 
        if (val < 0.5*maxVal) {
            halfR = xPos;
            break;
        }
    }
    fwhm = halfR - halfL;

    double sigConv = 0.0;
    sigConv = fwhm / 2.355;
    std::cout<<"\n\nPeak of conv: "<<xPeak<<";   FWHM: "<<fwhm<<";   sigConv: "<<sigConv<<"\n\n"<<std::endl;

    //store fitting results---
    FitResults results;


    //clone the original histogram---
    TH1F* histClone = (TH1F*)hist->Clone("histClone");

    //Plot the result within the specified range
    RooPlot* frame = x.frame(RooFit::Range(fitRangeMin, fitRangeMax));
    data.plotOn(frame, RooFit::DataError(RooAbsData::SumW2));
    convoluted.plotOn(frame, RooFit::Range(fitRangeMin, fitRangeMax));



    //evaluate fitting, chi square method--
    double chi2 = frame->chiSquare();//chi-squared value
    int ndf = convoluted.getParameters(data)->selectByAttrib("Constant", kFALSE)->getSize(); // number of degrees of freedom


    //Create a canvas and draw the frame
    TCanvas canvas("canvas", "Convoluted Landau + Gaussian Fit");
    histClone->Draw();
    frame->Draw("SAME");//Draw fitting function---



    //legend---------
    TLegend legend(0.5, 0.45, 0.9, 0.9);
    legend.SetTextSize(0.03);
    legend.AddEntry(histClone, "Original Histogram", "l");
    legend.AddEntry(frame->getObject(0), "Convoluted Landau + Gaussian Fit", "l");
    //Add fitting parameters to the legend
    RooLinkedListIter iter = convoluted.getParameters(data)->iterator();
    RooRealVar* var;
    while ((var = dynamic_cast<RooRealVar*>(iter.Next()))) {
        legend.AddEntry((TObject*)0, Form("%s = %.2f", var->GetName(), var->getVal()), "");

        //store fitting results to struct---------
        const char* varName = var->GetName();
        double varValue = var->getVal();
        if (strcmp(varName, "mean") == 0) {
                results.meanG = varValue;
        } 
        else if (strcmp(varName, "width") == 0) {
                results.widthG = varValue;
        } 
        else if (strcmp(varName, "mpv") == 0) {
                results.mpvL = varValue;
        } 
        else if (strcmp(varName, "sigma") == 0) {
                results.sigmaL = varValue;
        }


    }

    //Add the MPV of the convoluted function
    legend.AddEntry((TObject*)0, Form("MPV of Conv = %.2f", xPeak), "");
    //Add the MPV of the convoluted function
    legend.AddEntry((TObject*)0, Form("#sigma of Conv = %.2f", sigConv), "");

    //Add chi-squared information to the legend
//    legend.AddEntry((TObject*)0, Form("#chi^{2} / ndf = %.2f / %d", chi2, ndf), "");


    results.mpvC = xPeak; 
    results.sigC = sigConv;
    legend.Draw("SAME");

    canvas.Write();//Save canvas into the ROOT file---
//    frame->Write();//save fitting function---

    return results;

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
//ouput the most right bin is not zero entry
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
//Draw mpv values on initial TH2F plot------
void DrawScatterWithLine(TH2F* hist2D, const double* distances, const double* mpvConv, const double* sigConv, int nPoints, TFile* outputFile) {
    //Create a TGraph for scatter points
    //TGraph* scatterGraph = new TGraph(nPoints, distances, mpvConv);

    //Create a TGraphErrors object
    TGraphErrors* scatterGraph = new TGraphErrors(nPoints, distances, mpvConv, nullptr, sigConv);
    scatterGraph->SetLineColor(kRed);


    // Create a TGraph for linear connections
//    TGraph* lineGraph = new TGraph(nPoints, distances, mpvConv);
//    lineGraph->SetLineColor(kRed);  // Set line color to blue
//    lineGraph->SetLineWidth(1);      // Set line width

    // Create a canvas
    TCanvas canvas("mpvOnTh2F", "Fitting results with MPV values", 800, 600);

    // Draw the 2D histogram
    hist2D->Draw("colz");

    // Draw scatter points
    scatterGraph->SetMarkerStyle(20);  // Set marker style to a filled circle
    scatterGraph->SetMarkerSize(0.6);  // Set marker size
    scatterGraph->SetMarkerColor(kRed);
    scatterGraph->Draw("P");

    // Draw linear connections
//    lineGraph->Draw("L");

    // Save the canvas to the output ROOT file
    canvas.Write();

    // Clean up memory
    delete scatterGraph;
//    delete lineGraph;
}
//----------------------------------------------------------











//==========================================================
//Main function!---
void slice_fit_p1(){
    //change file name each time-----------------------
    string file_path = "../results/combine_2000results/";
    string file_suffix = "rsl99_2000num_e67_crtCut.root";
    string output_path = "../results/fit_Develop/pmt1/";
    string output_name = "fitCLG1";

    //store fitting results at txt file---
    std::ofstream outputTxt("../results/fit_Develop/pmt1/rsl99_fit.txt");

    //Choose the slice you want to look at!---
    //Define the X(distance) values where you want to extract data---
    Double_t distances[60];
    for (int i=0; i<60; ++i){
        distances[i] = 5*i + 2.5;
    }
    //--------------------------------------------------

    FitResults tmpResults; //used to tmporarily store fitResults---
    Double_t mpvConv[60];
    Double_t sigConv[60];


    //Open the input root file----------------------------------------------------
    TFile* inputFile = new TFile(TString(file_path)+TString(file_suffix), "READ");

    //Check if the file is open---
    if(!inputFile || inputFile->IsZombie()){
        std::cerr<<"Error: Cannot open input.root"<<std::endl;
        return;
    }

    //Access the TH2F from the file---
    TH2F* inputTH2F = (TH2F*)inputFile->Get("summedPMT1");

    if(!inputTH2F){
        std::cerr<<"Error: Cannot find TH2F in the input file"<< std::endl;
        inputFile->Close();
        return;
    }

    //Create an array of TH1F histograms to store extracted data---
    //num = 60 here!---
    const int num = sizeof(distances)/sizeof(distances[0]);
    TH1F* hists[num];


    //Extract values from TH2F and create TH1F---
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

    }

    //Open the output ROOT file for writing---
    TFile* outputFile = new TFile(TString(output_path)+TString(output_name)+"_"+TString(file_suffix), "RECREATE");

    double borderL = 0.0;
    double borderR = 0.0;

    //draw each histo in single canvas---
    for(int i=0; i<num; ++i){
        adjustXAxisBinWidth(hists[i], 1000.0);

        if(distances[i] <= 80){//hist has 100 bins
            combineBins(hists[i]);
            
            if(distances[i] <= 50){//hist has 50 bins
                combineBins(hists[i]);
            }
        }

        borderR = rightBorder(hists[i]);
        borderL = leftBorder(hists[i]); 

        tmpResults = CLG1(hists[i], outputFile, borderL, borderR);
//        tmpResults = CLG2(hists[i], outputFile, borderL, borderR);
        mpvConv[i] = tmpResults.mpvC * 1000.0;
        sigConv[i] = tmpResults.sigC * 1000.0;

//        hists[i]->Write();
    }
    
    //Draw mpvConv on TH2F---
    DrawScatterWithLine(inputTH2F, distances, mpvConv, sigConv, num, outputFile);


    //Create a subdirectory within the output file
    TDirectory* histDirectory = outputFile->mkdir("Initial_Hist");
    histDirectory->cd();
    for(int i=0; i<num; ++i){
        hists[i]->Write();
    }
    
    //Close the input and output files---
    outputFile->Close();
    inputFile->Close();



    //Store fitting files---------
    if (!outputTxt.is_open()) {
        std::cerr << "Error: Cannot open output file for writing" << std::endl;
        return;
    }
    //Write distances, mpvConv and sigConv
    for (int i = 0; i < num; ++i) {
        outputTxt << distances[i] << "\t" << mpvConv[i] << "\t" << sigConv[i] << std::endl;
    }
    outputTxt.close();

}
//----------------------------------------------------------




