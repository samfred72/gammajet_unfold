#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TF1.h"
#include <vector> 
#include <cmath>
#include "style.h"

double FindlargestValue(std::vector<TH1D*> histDetails){
    double maxVal = 0.0;
    for(TH1D* h : histDetails){
        double currentMax = h->GetMaximum();
        if (currentMax > maxVal) {
            maxVal = currentMax;
        }
    }
    return maxVal;
}

void singleGraph(TH1D* histDetails,const char* fileName  , TDirectory* directory, bool print, bool logg, bool EMCAL){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 700, 700);
    if(logg) c1->SetLogy();
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    
    histDetails->GetXaxis()->SetNdivisions(510);
    histDetails->GetYaxis()->SetNdivisions(510);
    histDetails->GetYaxis()->SetTitle("Ratio #frac{p_{T,Data}}{p_{T,MC}}");
    
    histDetails->SetLineColor(kRed+1);
    
    if(EMCAL) histDetails->Fit("pol2");
    else histDetails->Fit("pol1");
    histDetails->Draw("PE");

    
    TLegend* leg = new TLegend(0.65,0.80,0.85,0.85);
    leg->AddEntry(histDetails, "Data","l");
    leg->Draw();
    
    directory->cd();
    c1->Write("canvasRatio");
    histDetails->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.png", fileName)); 
    c1->Close();
}

void doubleGraph(TH1D* histDetails,TH1D* PythiaHIST,const char* fileName , const char* newTitle ,TDirectory* directory, bool print, bool logg ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 700, 700);
    if(logg) c1->SetLogy();
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    
    histDetails->GetXaxis()->SetNdivisions(510);
    histDetails->GetYaxis()->SetNdivisions(510);
    PythiaHIST->GetYaxis()->SetRangeUser(0,PythiaHIST->GetMaximum() * 1.5);

    
    histDetails->SetLineColor(kRed+1);
    PythiaHIST->SetLineColor(kBlue+1);
    
    PythiaHIST->Draw("HIST");
    histDetails->Draw("HIST SAME");
    
    TLegend* leg = new TLegend(0.65,0.80,0.85,0.85);
    leg->AddEntry(histDetails, "Data","l");
    leg->AddEntry(PythiaHIST, "Pythia","l");
    leg->Draw();

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.23, 0.85, newTitle);
    
    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.png", fileName)); 
    c1->Close();
}   

void doubleGraphERatio(std::vector<TH1D*> histDetails,std::vector<TH1D*> PythiaHIST,const char* fileName , const char* newTitle ,TDirectory* directory, bool print, bool logg ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 700, 700);
    if(logg) c1->SetLogy();
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    
    PythiaHIST[0]->GetXaxis()->SetNdivisions(510);
    PythiaHIST[0]->GetYaxis()->SetNdivisions(510);
    PythiaHIST[0]->GetYaxis()->SetRangeUser(0,FindlargestValue(histDetails) * 1.3);

    
    histDetails[0]->SetLineColor(kRed+1);
    histDetails[0]->SetMarkerColor(kRed+1);
    histDetails[0]->SetMarkerStyle(20);

    histDetails[1]->SetLineColor(kCyan+2);
    histDetails[1]->SetMarkerColor(kCyan+2);
    histDetails[1]->SetMarkerStyle(20);

    PythiaHIST[0]->SetLineColor(kRed+1);
    PythiaHIST[0]->SetMarkerColor(kRed+1);
    PythiaHIST[0]->SetMarkerStyle(24);

    PythiaHIST[1]->SetLineColor(kCyan+2);
    PythiaHIST[1]->SetMarkerColor(kCyan+2);
    PythiaHIST[1]->SetMarkerStyle(24);
    
    PythiaHIST[0]->Draw("E1");
    PythiaHIST[1]->Draw("E1 SAME");
    histDetails[0]->Draw("E1 SAME");
    histDetails[1]->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.55,0.65,0.75,0.85);
    leg->AddEntry(histDetails[0], "Leading Data","lep");
    leg->AddEntry(histDetails[1], "SubLeading Data","lep");
    leg->AddEntry(PythiaHIST[0], "Leading Pythia","lep");
    leg->AddEntry(PythiaHIST[1], "SubLeading Pythia","lep");
    leg->Draw();

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.23, 0.85, newTitle);
    
    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.png", fileName)); 
    c1->Close();
}   

void doubleGraphTripple(std::vector<TH1D*> histDetails,std::vector<TH1D*>  PythiaHIST,const char* fileName , const char* LatexLabel ,std::vector<const char*> newTitle ,TDirectory* directory, bool print, bool logg ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 1600, 400);
    c1->Divide(4,1);
    int index= 0;
    for(int i = 0 ; i < 4 ; i++){
        c1->cd(i + 1);
        if(i == 0){
            TLatex latex1;
            latex1.SetNDC();
            latex1.SetTextSize(0.06);
            latex1.DrawLatex(0.42, 0.85, LatexLabel);

            continue;
        }
        if(logg) c1->SetLogy();
        c1->SetTicks(1,1);
        c1->SetLeftMargin(0.15);
        
        histDetails[index]->GetXaxis()->SetNdivisions(510);
        histDetails[index]->GetYaxis()->SetNdivisions(510);
        PythiaHIST[index]->GetYaxis()->SetRangeUser(0,histDetails[index]->GetMaximum() * 1.3);

        
        histDetails[index]->SetLineColor(kRed+1);
        PythiaHIST[index]->SetLineColor(kBlue+1);
        
        PythiaHIST[index]->Draw("HIST");
        histDetails[index]->Draw("HIST SAME");
        
        TLegend* leg = new TLegend(0.65,0.80,0.85,0.85);
        leg->AddEntry(histDetails[index], "Data","l");
        leg->AddEntry(PythiaHIST[index], "Pythia","l");
        leg->Draw();

        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.03);
        latex1.DrawLatex(0.23, 0.85, newTitle[index]);
        index++;
    }    
    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.png", fileName)); 
    c1->Close();
}   

void DanPlots(TH1D* histDetails,const char* fileName , const char* newTitle , const char* LatexLabel,TDirectory* directory, bool print, bool logg ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 700, 700);
    if(logg) c1->SetLogy();
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    
    histDetails->SetTitle(newTitle);
    histDetails->GetXaxis()->SetNdivisions(510);
    histDetails->GetYaxis()->SetNdivisions(510);

    
    histDetails->SetLineColorAlpha(kAzure+1, 0.7);histDetails->SetLineWidth(2);
    
    histDetails->Draw("HIST");

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.28, 0.75,LatexLabel);
    
    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.png", fileName)); 
    c1->Close();
}   

void xJGraph(TH1D* histDetails,const char* fileName , const char* LLabel , TDirectory* directory, bool print,bool logg ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 900, 900);
    if(logg) c1->SetLogy();
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    
    
    histDetails->GetXaxis()->SetNdivisions(510);
    histDetails->GetYaxis()->SetNdivisions(510);

    
    histDetails->SetLineColor(kRed+1);
    
    histDetails->Draw("HIST");
    
    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.20, 0.77,Form(
        "#splitline{#Delta#phi_{1,2} > 7#pi/8}"
        "{#splitline{p_{T,1} > 30 GeV}"
        "{x_j =%s }}", LLabel));
    
    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.png", fileName)); 
    c1->Close();
}   

void singleGraph2D(TH2D* histDetails,const char* fileName , const char* LatexLabel , TDirectory* directory,bool print, bool logg ){
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 900, 900);
    if(logg) c1->SetLogz();
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    c1->SetTopMargin(0.22);
    
    histDetails->GetXaxis()->SetNdivisions(510);
    histDetails->GetYaxis()->SetNdivisions(510);
    
    histDetails->Draw("COLZ");

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.42, 0.95, LatexLabel);
    
    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    c1->Close();
}

void trippletGraph2D(std::vector<TH2D*> histDetails,const char* fileName , const char* LatexLabel , std::vector<const char*> range, TDirectory* directory, bool print, bool logg ){
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 1800, 400);
    c1->Divide(4);
    int index = 0;
    for(int i = 0 ; i < 4 ; i++){
        c1->cd(i + 1); 
        if(i == 0 ){
            TLatex latex1;
            latex1.SetNDC();
            latex1.SetTextSize(0.06);
            latex1.DrawLatex(0.42, 0.85, LatexLabel);

            continue;
        }
        gStyle->SetPaintTextFormat("4.3f"); // Set to 3 decimal places
        if(logg) c1->SetLogz();
        c1->SetTicks(1,1);
        gPad->SetRightMargin(0.17);
        histDetails[index]->GetXaxis()->SetNdivisions(510);
        histDetails[index]->GetYaxis()->SetNdivisions(510);
        
        histDetails[index]->SetMarkerColor(kWhite); // Good for dark COLZ palettes
        histDetails[index]->SetMarkerSize(1.2);     // Shrink text further if 3 decimals
        histDetails[index]->Draw("COLZ TEXT");
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.05);
        latex1.DrawLatex(0.28, 0.85, range[index]);
        index++;
    }
    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    c1->Close();
}

void singleGraph3D(TH3D* histDetails,const char* fileName , TDirectory* directory, bool print, bool logg ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 600, 600);
    if(logg) {c1->SetLogz(); c1->SetLogy();}
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    
    histDetails->GetXaxis()->SetNdivisions(510);
    histDetails->GetYaxis()->SetNdivisions(510);
    histDetails->GetZaxis()->SetNdivisions(510);

    histDetails->Draw("BOX3");
    
    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.png", fileName)); 
    c1->Close();

}

void CompGraphxJ(std::vector<TH1D*> histDetails,const char* fileName , const char* newTitle ,const char* Llabel1,const char* Llabel2,const char* LatexLabel ,TDirectory* directory, bool print, bool logg ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 900, 900);
    
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);

    /*
    TF1* f1 = new TF1("f1", "gaus", 0.05, 1.0);
    histDetails[0]->Fit(f1, "R");   

    TF1* f2 = new TF1("f2", "gaus", 0.05, 1.2);
    histDetails[1]->Fit(f2, "R");   
    */

    if(logg) c1->SetLogy();
    histDetails[0]->SetTitle(newTitle);
    histDetails[0]->GetXaxis()->SetNdivisions(510);
    histDetails[0]->GetYaxis()->SetNdivisions(510);
    histDetails[0]->GetYaxis()->SetRangeUser(0,histDetails[1]->GetMaximum()*1.2);

    
    histDetails[0]->SetLineColorAlpha(kRed+1,0.7); histDetails[0]->SetLineWidth(2);
    histDetails[1]->SetLineColorAlpha(kBlue+1,0.7); histDetails[1]->SetLineWidth(2);
    /*
    f1->SetLineColorAlpha(kRed+1,0.5);
    f2->SetLineColorAlpha(kBlue+1,0.5);
    */
    
    histDetails[0]->Draw("HIST E1");
    histDetails[1]->Draw("HIST E1 SAME");
    

    TLegend* leg = new TLegend(0.2,0.7,0.35,0.85);
    leg->SetBorderSize(0);
    leg->AddEntry(histDetails[0], Llabel1,"l");
    leg->AddEntry(histDetails[1], Llabel2,"l");
    leg->Draw();
    
    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.64, 0.83,LatexLabel);
    double mean1 = histDetails[0]->GetMean();double sigma1 = histDetails[0]->GetMeanError();
    double mean2 = histDetails[1]->GetMean();double sigma2 = histDetails[1]->GetMeanError();
    latex1.SetTextSize(0.025);
    latex1.DrawLatex(0.20, 0.65,Form("#splitline{#LT x_{J,1} #GT = %.3f #pm %.3f}{#LT x_{J,2} #GT = %.3f #pm %.3f}", mean1, sigma1, mean2, sigma2));
    std::cout << Form("Standard Deviation %.3f and %.3f", histDetails[0]->GetStdDev(), histDetails[1]->GetStdDev()) << std::endl << std::endl << std::endl;  
    

    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.png", fileName)); 
    c1->Close();
}

void CompGraphxJ3x5(std::vector<TH1D*> histDetails, const char* fileName , std::vector<const char*> newTitles ,const char* Llabel1,const char* Llabel2,std::vector<const char*> LatexLabels ,TDirectory* directory, bool print, bool logg ){
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 2400, 1200);
    
    c1->SetTicks(1,1);
    c1->Divide(6,3,0.001,0.001);
    c1->SetLeftMargin(0.15);
    
    int p = 0; int k = 0;
    int index = 0;
    if(logg) c1->SetLogy();
    for(int count = 1; count <= 18; count++ ){
        c1->cd(count);
        if((count) % 6== 0){
            TLatex latex1;
            latex1.SetNDC();
            latex1.SetTextSize(0.05);
            latex1.DrawLatex(0.06, 0.75,LatexLabels[p]);
            
            TLegend* leg = new TLegend(0.50,0.60,0.85,0.85);
            leg->SetBorderSize(0);
            leg->AddEntry(histDetails[0], Llabel1,"l");
            leg->AddEntry(histDetails[1], Llabel2,"l");
            leg->Draw();

            p++; k=0;
            continue;
        }
        
    
        histDetails[index]->SetTitle(newTitles[k]);
        histDetails[index]->GetXaxis()->SetNdivisions(510);
        histDetails[index]->GetYaxis()->SetNdivisions(510);
        histDetails[index]->GetYaxis()->SetRangeUser(0,histDetails[index+1]->GetMaximum()*1.2);

        histDetails[index]->SetLineColorAlpha(kRed+1,0.7); histDetails[index]->SetLineWidth(2);
        histDetails[index+1]->SetLineColorAlpha(kBlue+1,0.7); histDetails[index+1]->SetLineWidth(2);

        histDetails[index]->Draw("HIST E1");
        histDetails[index+1]->Draw("HIST E1 SAME");
        
        double mean1 = histDetails[index]->GetMean();double sigma1 = histDetails[index]->GetMeanError();
        double mean2 = histDetails[index+1]->GetMean();double sigma2 = histDetails[index+1]->GetMeanError();
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.03);
        latex1.DrawLatex(0.12, 0.80,Form("#splitline{#LT x_{J,1} #GT = %.3f #pm %.3f}{#LT x_{J,2} #GT = %.3f #pm %.3f}", mean1, sigma1, mean2, sigma2));
        std::cout << Form("Count: %d, Standard Deviation %.3f and %.3f",count,  histDetails[index]->GetStdDev(), histDetails[index+1]->GetStdDev()) << std::endl << std::endl << std::endl;  
        index+=2; k++;
    }

    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    c1->Close();
}

void Comp3x5Pyth_Data(std::vector<TH1D*> histDetails,std::vector<TH1D*> histPythiaDetails, const char* fileName , std::vector<const char*> newTitles ,const char* Llabel1,const char* Llabel2,std::vector<const char*> LatexLabels ,TDirectory* directory, bool print, bool vectorSum,bool logg ){
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 2400, 1200);
    
    
    c1->Divide(6,3,0.001,0.001);
    c1->SetLeftMargin(0.30);

    int p = 0; int k = 0;
    int index = 0;
    if(vectorSum) index = 1; 
    if(logg) c1->SetLogy();
    for(int count = 1; count <= 18; count++ ){
        c1->cd(count);
        gPad->SetTicks(1,1);
        if((count) % 6== 0){
            TLatex latex1;
            latex1.SetNDC();
            latex1.SetTextSize(0.05);
            latex1.DrawLatex(0.06, 0.75,LatexLabels[p]);
            latex1.DrawLatex(0.28, 0.25,"multijet x_{J} = #frac{p_{T,1}}{p_{T,2+3}}");
            
            TLegend* leg = new TLegend(0.45,0.55,0.85,0.85);
            leg->SetBorderSize(0);
            leg->AddEntry(histDetails[index-2], Llabel1,"l");
            leg->AddEntry(histPythiaDetails[index-2], Llabel2,"l");
            leg->Draw();

            p++; k = 0;
            continue;
        }
    
        histDetails[index]->GetXaxis()->SetNdivisions(510);histDetails[index]->GetXaxis()->SetRangeUser(0,1.5);
        histDetails[index]->GetYaxis()->SetNdivisions(510); 
        histDetails[index]->GetYaxis()->SetRangeUser(0,histDetails[index]->GetMaximum()*1.3);

        histDetails[index]->SetLineColorAlpha(kOrange+6,0.85); histDetails[index]->SetLineWidth(1);
        histPythiaDetails[index]->SetLineColorAlpha(kCyan+2,1); 

        histDetails[index]->Draw("HIST E1");
        histPythiaDetails[index]->Draw("HIST E1 SAME");
        
        double mean1 = histDetails[index]->GetMean();double sigma1 = histDetails[index]->GetMeanError();
        double mean2 = histPythiaDetails[index]->GetMean();double sigma2 = histPythiaDetails[index]->GetMeanError(); 
        double ratio = mean1/mean2; double ratioError = ratio * sqrt(pow(sigma1/mean1, 2) + pow(sigma2/mean2, 2));
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.04);
        latex1.DrawLatex(0.13, 0.82,newTitles[k]);
        latex1.DrawLatex(0.13, 0.72,Form("#splitline{#LT x_{J} #GT_{Data} = %.3f #pm %.3f}{#LT x_{J} #GT_{Pythia} = %.3f #pm %.3f}", mean1, sigma1, mean2, sigma2));
        std::cout << Form("Count: %d, Standard Deviation %.3f and %.3f",count,  histDetails[index]->GetStdDev(), histPythiaDetails[index]->GetStdDev()) << std::endl << std::endl << std::endl;  
        index+=2; k++;
    }

    directory->cd();
    c1->Write(fileName);
    histPythiaDetails[2]->Print("all");
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    c1->Close();
}

void Comp1x5Pyth_Data(std::vector<TH1D*> histDetails,std::vector<TH1D*> histPythiaDetails, const char* fileName , std::vector<const char*> newTitles ,const char* Llabel1,const char* Llabel2,std::vector<const char*> LatexLabels ,TDirectory* directory, bool print, bool vectorSum,bool logg ){
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 2700, 400);
    
    c1->SetTicks(1,1);
    c1->Divide(7,1,0.001,0.001);
    c1->SetLeftMargin(0.15);

    int k = 0;
    int index = 0;
    if(vectorSum) index = 21; 
    if(logg) c1->SetLogy();
    for(int count = 1; count <= 7; count++ ){
        c1->cd(count);
        gPad->SetTicks(1,1);
        if(count == 1){
            TLatex latex1;
            latex1.SetNDC();
            latex1.SetTextSize(0.05);
            if(vectorSum){
                latex1.DrawLatex(0.06, 0.75,LatexLabels[2]);
                latex1.DrawLatex(0.28, 0.25,"multijet x_{J} = #frac{p_{T,1}}{p_{T,2+3}}");
            }else{
            latex1.DrawLatex(0.06, 0.75,"#splitline{#color[2]{Combination:}}{"
            "#splitline{|#eta_{1,2,3}| < 0.7}{"
            "#splitline{#Delta#phi_{1,2} > 3#pi/4}{"
            "#splitline{p_{T,2} > 5 GeV}{"
            "#scale[0.75]{Scaled Trigger 22}}"
            "}}}");
            }
            TLegend* leg = new TLegend(0.45,0.55,0.85,0.85);
            leg->SetBorderSize(0);
            leg->AddEntry(histDetails[index-2], Llabel1,"l");
            leg->AddEntry(histPythiaDetails[index-2], Llabel2,"l");
            leg->Draw();

            k = 0;
            continue;
        }
    
        histDetails[index]->GetXaxis()->SetNdivisions(510);
        histDetails[index]->GetYaxis()->SetNdivisions(510);
        histDetails[index]->GetYaxis()->SetRangeUser(0,histDetails[index]->GetMaximum()*1.1);

        histDetails[index]->SetLineColorAlpha(kOrange+6,0.7); histDetails[index]->SetLineWidth(2);
        histPythiaDetails[index]->SetLineColorAlpha(kCyan+2,0.7); histPythiaDetails[index]->SetLineWidth(2);

        histDetails[index]->Draw("HIST E1");
        histPythiaDetails[index]->Draw("HIST E1 SAME");
        
        double mean1 = histDetails[index]->GetMean();double sigma1 = histDetails[index]->GetMeanError();
        double mean2 = histPythiaDetails[index]->GetMean();double sigma2 = histPythiaDetails[index]->GetMeanError(); 
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.03);
        latex1.DrawLatex(0.13, 0.82,newTitles[k]);
        latex1.DrawLatex(0.13, 0.72,Form("#splitline{#LT x_{J,1} #GT = %.3f #pm %.3f}{#LT x_{J,2} #GT = %.3f #pm %.3f}", mean1, sigma1, mean2, sigma2));
        std::cout << Form("Count: %d, Standard Deviation %.3f and %.3f",count,  histDetails[index]->GetStdDev(), histPythiaDetails[index]->GetStdDev()) << std::endl << std::endl << std::endl;  
        index+=2; k++;
    }

    directory->cd();
    c1->Write(fileName);
    histPythiaDetails[2]->Print("all");
    histDetails[2]->Print("all");
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    c1->Close();
}

void Comp1x7Pyth_Data(std::vector<TH1D*> histDetails,std::vector<TH1D*> histPythiaDetails, const char* fileName , std::vector<const char*> newTitles ,const char* Llabel1,const char* Llabel2,const char* LatexLabel ,TDirectory* directory, bool print ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 1600, 800);
    
    c1->SetTicks(1,1);
    c1->Divide(4,2,0.001,0.001);
    c1->SetLeftMargin(0.15);

    int index = 0;
    for(int count = 1; count <= 8; count++ ){
        c1->cd(count);
        gPad->SetTicks(1,1);
        if(count == 1){
            TLatex latex1;
            latex1.SetNDC();
            latex1.SetTextSize(0.05); 

            latex1.DrawLatex(0.06, 0.75,LatexLabel);
            latex1.DrawLatex(0.28, 0.25,"multijet x_{J} = #frac{p_{T,1}}{p_{T,2+3}}");
            
            histDetails[index]->SetLineColorAlpha(kOrange+6,0.7); histDetails[index]->SetLineWidth(1);
            histPythiaDetails[index]->SetLineColorAlpha(kCyan+2,0.9); histPythiaDetails[index]->SetLineWidth(1);

            TLegend* leg = new TLegend(0.50,0.55,0.90,0.95);
            leg->SetBorderSize(0);
            leg->AddEntry(histDetails[index], Llabel1,"l");
            leg->AddEntry(histPythiaDetails[index], Llabel2,"l");
            leg->Draw();
            continue;
        }
        
        if(count == 5){ //making it evenly stacked.
            continue;
        }
    
        histDetails[index]->GetXaxis()->SetNdivisions(510);
        histDetails[index]->GetYaxis()->SetNdivisions(510);
        histDetails[index]->GetYaxis()->SetRangeUser(0,histPythiaDetails[index]->GetMaximum()*1.4);
        

        histDetails[index]->SetLineColorAlpha(kOrange+6,0.7); 
        histPythiaDetails[index]->SetLineColorAlpha(kCyan+2,0.9); 

        histDetails[index]->Draw("HIST E1");
        histPythiaDetails[index]->Draw("HIST E1 SAME");
        
        
        double mean1 = histDetails[index]->GetMean();double sigma1 = histDetails[index]->GetMeanError();
        double mean2 = histPythiaDetails[index]->GetMean();double sigma2 = histPythiaDetails[index]->GetMeanError(); 
        double ratio = mean1/mean2; double ratioError = ratio * sqrt(pow(sigma1/mean1, 2) + pow(sigma2/mean2, 2));
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.05);
        latex1.DrawLatex(0.13, 0.82,newTitles[index]);
        latex1.DrawLatex(0.13, 0.72,Form("#splitline{#LT x_{J} #GT_{Data} = %.3f #pm %.3f}{#LT x_{J} #GT_{Pythia} = %.3f #pm %.3f}", mean1, sigma1, mean2, sigma2));
        latex1.DrawLatex(0.13, 0.62,Form("Ratio = %.3f #pm %.3f", ratio, ratioError));
        std::cout << Form("Count: %d, Standard Deviation %.3f and %.3f",count-1,  histDetails[index]->GetStdDev(), histPythiaDetails[index]->GetStdDev()) << std::endl;  
        std::cout << Form("the total number of events in %d: Data = %f, Pythia = %f , The hand calculated uncertainties are %.3f and %.3f",count-1, histDetails[index]->GetEntries(), 
        histPythiaDetails[index]->GetEntries(), histDetails[index]->GetStdDev() / sqrt(histDetails[index]->GetEntries()), histPythiaDetails[index]->GetStdDev() / sqrt(histPythiaDetails[index]->GetEffectiveEntries())) << std::endl << std::endl << std::endl;
        index++; 
    }

    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    c1->Close();
}

void CompGraphpT(std::vector<TH1D*> histDetails,const char* fileName , const char* newTitle ,const char* Llabel1,const char* Llabel2,const char* Llabel3, const char* Llabel4, const char* LatexLabel ,TDirectory* directory, bool print, bool logg ){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 900, 900);
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    
    if(logg) c1->SetLogy();

    
    TH1F *hr = c1->DrawFrame(0, 7, 70, FindlargestValue(histDetails)*1.3);
    hr->GetXaxis()->SetTitle("Calibrated p_{T} (GeV)");
    hr->GetYaxis()->SetTitle("Counts");
    hr->GetXaxis()->SetNdivisions(510);
    hr->GetYaxis()->SetNdivisions(510);
    hr->SetMinimum(0.5);
    hr->Draw();
    
    
    histDetails[0]->SetLineColorAlpha(kRed+1,0.7); histDetails[0]->SetLineWidth(2);
    histDetails[1]->SetLineColorAlpha(kBlue+1,0.7); histDetails[1]->SetLineWidth(2);
    histDetails[2]->SetLineColorAlpha(kGreen+4,1); histDetails[2]->SetLineWidth(2); histDetails[2]->SetLineStyle(7);
    histDetails[3]->SetLineColorAlpha(kViolet+5,0.7); histDetails[3]->SetLineWidth(2); histDetails[3]->SetLineStyle(9);
    
    
    histDetails[0]->Draw("HIST SAME");
    histDetails[3]->Draw("HIST SAME ");
    histDetails[1]->Draw("HIST SAME");
    histDetails[2]->Draw("HIST SAME ");
    

    TLegend* leg = new TLegend(0.65,0.45,0.85,0.60);
    leg->SetBorderSize(0);
    leg->SetMargin(0.4);
    leg->AddEntry(histDetails[0], Llabel1,"L");
    leg->AddEntry(histDetails[1], Llabel2,"L");
    leg->AddEntry(histDetails[2], Llabel3,"L");
    leg->AddEntry(histDetails[3], Llabel4,"L");
    leg->Draw();
    
    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.60, 0.85,LatexLabel);

    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    c1->Close();
}

void dPhiGraph(std::vector<TH1D*> histDetails,const char* fileName , const char* newTitle ,std::vector<const char*> Llabels, const char* LatexLabel, TDirectory* directory, bool print){
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Efficiency", 900, 900);
    c1->SetTicks(1,1);
    c1->SetLeftMargin(0.15);
    
    histDetails[0]->SetTitle(newTitle);
    histDetails[0]->GetXaxis()->SetNdivisions(510);
    histDetails[0]->GetYaxis()->SetNdivisions(510);
    histDetails[0]->GetYaxis()->SetRangeUser(0,FindlargestValue(histDetails)*1.2);

    
    histDetails[0]->SetLineColorAlpha(kRed+1,0.7);histDetails[0]->SetLineWidth(2);
    histDetails[1]->SetLineColorAlpha(kBlue+1,0.7); histDetails[1]->SetLineWidth(2);
    histDetails[2]->SetLineColorAlpha(kGreen+1,0.7); histDetails[2]->SetLineWidth(2);
    histDetails[3]->SetLineColorAlpha(kViolet+1,0.7); histDetails[3]->SetLineWidth(2);
    
    
    histDetails[0]->Draw("HIST");
    histDetails[1]->Draw("HIST SAME");
    histDetails[2]->Draw("HIST SAME");
    histDetails[3]->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.59,0.65,0.80,0.85);
    leg->SetTextSize(0.028); 
    leg->SetMargin(0.2); 
    leg->SetBorderSize(0);
    leg->AddEntry(histDetails[0], Llabels[0],"l");
    leg->AddEntry(histDetails[1], Llabels[1],"l");
    leg->AddEntry(histDetails[2], Llabels[2],"l");
    leg->AddEntry(histDetails[3], Llabels[3],"l");
    leg->Draw();

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.20, 0.77,LatexLabel);

    directory->cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    c1->Close();
}

void triggerEfficiency(std::vector<TH1D*> histOG, TH1D* PythiaHist, const char* fileName, TDirectory* directory, bool print){
    TCanvas *c1 = new TCanvas("c1", "Distribution", 750, 700);
    c1-> SetTicks(1,1);c1->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);

    histOG[0]->GetXaxis()->SetNdivisions(510);
    histOG[0]->GetYaxis()->SetNdivisions(510);
    histOG[0]->GetYaxis()->SetRangeUser(0,2);

    histOG[0]->SetMarkerStyle(20);
    histOG[0]->SetMarkerColor(kRed+1);
    histOG[0]->SetLineColor(kRed+1);

    histOG[1]->SetMarkerStyle(21);
    histOG[1]->SetMarkerColor(kAzure+2);
    histOG[1]->SetLineColor(kAzure+2);

    histOG[2]->SetMarkerStyle(22);
    histOG[2]->SetMarkerColor(kViolet+2);
    histOG[2]->SetLineColor(kViolet+2);

    PythiaHist->SetMarkerStyle(20);
    PythiaHist->SetMarkerColor(kBlack);
    PythiaHist->SetLineColor(kBlack);

    histOG[0]->Draw("E");
    histOG[1]->Draw("E SAME");
    histOG[2]->Draw("E SAME");
    PythiaHist->Draw("E SAME");

    TLegend* leg = new TLegend(0.68,0.75,0.84,0.85);
    leg->AddEntry(PythiaHist, "Pythia", "p");
    leg->AddEntry(histOG[0], "Scaled Trigger 18", "p");
    leg->AddEntry(histOG[1], "Scaled Trigger 22", "p");
    leg->AddEntry(histOG[2], "Scaled Trigger 34", "p");
    leg->Draw();
    
    TLine line(
        histOG[0]->GetXaxis()->GetXmin(), 1.0,
        histOG[0]->GetXaxis()->GetXmax(), 1.0
    );
    line.SetLineStyle(2);
    line.Draw();

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03);
    latex1.DrawLatex(0.20, 0.74,
    "#splitline{#color[2]{Run-24 #it{p+p}}}"
    "{#splitline{Dijet requirement: #Delta #phi > 3#pi/4}"
    "{p_{T,1} , p_{T}^{Reco} > 15 GeV}"
    "}");

    directory-> cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("%s.pdf", fileName)); 
    gStyle->SetOptStat(0);
    c1->Close();
}

void makeTGraph(TGraphErrors* graph, TGraphErrors* pythiaGraph, const char* fileName, const char* LatexLabel,  TDirectory* directory, bool print){
    TCanvas *c1 = new TCanvas("c1", "Distribution", 800, 800);
    c1-> SetTicks(1,1);c1->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    
    graph->GetYaxis()->SetTitle("#LT x_{J} #GT"); graph->GetXaxis()->SetTitle("p_{T,1} (GeV)");  
    graph->GetYaxis()->SetRangeUser(1,2.1);
    graph->SetTitle(" "); 
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlack);
    
    pythiaGraph->SetMarkerStyle(22);
    pythiaGraph->SetMarkerColor(kViolet + 7);
    pythiaGraph->SetLineColor(kViolet + 7);

    graph->Draw("AP");
    pythiaGraph->Draw("P SAME");

    TLegend* leg = new TLegend(0.62,0.10,0.84,0.25);
    leg->SetBorderSize(0);
    leg->AddEntry(graph, "Data", "p");
    leg->AddEntry(pythiaGraph, "RECO Pythia", "p");
    leg->Draw();

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03); 
    latex1.DrawLatex(0.21, 0.85,LatexLabel);
    
    directory-> cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    gStyle->SetOptStat(0);
    c1->Close();
}

void RatioTGraph(TGraphErrors* graph,  const char* fileName, const char* LatexLabel,  TDirectory* directory, bool print){
    TCanvas *c1 = new TCanvas("c1", "Distribution", 800, 800);
    c1-> SetTicks(1,1);c1->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    
    graph->GetYaxis()->SetRangeUser(0.85,1.1);  
    graph->SetTitle(" ");  
    graph->GetYaxis()->SetTitle("#LT x_{J} #GT_{Data}/#LT x_{J} #GT_{Pythia}"); graph->GetXaxis()->SetTitle("Calibrated p_{T,1} (GeV)");  
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlack);

    graph->Draw("AP");

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03); 
    latex1.DrawLatex(0.20, 0.85,LatexLabel);

    
    directory-> cd();
    graph->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    gStyle->SetOptStat(0);
    c1->Close();
}

void Ratio3TGraph(std::vector<TGraphErrors*> graph,  const char* fileName, const char* LatexLabel,  TDirectory* directory, bool print){
    TCanvas *c1 = new TCanvas("c1", "Distribution", 800, 800);
    c1-> SetTicks(1,1);c1->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    
    graph[0]->GetYaxis()->SetRangeUser(0.6,1.3);  
    graph[0]->SetTitle(" ");  
    graph[0]->GetYaxis()->SetTitle("#LT x_{J} #GT_{Data}/#LT x_{J} #GT_{Pythia}"); graph[0]->GetXaxis()->SetTitle("Calibrated p_{T,1} (GeV)");  
    
    graph[0]->SetMarkerStyle(20);
    graph[0]->SetMarkerColor(kGreen+2);
    graph[0]->SetLineColorAlpha(kGreen+2,0.7); 
    graph[0]->SetLineWidth(1);

    graph[1]->SetMarkerStyle(22);
    graph[1]->SetMarkerColor(kRed+2);
    graph[1]->SetLineColorAlpha(kRed+2,0.9); 
    graph[1]->SetLineWidth(1);

    graph[2]->SetMarkerStyle(23);
    graph[2]->SetMarkerColor(kViolet+6); 
    graph[2]->SetLineColorAlpha(kViolet+6,0.7); 
    graph[2]->SetLineWidth(1);

    graph[0]->Draw("AP");
    graph[1]->Draw("P SAME");
    graph[2]->Draw("P SAME");

    TLegend* leg = new TLegend(0.20,0.30,0.40,0.45);
    leg->SetBorderSize(0);
    leg->AddEntry(graph[0], "Full Data","p");
    leg->AddEntry(graph[1], "EMFRAC^{lead} #geq 0.7","p");
    leg->AddEntry(graph[2], "EMFRAC^{lead} < 0.7","p");
    leg->Draw();

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03); 
    latex1.DrawLatex(0.20, 0.85,LatexLabel);

        

        
    directory-> cd();
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    gStyle->SetOptStat(0);
    c1->Close();
}

void EMCALtGraph(TGraph* graph, TGraph* pythiaGraph, const char* fileName, const char* LatexLabel,  TDirectory* directory, bool print){
    TCanvas *c1 = new TCanvas("c1", "Distribution", 800, 800);
    c1-> SetTicks(1,1);c1->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    if (graph->GetN() == 0 && pythiaGraph->GetN() == 0) {
        std::cerr << "Both graphs are empty." << std::endl;
        return;
    }else{
        std::cout << Form("%d %d is the number of e",graph->GetN(),pythiaGraph->GetN()) << std::endl;   
    }
    
    graph->GetYaxis()->SetTitle("EMFRAC (units)"); graph->GetXaxis()->SetTitle("Calibrated p_{T}^{Lead}");  

    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlack);
    
    pythiaGraph->SetMarkerStyle(22);
    pythiaGraph->SetMarkerColor(kOrange + 7);
    pythiaGraph->SetLineColor(kOrange + 7);

    graph->Draw("AP");
    pythiaGraph->Draw("P SAME");

    TLegend* leg = new TLegend(0.67,0.40,0.84,0.55);
    leg->SetBorderSize(0);
    leg->AddEntry(graph, "Data", "p");
    leg->AddEntry(pythiaGraph, "RECO Pythia", "p");
    leg->Draw();

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.03); 
    latex1.DrawLatex(0.75, 0.85,LatexLabel);
    
    directory-> cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    gStyle->SetOptStat(0);
    c1->Close();
}

void makeTGraphCOMBINED(std::vector<TGraphErrors*> graph, std::vector<TGraphErrors*> pythiaGraph,  const char* fileName, TDirectory* directory, bool print){
    TCanvas *c1 = new TCanvas("c1", "Distribution", 800, 800);
    c1-> SetTicks(1,1);c1->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    
    TLegend* leg = new TLegend(0.63,0.60,0.84,0.85);
    leg->SetBorderSize(0);
    std::vector<int> markerColors = {kBlack, kBlue+1, kGreen+1, kViolet+1, kCyan+2, kOrange+7};
    std::vector<int> cuts = {3,5,7};
    graph[0]->SetTitle("#LT x_{J} #GT multijet vs Lead Jet p_{T}"); graph[0]->GetYaxis()->SetRangeUser(0.68,0.95); graph[0]->GetXaxis()->SetRangeUser(15,65);
    graph[0]->GetYaxis()->SetTitle("#LT x_{J} #GT"); graph[0]->GetXaxis()->SetTitle("p_{T,1} (GeV)");

    for(int i = 0 ; i < graph.size() ; i++){
        graph[i]->SetMarkerStyle(20+ i );
        graph[i]->SetMarkerColor(markerColors[i]);
        graph[i]->SetLineColor(markerColors[i]);

        pythiaGraph[i]->SetMarkerStyle(22+ i );
        pythiaGraph[i]->SetMarkerColor(markerColors[i+3]);
        pythiaGraph[i]->SetLineColor(markerColors[i+3]);
        
        if (i == 0 ) graph[i]->Draw("AP");
            else graph[i]->Draw("P SAME");
        pythiaGraph[i]->Draw("P SAME");

        leg->AddEntry(graph[i], Form("Data xJ, p_{T,3} > %d", cuts[i] ), "p");
        leg->AddEntry(pythiaGraph[i], Form("Pythia, p_{T,3} > %d", cuts[i] ), "p");
        leg->Draw();
    }

    directory-> cd();
    c1->Write(fileName);
    if (print) c1->SaveAs(Form("images/%s.pdf", fileName)); 
    gStyle->SetOptStat(0);
    c1->Close();
}