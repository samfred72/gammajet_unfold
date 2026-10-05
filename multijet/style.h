#ifndef STYLE_H
#define STYLE_H

#include <vector>
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#pragma once  

void singleGraph(TH1D* histDetails,const char* fileName  , TDirectory* directory, bool print = false, bool logg = false, bool EMCAL = false);

void xJGraph(TH1D* histDetails,const char* fileName , const char* LLabel , TDirectory* directory, bool print = false,bool logg = false); 

void doubleGraph(TH1D* histDetails,TH1D* PythiaHIST,const char* fileName , const char* newTitle ,TDirectory* directory, bool print = false, bool logg = false);

void doubleGraphERatio(std::vector<TH1D*> histDetails,std::vector<TH1D*> PythiaHIST,const char* fileName , const char* newTitle ,TDirectory* directory, bool print=false, bool logg =false);

void doubleGraphTripple(std::vector<TH1D*> histDetails,std::vector<TH1D*>  PythiaHIST,const char* fileName , const char* LatexLabel ,std::vector<const char*> newTitle ,TDirectory* directory, bool print=false , bool logg =false);

void DanPlots(TH1D* histDetails,const char* fileName , const char* newTitle ,const char* LatexLabel , TDirectory* directory, bool print = false, bool logg = false);

void singleGraph2D(TH2D* histDetails,const char* fileName , const char* LatexLabel , TDirectory* directory, bool print = false, bool logg = false);

void trippletGraph2D(std::vector<TH2D*> histDetails,const char* fileName , const char* LatexLabel , std::vector<const char*> range, TDirectory* directory, bool print = false, bool logg = false);

void singleGraph3D(TH3D* histDetails,const char* fileName , TDirectory* directory, bool print = false, bool logg = false);

void CompGraphxJ(std::vector<TH1D*> histDetails,const char* fileName , const char* newTitle ,const char* Llabel1,const char* Llabel2,const char* LatexLabel , TDirectory* directory, bool print = false, bool logg = false);

void CompGraphxJ3x5(std::vector<TH1D*> histDetails,const char* fileName , std::vector<const char*> newTitle ,const char* Llabel1,const char* Llabel2,std::vector<const char*> LatexLabels  ,TDirectory* directory, bool print = false, bool logg = false);

void Comp3x5Pyth_Data(std::vector<TH1D*> histDetails,std::vector<TH1D*> histPythiaDetails, const char* fileName , std::vector<const char*> newTitles ,const char* Llabel1,const char* Llabel2,std::vector<const char*> LatexLabels ,TDirectory* directory, bool print = false, bool vectorSum = false,bool logg = false);

void Comp1x5Pyth_Data(std::vector<TH1D*> histDetails,std::vector<TH1D*> histPythiaDetails, const char* fileName , std::vector<const char*> newTitles ,const char* Llabel1,const char* Llabel2,std::vector<const char*> LatexLabels ,TDirectory* directory, bool print=false, bool vectorSum=false,bool logg =false);

void Comp1x7Pyth_Data(std::vector<TH1D*> histDetails,std::vector<TH1D*> histPythiaDetails, const char* fileName , std::vector<const char*> newTitles ,const char* Llabel1,const char* Llabel2,const char* LatexLabel  ,TDirectory* directory, bool print = false);

void CompGraphpT(std::vector<TH1D*> histDetails,const char* fileName , const char* newTitle ,const char* Llabel1,const char* Llabel2,const char* Llabel3, const char* Llabel4, const char* LatexLabel , TDirectory* directory, bool print = false, bool logg = false);

void dPhiGraph(std::vector<TH1D*> histDetails,const char* fileName , const char* newTitle ,std::vector<const char*> Llabels, const char* LatexLabel, TDirectory* directory, bool print = false);

void triggerEfficiency(std::vector<TH1D*> histOG, TH1D* PythiaHist, const char* fileName, TDirectory* directory, bool print = false);

void makeTGraph(TGraphErrors* graph, TGraphErrors* pythiaGraph, const char* fileName, const char* LatexLabel, TDirectory* directory, bool print = false);

void RatioTGraph(TGraphErrors* graph,  const char* fileName, const char* LatexLabel,  TDirectory* directory, bool print =false);

void Ratio3TGraph(std::vector<TGraphErrors*> graph,  const char* fileName, const char* LatexLabel,  TDirectory* directory, bool print= false);

void EMCALtGraph(TGraph* graph, TGraph* pythiaGraph, const char* fileName, const char* LatexLabel,  TDirectory* directory, bool print= false);

void makeTGraphCOMBINED(std::vector<TGraphErrors*> graph, std::vector<TGraphErrors*> pythiaGraph,  const char* fileName, TDirectory* directory, bool print = false); 

#endif
