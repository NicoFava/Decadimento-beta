#ifndef FUNZIONI_H
#define FUNZIONI_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <numeric>
#include <filesystem>
#include <algorithm>
#include <set>
#include <chrono>
#include <unordered_set>
#include <iomanip>

// Librerie di ROOT
#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TVector3.h"
#include "TLorentzVector.h"

const double M_NEUTRON = 939.565; // MeV
const double M_PROTON  = 938.272; // MeV
const double M_ELECTRON= 0.511;   // MeV

using namespace std;

struct Event {
    TLorentzVector p_proton;
    TLorentzVector p_electron;
    TLorentzVector p_neutrino;
};

// Dichiarazione delle funzioni
bool IsInPhysicalRegion(double Ee, double Ev);
void run_dalitz();
Event GenerateDecayEvent(TRandom3* gRandom, double Ee, double Ev);
TVector3 GenerateRandomDirection(TRandom3* gRandom);
void RotateVectorToPlane(TLorentzVector& vec, const TVector3& normal);

#endif