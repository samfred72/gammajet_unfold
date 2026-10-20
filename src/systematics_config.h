// Systematic sources, shared by macros/compute_systematics.C (the combination) and
// drawing/draw_systematics.C (the pages). Header-only; not part of libgammajet_unfold.
#ifndef systematics_config_h
#define systematics_config_h

#include "ana.h"
#include <map>
#include <set>
#include <string>
#include <vector>
#include "TColor.h"

const int niterate = 2;
// herwig (generator) source switched off: the current Herwig production has wrong settings. To be
// restored once Herwig is reproduced - set true (no other change needed).
const bool includeHerwig = false;
// ana::systags minus nominal, plus the sources that are not reprocessings. A new systag needs
// systColors/systSources entries (.at() throws if missing).
const vector<string> systematics = [] {
  vector<string> v;
  for (const string & s : ana::systags) if (s != "nominal") v.push_back(s);
  if (includeHerwig) v.push_back("herwig");
  v.insert(v.end(), {"niterLow", "niterHigh", "priorSensitivity"});
  return v;
}();
const map<string,int> systColors = {
  {"JERhigh",      kRed+1},
  {"JERlow",       kOrange+7},
  {"emscale_high", kAzure+1},
  {"emscale_low",  kTeal+2},
  {"EMRhigh",      kGray+2},
  {"EMRlow",       kGray+1},
  {"jes_high",     kPink+6},
  {"jes_low",      kYellow+2},
  {"threejet",     kGreen+2},
  {"narrowBDT",    kBlue+1},
  {"narrowISO",    kMagenta+1},
  {"narrowBDTbkg", kBlue+3},
  {"narrowISObkg", kMagenta+3},
  {"wideISObkg",   kMagenta-4},
  {"timingwide",   kOrange+3},
  {"herwig",       kSpring+5},
  {"niterLow",     kCyan+2},
  {"niterHigh",    kViolet+1},
  {"priorSensitivity", kOrange+2},
};

// Per source: response generator (sim), Data reprocessing (dataSystag) and niter.
// priorSensitivity has no entry (special-cased below).
struct SystSource { string sim; string dataSystag; int niter; };
const map<string, SystSource> systSources = {
  {"JERhigh",      {"pythia", "JERhigh",      niterate}},
  {"JERlow",       {"pythia", "JERlow",       niterate}},
  {"emscale_high", {"pythia", "emscale_high", niterate}},
  {"emscale_low",  {"pythia", "emscale_low",  niterate}},
  {"EMRhigh",      {"pythia", "EMRhigh",      niterate}},
  {"EMRlow",       {"pythia", "EMRlow",       niterate}},
  {"jes_high",     {"pythia", "jes_high",     niterate}},
  {"jes_low",      {"pythia", "jes_low",      niterate}},
  {"threejet",     {"pythia", "threejet",     niterate}},
  {"narrowBDT",    {"pythia", "narrowBDT",    niterate}},
  {"narrowISO",    {"pythia", "narrowISO",    niterate}},
  {"narrowBDTbkg", {"pythia", "narrowBDTbkg", niterate}},
  {"narrowISObkg", {"pythia", "narrowISObkg", niterate}},
  {"wideISObkg",   {"pythia", "wideISObkg",   niterate}},
  {"timingwide",   {"pythia", "timingwide",   niterate}},
  {"herwig",       {"herwig", "nominal",      niterate}},
  {"niterLow",     {"pythia", "nominal",      niterate-1}},
  {"niterHigh",    {"pythia", "nominal",      niterate+1}},
};

// High/low sources split per bin by sign (from ana::asymmetricSystagPairs).
const set<string> asymmetricSystematics = [] {
  set<string> s;
  for (const auto & pr : ana::asymmetricSystagPairs) { s.insert(pr.first); s.insert(pr.second); }
  return s;
}();

// Colorblind-safe 8-hue palette, fixed order.
const int colorBlue      = TColor::GetColor("#2a78d6");
const int colorOrange    = TColor::GetColor("#eb6834");
const int colorAqua      = TColor::GetColor("#1baf7a");
const int colorUnfolding = TColor::GetColor("#eda100");
const int colorMagenta   = TColor::GetColor("#e87ba4");
const int colorPurity    = TColor::GetColor("#008300");
const int colorRed       = TColor::GetColor("#e34948");
const int colorGrey      = TColor::GetColor("#767676");
const int colorBrown     = TColor::GetColor("#8a5a2b"); // ninth source (timing); outside the 8-hue set

// Final-page grouping: a high/low pair shares one color and legend entry; a symmetric source is
// drawn with its negation. Purity and Unfolding members are drawn only as the combined curves.
struct DisplayGroup { string label; int color; vector<string> members; bool symmetric; };
const vector<DisplayGroup> displayGroups = [] {
  vector<DisplayGroup> v = {
  {"JER",       colorBlue,    {"JERhigh", "JERlow"},          false},
  {"JES",       colorOrange,  {"jes_high", "jes_low"},        false},
  {"EM scale",  colorAqua,    {"emscale_high", "emscale_low"},false},
  {"EMR",       colorGrey,    {"EMRhigh", "EMRlow"},          false},
  {"threejet",  colorMagenta, {"threejet"},                   true},
  {"timing",    colorBrown,   {"timingwide"},                 true},
  };
  if (includeHerwig) v.push_back({"herwig", colorRed, {"herwig"}, true});
  return v;
}();
// Combined into one symmetric Purity source.
const set<string> purityMembers = {
  "narrowBDT", "narrowISO", "narrowBDTbkg", "narrowISObkg", "wideISObkg"
};
// Combined into the asymmetric Unfolding source (see unfoldingUncUp/Down).
const set<string> unfoldingMembers = {
  "niterLow", "niterHigh", "priorSensitivity"
};

#endif
