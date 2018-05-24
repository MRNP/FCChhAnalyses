
#include <ROOT/TDataFrame.hxx>
#include "TLorentzVector.h"
#include <TSystem.h>

// FCC event datamodel includes
//#include "datamodel/MCParticleCollection.h"
//#include "datamodel/MCParticleData.h"
#include "datamodel/ParticleData.h"
#include "datamodel/TaggedParticleData.h"

/*
 * Running from lxplus:
 *   source /cvmfs/sft.cern.ch/lcg/views/dev3/latest/x86_64-slc6-gcc62-opt/setup.sh
 *   g++ -g -o analysis heppyAnalysis.cxx `root-config --cflags --glibs` -lROOTVecOps -lROOTDataFrame
 */




// Test to reproduce Heppy analysis
int main(int argc, char* argv[]){


   // get input file from command line
   std::string fname;
   if (argc > 1) {
    fname = argv[1]; 
   } else { // no command line arguments given, use some default
     fname = "FCCDelphesOutput.root";
   }
   TFile f1(fname.c_str());

   std::cout << "Read file " << fname << std::endl;

   // fcc edm libraries
   gSystem->Load("libdatamodel.so");
   //f1.MakeProject("dictsDelphes", "*", "RECREATE+");
   //gSystem->Load("dictsDelphes/dictsDelphes.so");

   
   std::cout << "Creating TDataFrame ..." << std::endl;
   ROOT::Experimental::TDataFrame df("events", fname);

   // Range issue: 
   //auto d_0_30 = df.Range(0, 10);
   //

   auto selectLeptons = [](std::vector<fcc::ParticleData> in, std::vector<fcc::TaggedParticleData> iso) {
    std::vector<fcc::ParticleData> result;
    result.reserve(in.size());
    for (int i = 0; i < in.size(); ++i) {
      auto & p = in[i];
      if (std::sqrt(std::pow(p.core.p4.px,2) + std::pow(p.core.p4.py,2)) > 20) {
        if (iso[i].tag  < 0.4) {
          result.emplace_back(p);

        }
      }
    }
    return result;
   };



   auto getPts2 = [](std::vector<fcc::ParticleData> in){
     std::vector<float> result;
       for (int i = 0; i < in.size(); ++i) {
         result.push_back(sqrt(in[i].core.p4.px * in[i].core.p4.px + in[i].core.p4.py * in[i].core.p4.py));
       }
       return result;
   };

   auto id2 = [](std::vector<fcc::ParticleData> x) {
     return x;
   };

   auto mergeElectronsAndMuons = [](std::vector<fcc::ParticleData> x, std::vector<fcc::ParticleData> y) {
     std::vector<fcc::ParticleData> result;
     result.reserve(x.size() + y.size());
     result.insert( result.end(), x.begin(), x.end() );
     result.insert( result.end(), y.begin(), y.end() );
     return result;

   };

  auto LeptonicZBuilder = [](std::vector<fcc::ParticleData> leptons) {

        std::vector<fcc::ParticleData> result;
        int n = leptons.size();
        if (n >2) {
          std::vector<bool> v(n);
          std::fill(v.end() - 2, v.end(), true);
          do {
            fcc::ParticleData zed;
            zed.core.pdgId = 23;
            TLorentzVector zed_lv; 
            for (int i = 0; i < n; ++i) {
                if (v[i]) {
                  zed.core.charge += leptons[i].core.charge;
                  TLorentzVector lepton_lv;
                  lepton_lv.SetXYZM(leptons[i].core.p4.px, leptons[i].core.p4.py, leptons[i].core.p4.pz, leptons[i].core.p4.mass);
                  zed_lv += lepton_lv;
                }
            }
            zed.core.p4.px = zed_lv.Px();
            zed.core.p4.py = zed_lv.Py();
            zed.core.p4.pz = zed_lv.Pz();
            zed.core.p4.mass = zed_lv.M();
            result.emplace_back(zed);
          
          } while (std::next_permutation(v.begin(), v.end()));
        }

    return result;
  };

  auto LeptonicHiggsBuilder = [](std::vector<fcc::ParticleData> leptons) {

        std::vector<fcc::ParticleData> result;
        int n = leptons.size();
        if (n >2) {
          std::vector<bool> v(n);
          std::fill(v.end() - 2, v.end(), true);
          do {
            fcc::ParticleData zed;
            zed.core.pdgId = 25;
            TLorentzVector zed_lv; 
            for (int i = 0; i < n; ++i) {
                if (v[i]) {
                  zed.core.charge += leptons[i].core.charge;
                  TLorentzVector lepton_lv;
                  lepton_lv.SetXYZM(leptons[i].core.p4.px, leptons[i].core.p4.py, leptons[i].core.p4.pz, leptons[i].core.p4.mass);
                  zed_lv += lepton_lv;
                }
            }
            zed.core.p4.px = zed_lv.Px();
            zed.core.p4.py = zed_lv.Py();
            zed.core.p4.pz = zed_lv.Pz();
            zed.core.p4.mass = zed_lv.M();
            result.emplace_back(zed);

          
          } while (std::next_permutation(v.begin(), v.end()));
        }

    if (result.size() > 1) {
    auto  higgsresonancesort = [] (fcc::ParticleData i ,fcc::ParticleData j) { return (abs( 125. -i.core.p4.mass)<abs(125.-j.core.p4.mass)); };
    std::sort(result.begin(), result.end(), higgsresonancesort);

    std::vector<fcc::ParticleData>::const_iterator first = result.begin();
    std::vector<fcc::ParticleData>::const_iterator last = result.begin() + 1;
    std::vector<fcc::ParticleData> onlyBestHiggs(first, last);
    return onlyBestHiggs;
    } else {
    return result;
    }
  };




   std::cout << "Apply simple selectors ..." << std::endl;
   auto selectors =  df
                      .Define("selected_electrons", selectLeptons, {"electrons", "electronITags"})
                      .Define("selected_muons", selectLeptons, {"muons", "muonITags"})
                      .Define("selected_leptons", mergeElectronsAndMuons, {"selected_electrons", "selected_muons"})
                      .Define("zeds", LeptonicZBuilder, {"selected_leptons"})
                      .Define("selected_leptons_pt", getPts2, {"selected_leptons"})
                      .Define("zeds_pt", getPts2, {"zeds"})
                      .Define("higgs", LeptonicHiggsBuilder, {"zeds"})
                      .Define("higgs_pt", getPts2, {"higgs"})


                    ;
  auto nentries = selectors.Count();
  std::cout << "Count events: " <<  *nentries << std::endl;
  std::cout << "Writing snapshot to disk ..." << std::endl;
  selectors.Snapshot("events", "tree.root",
    { 
      
      "zeds",
      "zeds_pt",
      "selected_muons",
      "selected_leptons",
      "selected_electrons",
      "selected_leptons_pt",
      "higgs",
      "higgs_pt"
      }
    );

   return 0;
}
