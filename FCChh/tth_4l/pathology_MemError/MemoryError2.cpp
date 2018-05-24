
#include <ROOT/TDataFrame.hxx>
#include "TLorentzVector.h"
#include <TSystem.h>
#include <ROOT/TVec.hxx>

// FCC event datamodel includes
//#include "datamodel/MCParticleCollection.h"
//#include "datamodel/MCParticleData.h"
#include "datamodel/ParticleData.h"
#include "datamodel/TaggedParticleData.h"

/*
 * Running from lxplus:
 *   
 *   # setup environment as in ../init.sh to get fcc libs
 *   # use Makefile to compile
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



   auto mergeElectronsAndMuons = [](std::vector<fcc::ParticleData> x, std::vector<fcc::ParticleData> y) {
     std::vector<fcc::ParticleData> result;
     result.reserve(x.size() + y.size());
     result.insert( result.end(), x.begin(), x.end() );
     result.insert( result.end(), y.begin(), y.end() );
     return result;

   };

   std::cout << "Apply simple selectors ..." << std::endl;
   auto selectors =  df
                      .Define("selected_electrons", selectLeptons, {"electrons", "electronITags"})
                      .Define("selected_muons", selectLeptons, {"muons", "muonITags"})
                      .Define("selected_leptons", mergeElectronsAndMuons, {"electrons", "muons"})

                    ;
  auto nentries = selectors.Count();
  std::cout << *nentries << std::endl;
   std::cout << "Writing snapshot to disk ..." << std::endl;
   selectors.Snapshot("events", "tree.root",
    { 
      
      "selected_leptons",
      "selected_electrons",
      "selected_muons",
      });

   return 0;
}
