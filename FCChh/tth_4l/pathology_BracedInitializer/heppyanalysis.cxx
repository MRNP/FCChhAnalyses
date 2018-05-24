
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


   auto getPts2 = [](std::vector<fcc::ParticleData> in){
     std::vector<float> result;
       for (int i = 0; i < in.size(); ++i) {
         result.push_back(sqrt(in[i].core.p4.px * in[i].core.p4.px + in[i].core.p4.py * in[i].core.p4.py));
       }
       return result;
   };




   std::cout << "Apply simple selectors ..." << std::endl;
   auto selectors =  df
                      .Define("electrons_pt", getPts2, {"electrons"}) /// here one column name in braces is fine

                    ;
  auto nentries = selectors.Count();
  std::cout << *nentries << std::endl;
   std::cout << "Writing snapshot to disk ..." << std::endl;
   selectors.Snapshot("events", "tree.root",
    { "electrons_pt" }); ///// <<<<<<<  this wo't wor

   return 0;
}
