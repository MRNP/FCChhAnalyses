
#include <ROOT/TDataFrame.hxx>
#include "TVector3.h"
#include <TCanvas.h>
#include <TApplication.h>
#include <TSystem.h>
#include <ROOT/TVec.hxx>

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
   f1.MakeProject("dictsDelphes", "*", "RECREATE+");
   gSystem->Load("dictsDelphes/dictsDelphes.so");
   //gSystem->Load("libdatamodel.so");

   using ints   = ROOT::Experimental::VecOps::TVec<int>;
   using floats = ROOT::Experimental::VecOps::TVec<float>;
   using u_ints = ROOT::Experimental::VecOps::TVec<unsigned int>;
   
   std::cout << "Creating TDataFrame ..." << std::endl;
   ROOT::Experimental::TDataFrame df("events", fname);
   //auto d_0_30 = df.Range(0, 10);

   auto getPts = [](floats &pxs, floats &pys){
       auto all_pts = sqrt(pxs * pxs + pys * pys);
       return all_pts;
   };

   std::cout << "Apply simple selectors ..." << std::endl;
   auto selectors = df.Define("m_pts", getPts, {"muons.core.p4.px", "muons.core.p4.py"})
                      .Define("e_pts", getPts, {"electrons.core.p4.px", "electrons.core.p4.py"})
                      .Define("e_pts_30", "e_pts[e_pts > 30]") // selector: sel_electrons
                      .Define("m_pts_30", "m_pts[m_pts > 30]") // selector: sel_muons
                    ;

   // TODO: Get Isolation from ITags and use in selectors
   // (I'm not sure if we need to deal with sumpt, since it is multiplied and then divided, have to look into that ...

   //auto get_sumpt = [](floats &pts){
   //    float sum = 0;
   //    for(auto pt : pts)
   //        sum += pt;
   //    return sum;
   //}//;

   // Define total pt for the particles
   // Heppy: ptc.iso.sumpt
   // iso-> Isolation: '''Holds the results of an isolation calculation.'''
   ////auto ptsums =  selectors.Define("m_sumpt", get_sumpt, {"m_pts"})
   //                        .Define("e_sumpt", get_sumpt, {"e_pts"})
   //                ;
   //
   //auto getVects = [](floats &pxs, floats &pys, floats &pzs){
   //   ROOT::Experimental::VecOps::TVec<TVector3> all_vects;
   //   for(int i=0; i<pxs.size(); ++i){
   //      all_vects.emplace_back(TVector3(pxs[i], pys[i], pzs[i]));
   //   }
   //   return all_vects;
   //};

   //// New df to quickly identify problems
   //auto vectors = df.Define("j_vects", getVects, {"pfjets04.core.p4.px", "pfjets04.core.p4.py", "pfjets04.core.p4.pz"})
   //                 .Define("e_vects", getVects, {"electrons.core.p4.px", "electrons.core.p4.py", "electrons.core.p4.pz"})
   //                 .Define("m_vects", getVects, {"muons.core.p4.px", "muons.core.p4.py", "muons.core.p4.pz"})
   //              ;

  // TODO: port LeptonicHiggsbuilder to c++

  // TODO: Event filtering based on jets and higgs candidates

   // TODO: figure out the problem here
   //using TVector3s ROOT::Experimental::VecOps::TVec<TVector3>; 
   //auto getDeltas = [](TVector3s ptcs_a, TVector3s ptcs_b){
   //   for(auto v : ptcs_a)
   //}

   // Not working (see example.root)
   // Different vector sizes
   //   m_sumpt : 100 elements
   //   m_pts : 86 elements
   //auto muons_h = ptsums.Define("good_muons", "m_sumpt[(m_sumpt / m_pts) < 0.1]").Histo1D("good_muons");

   //TCanvas c12;
   //muons_h->Draw();
   //c12.Draw();
   //
   //

   // TODO: Output to flat tree
   //
   std::cout << "Writing snapshot to disk ..." << std::endl;
   selectors.Snapshot("events", "tree.root", {"e_pts_30", "e_pts_30"} );
   return 0;
}
