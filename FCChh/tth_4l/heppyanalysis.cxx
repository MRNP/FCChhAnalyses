
#include <ROOT/TDataFrame.hxx>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TCanvas.h>
#include <TApplication.h>
#include <TSystem.h>
#include <ROOT/TVec.hxx>
//#include "datamodel/MCParticleCollection.h"
//#include "datamodel/ParticleData.h"

/*
 * Running from lxplus:
 *   source /cvmfs/sft.cern.ch/lcg/views/dev3/latest/x86_64-slc6-gcc62-opt/setup.sh
 *   g++ -g -o analysis heppyAnalysis.cxx `root-config --cflags --glibs` -lROOTVecOps -lROOTDataFrame
 */

   bool higgsresonancesort (TLorentzVector i ,TLorentzVector j) { return (abs( 125. -i.M())<abs(125.-j.M())); };
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
   //gSystem->Load("libdatamodel.so");
   //f1.MakeProject("dictsDelphes", "*", "RECREATE+");
   //gSystem->Load("dictsDelphes/dictsDelphes.so");

   using ints   = ROOT::Experimental::VecOps::TVec<int>;
   //using fourvectors   = ROOT::Experimental::VecOps::TVec<TLorentzVector>;
   using fourvectors   = std::vector<TLorentzVector>;
   using floats = ROOT::Experimental::VecOps::TVec<float>;
   using u_ints = ROOT::Experimental::VecOps::TVec<unsigned int>;
   //using particles = ROOT::Experimental::VecOps::TVec<fcc::ParticleData>;
   
   std::cout << "Creating TDataFrame ..." << std::endl;
   ROOT::Experimental::TDataFrame df("events", fname);
   //auto d_0_30 = df.Range(0, 10);
   //
   
   auto getInvMass = [] (floats& pxs, floats& pys, floats& pzs, floats& pes) {
     return sqrt(pes*pes - pxs*pxs - pys*pys - pzs * pzs);
     };

   auto getPts = [](floats &pxs, floats &pys){
       auto all_pts = sqrt(pxs * pxs + pys * pys);
       //std::cout << pxs[0] << "\t" << pxs.size() << std::endl;
       return all_pts;
   };

   auto getVecSize = [](ints& pxs) -> Int_t{
       return (Int_t)pxs.size();
   };

   auto getFourVectorWithCut = [](floats& pxs, floats& pys, floats& pzs, floats& pes, floats& iso) {
      fourvectors fv;
       for (int i = 0; i < pxs.size(); ++i) {
         if ((sqrt(pxs[i] * pxs[i] + pys[i] * pys[i]) > 20) && (iso[i] < 0.4)) {

          fv.emplace_back();
          fv.back().SetXYZM(pxs[i], pys[i], pzs[i], pes[i]);
         }
       }
       return fv;
   };
   auto getChargewithCut = [](ints& charge, floats& pxs, floats& pys, floats& iso) {
     ints chargewithcut;
     for (int i = 0; i < charge.size(); i++) {
       if ((sqrt(pxs[i]*pxs[i] + pys[i]*pys[i]) > 20.) && (iso[i] < 0.4)) {
         chargewithcut.emplace_back(charge[i]);
       }

     }

    return chargewithcut;
   };

   auto getPtFromFourVector = [](fourvectors & fv) ->floats {
     floats pts;
      for (int i = 0; i < fv.size(); ++i) {
        pts.emplace_back(fv[i].Pt());
      }
      return pts;

   };
   auto getMassFromFourVector = [](fourvectors & fv) ->floats {
     floats pts;
      for (int i = 0; i < fv.size(); ++i) {
        pts.emplace_back(fv[i].M());
      }
      return pts;

   };



   auto getPts2 = [](std::vector<float> &pxs, std::vector<float> &pys){
       floats all_pts;
       for (int i = 0; i < pxs.size(); ++i) {
         std::cout << pxs[i] << std::endl;
         all_pts.push_back(sqrt(pxs[i] * pxs[i] + pys[i] * pys[i]));
       }
       return all_pts;
   };

   auto id = [](floats &pxs){
       return pxs;
   };


   auto buildResonances = [](fourvectors& es, fourvectors & ms, ints& echarge, ints& mcharge) {
     fourvectors higgses;
     //// eeee
     //
     /* 
     if (es.size() >  3) {
       for (int i = 0; i < es.size() - 3; ++i) {
         for(int j=i+1; j < es.size() - 2; ++j) {
           for(int k=j+1; k < es.size() - 1; ++k) {
             for(int l=k+1; l < es.size(); ++l) {
               if(echarge[i] + echarge[j] + echarge[k] + echarge[l] == 0) {
                 higgses.push_back(TLorentzVector(es[i] + es[j] + es[k] + es[l]));
               }
             }
           }
         }
       }
     }
    //// mumumumu
     if (ms.size() >  3) {
       for (int i = 0; i < ms.size() - 3; ++i) {
         for(int j=i+1; j < ms.size() - 2; ++j) {
           for(int k=j+1; k < ms.size() - 1; ++k) {
             for(int l=k+1; l < ms.size(); ++l) {
               if(mcharge[i] + mcharge[j] + mcharge[k] + mcharge[l] == 0) {
                 higgses.push_back(TLorentzVector(ms[i] + ms[j] + ms[k] + ms[l]));
               }
             }
           }
         }
       }
     }
    //// eemumu
     if (ms.size() >  1 && es.size() > 1) {
       for (int i = 0; i < es.size() - 1; ++i) {
         for(int j=i+1; j < es.size(); ++j) {
           for(int k=0; k < ms.size() - 1; ++k) {
             for(int l=k+1; l < ms.size(); ++l) {
               if((echarge[i] + echarge[j] == 0) && (mcharge[k] + mcharge[l] == 0)) { 
               higgses.push_back(TLorentzVector(es[i] + es[j] + ms[k] + ms[l]));
               }
             }
           }
         }
       }
     }
     */
    //// mumumumu
     if (ms.size() >  3) {
       for (int i = 0; i < ms.size() - 3; ++i) {
         for(int j=i+1; j < ms.size() - 2; ++j) {
           for(int k=j+1; k < ms.size() - 1; ++k) {
             for(int l=k+1; l < ms.size(); ++l) {
               //if(mcharge[i] + mcharge[j] + mcharge[k] + mcharge[l] == 0) {
                 higgses.push_back(TLorentzVector(ms[i] + ms[j] + ms[k] + ms[l]));
               //}
             }
           }
         }
       }
     }
     std::sort(higgses.begin(), higgses.end(), higgsresonancesort);
     fourvectors finalhiggses;
     if (higgses.size() > 0) {
       finalhiggses.push_back(TLorentzVector(higgses[0]));
     }
     return finalhiggses;

   };

   std::cout << "Apply simple selectors ..." << std::endl;
   auto selectors = df//Define("m_pts", getPts, {"muons.core.p4.px", "muons.core.p4.py"})
                      //.Define("e_pts", getPts, {"electrons.core.p4.px", "electrons.core.p4.py"})
                      //.Define("e_pts_30", "e_pts[(e_pts > 20)]") // selector: sel_electrons
                      //.Define("m_pts_30", "m_pts[m_pts > 20]") // selector: sel_muons
                      //.Define("m_pts_sel1", "m_pts[m_pts > 20]") // selector: sel_muons
                      //.Define("e_size", getVecSize, {"m_pts_30"})
                      .Define("e_iso", id, {"electronITags.tag"})
                      .Define("m_iso", id, {"muonITags.tag"})
                      .Define("e_tlv", getFourVectorWithCut, {"electrons.core.p4.px", "electrons.core.p4.py","electrons.core.p4.pz","electrons.core.p4.mass", "e_iso"})
                      .Define("m_tlv", getFourVectorWithCut, {"muons.core.p4.px", "muons.core.p4.py","muons.core.p4.pz","muons.core.p4.mass", "m_iso"})
                      .Define("echarge", getChargewithCut, {"electrons.core.charge", "electrons.core.p4.px", "electrons.core.p4.py", "e_iso"})
                      .Define("mcharge", getChargewithCut, {"muons.core.charge", "muons.core.p4.px", "muons.core.p4.py", "m_iso"})
                      .Define("esize", getVecSize, {"echarge"})
                      .Define("msize", getVecSize, {"mcharge"})
                      .Define("higgses", buildResonances, {"e_tlv", "m_tlv", "echarge", "mcharge"})
                      .Define("e_pt_sel", getPtFromFourVector, {"e_tlv"})
                      .Define("mptsel", getPtFromFourVector, {"m_tlv"})
                      .Define("higgses_pt", getPtFromFourVector, {"higgses"})
                      .Define("higgses_m", getMassFromFourVector, {"higgses"})
                      .Define("e_m", getMassFromFourVector, {"m_tlv"})
   .Define("ex", id, {"muons.core.p4.px"}).Define("ey", id, {"muons.core.p4.py"}).Define("ez", id, {"muons.core.p4.pz"}).Define("em", id, {"muons.core.p4.mass"})

                    ;
  auto nentries = selectors.Count();
  std::cout << *nentries << std::endl;
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
   //selectors.Snapshot("events", "tree.root",{ "e_pts", "m_pts_30", "e_tlv"} );
   //selectors.Snapshot("events", "tree.root",{ "e_pt_sel","mptsel", "higgses_pt"});
   //selectors.Snapshot("events", "tree.root",{ "higgses_m","higgses_pt", "mptsel", "ex", "ey", "ez", "em"});
   selectors.Snapshot("events", "tree.root",{"m_iso","e_m", "higgses_m", "higgses_pt", "ex", "ey", "ez", "em"});

   return 0;
}
