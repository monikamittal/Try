// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Math/MathUtils.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {

  class CMS_FSQ_13_035_test : public Analysis {
  public:
    /// Constructor
    CMS_FSQ_13_035_test() : Analysis("CMS_FSQ_13_035_test") {    }
  
    /// Book histograms and initialise projections before the run
    void init() {

   //ofstream mout;
   // mout.open("rivet_out.txt",std::ios:app); 
    //Select all final state particles
    const FinalState fs(-5,5);
    addProjection(fs, "FS");

    vector<pair<PdgId,PdgId> > vidsZ;
      vidsZ.push_back(make_pair(ELECTRON, POSITRON));
      vidsZ.push_back(make_pair(MUON, ANTIMUON));

    //ZFinder with finalstate that include FSR // with true for clusterd photons
    FinalState fsZ(-5,5);
    ZFinder invfsZ(fsZ, -5., 5., 20*GeV, MUON, 60.0*GeV, 120.0*GeV, 0.1, true, false );
    addProjection(invfsZ, "INVFSZ"); 
     
    // FS modifier to exclude classes of particles from final state
    VetoedFinalState vfs(fsZ);
    vfs.addVetoOnThisFinalState(invfsZ); //Veto particles from a supplied final state
    addProjection(vfs, "VFS");
    addProjection(FastJets(vfs, FastJets::ANTIKT, 0.5), "Jets"); 

    /// @todo Initialise and register projections here
    //_h_Z_pT_normalised = bookHistogram1D(1, 1, 1);
    // _h_Mj1j2_normalised = bookHistogram1D(1, 1, 1); //(1,1,1) here gives reference to data histos d01-x01-y01
     _h_Mjj = bookHistogram1D( "_h_Mjj",5000,0,5000);
     _h_Nevents = bookHistogram1D( "_h_Nevents",20000,0,3000000);
     _h_750_JetMultiplicity = bookHistogram1D( "_h_750_JetMultiplicity",5,0,5);
     _h_750_Pt3rdJet =  bookHistogram1D("_h_750_Pt3rdJet",280,0.,280.);
     _h_750_SumPt = bookHistogram1D("_h_750_SumPt",250,0.,250.);
     _h_750_y3star = bookHistogram1D("_h_750_y3star",5,0,5);
     _h_1000_JetMultiplicity = bookHistogram1D( "_h_1000_JetMultiplicity",5,0,5);
     _h_1000_Pt3rdJet =  bookHistogram1D("_h_1000_Pt3rdJet",280,0.,280.);
     _h_1000_SumPt = bookHistogram1D("_h_1000_SumPt",250,0.,250.);
     _h_1000_y3star = bookHistogram1D("_h_1000_y3star",5,0,5);
     _h_1250_JetMultiplicity = bookHistogram1D( "_h_1250_JetMultiplicity",5,0,5);
     _h_1250_Pt3rdJet =  bookHistogram1D("_h_1250_Pt3rdJet",280,0,280);
     _h_1250_SumPt = bookHistogram1D("_h_1250_SumPt",250,0,250);
     _h_1250_y3star = bookHistogram1D("_h_1250_y3star",5,0,5);
     _h_200_JetMultiplicity = bookHistogram1D( "_h_200_JetMultiplicity",5,0,5);
     _h_200_Pt3rdJet =  bookHistogram1D("_h_200_Pt3rdJet",280,0,280);
     _h_200_SumPt = bookHistogram1D("_h_200_SumPt",250,0,250);
     _h_200_y3star = bookHistogram1D("_h_200_y3star",5,0,5);
     _h_Deta = bookHistogram1D("_h_Deta",20,0.,9.);
     _h_Dphi = bookHistogram1D("_h_Dphi",20,0,3.1415);
     _h_RPthard = bookHistogram1D("_h_RPthard",50,0.,1.);
     _h_zstar   =bookHistogram1D("_h_zstar",20,0,3);
    }

    bool ApplyElectronCutsForZee(double pt1, double pt2, double eta1, double eta2){
      bool isFid1 = ((fabs(eta1)<1.4442)||((fabs(eta1)>1.566)&&(fabs(eta1)<2.4)));
      bool isFid2 = ((fabs(eta2)<1.4442)||((fabs(eta2)>1.566)&&(fabs(eta2)<2.4)));
      if (isFid1 && isFid2 && pt1>20.0*GeV && pt2 >20.0*GeV) return true;
      else return false;
    }
    
    bool ApplyMuonCutsForZmm(double pt1, double pt2, double eta1, double eta2){
      bool isFid1 = ((fabs(eta1)<2.1));
      bool isFid2 = ((fabs(eta2)<2.1));
      if (isFid1 && isFid2 && pt1>20.0*GeV && pt2 >20.0*GeV) return true;
      else return false;
    } 
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      
      const ZFinder& invfsZ = applyProjection<ZFinder>(event, "INVFSZ");  
      bool isZmm(false);
      // bool isZee(false);
      bool isZ(false);
      
      isZ = (!(invfsZ.empty()));  
      const ParticleVector&  ZDecayProducts =  invfsZ.constituents();
      if (ZDecayProducts.size() < 2) vetoEvent;
      
      double pt1=-9999., pt2=-9999.;
      double phi1=-9999., phi2=-9999.;
      double eta1=-9999., eta2=-9999.;
      double dimuon_mass=-9999.;
      double Zmass = 91.2*GeV;
      
      if (isZ){
	pt1 = ZDecayProducts[0].momentum().pT();
	pt2 = ZDecayProducts[1].momentum().pT();
	eta1 = ZDecayProducts[0].momentum().eta();
	eta2 = ZDecayProducts[1].momentum().eta();
	phi1 = ZDecayProducts[0].momentum().phi();
	phi2 = ZDecayProducts[1].momentum().phi();
      }
      
      isZmm = isZ && ((fabs(ZDecayProducts[0].pdgId()) == 13) && (fabs(ZDecayProducts[1].pdgId()) == 13));
      // isZee = isZ && ((fabs(ZDecayProducts[0].pdgId()) == 11) && (fabs(ZDecayProducts[1].pdgId()) == 11));
      
      if(!((isZmm))){
	cout << "vetoEvent" << endl;
	vetoEvent;
      }  
      
      FourMomentum Z;
      Z = ZDecayProducts[0].momentum()+ZDecayProducts[1].momentum();	
      dimuon_mass= Z.mass();
      
      
      bool passBosonConditions = false;
      if(isZmm)passBosonConditions = ApplyMuonCutsForZmm(pt1,pt2,eta1,eta2);
      //    if(isZee)passBosonConditions = ApplyElectronCutsForZee(pt1,pt2,eta1,eta2);
      if(!passBosonConditions)vetoEvent;
      
      if( fabs(dimuon_mass - Zmass) > 15.0)vetoEvent; 
      
      //Obtain the jets.
      const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(15.0*GeV);
      Jets jet_list; 
      foreach (const Jet& j, jets) {
	const double jeta = j.momentum().eta();
	const double jphi = j.momentum().phi();
	const double jpt = j.momentum().pT();
	if((fabs(jeta) < 4.7) && (jpt >15.0*GeV)){
	      if ((deltaR(eta1, phi1, jeta, jphi) > 0.3) && (deltaR(eta2, phi2, jeta, jphi) > 0.3)){jet_list.push_back(j);}
//     continue;
	 }
       }
      //only events with 2jets i.e Njet==2 && Njet>=2
      if (jet_list.size() < 2 )vetoEvent;
      double jet1Pt = jet_list[0].momentum().pT();
      double jet2Pt = jet_list[1].momentum().pT();
      double jet1Eta = jet_list[0].momentum().eta();
      double jet2Eta = jet_list[1].momentum().eta();
      double jet1Phi = jet_list[0].momentum().phi();
      double jet2Phi = jet_list[1].momentum().phi();
      
      double eta_jmin = jet1Eta;
      double eta_jmax = jet2Eta;
      if( jet1Eta > jet2Eta ) {
	eta_jmin = jet2Eta;
	eta_jmax = jet1Eta;
      }
      
      
      
      if ( jet1Pt <50.0*GeV || jet2Pt <30.0*GeV || jet1Eta >4.7 || jet2Eta >4.7 ) vetoEvent;
      
      FourMomentum dijet;
      dijet = jet_list[0].momentum()+jet_list[1].momentum();
      double dijet_mass = dijet.mass();     //use to apply different cut on dijet mass
      double detaj1j2 = fabs(jet1Eta-jet2Eta);
      double dphij1j2 = fabs(jet1Phi-jet2Phi);
     _h_Deta->fill( detaj1j2,weight);
     _h_Dphi->fill(dphij1j2,weight);
     _h_Mjj->fill(dijet_mass,weight);

     double zstar=fabs(Z.rapidity()-(0.5*(jet_list[0].momentum().rapidity()+jet_list[1].momentum().rapidity())))/detaj1j2;
     double pthard =(Z + jet_list[0].momentum()+jet_list[0].momentum()).pT();
     double RPthard=pthard/(Z.pT()+jet1Pt+jet2Pt); 
     
      _h_RPthard->fill(RPthard,weight);
      _h_zstar->fill(zstar,weight);
      
      int njets_btw =-999;
      double pt_3rd =-999, SumPt =-999, y3star=-999;
      //If Dijet mass is greater than 750
      
      int iter_=0;
      if ( dijet_mass > 750.0*GeV) {
        if(iter_ ==0){ njets_btw = 0;}
        for( unsigned int i = 2; i < jet_list.size(); ++i) {
	  double PtJet_btw = jet_list[i].momentum().pT();
	  double EtaJet_btw =jet_list[i].momentum().eta();
          if(PtJet_btw > 15.0*GeV){
            if ((EtaJet_btw >eta_jmin) && (EtaJet_btw < eta_jmax) ){
	    iter_++;
	    if (iter_ ==1 ){ SumPt = 0.0;
	      pt_3rd = jet_list[i].momentum().pT();
	      y3star = fabs( Z.rapidity() - 0.5*(jet_list[0].momentum().rapidity() + jet_list[1].momentum().rapidity()));
	      
	    }
	    SumPt+=  jet_list[i].momentum().pT();
	    njets_btw++;
	  }
	}
       }
      }

//cout <<"EventWeight: " << weight <<"  njets :   "<< njets_btw << " pt_3rd :   "  <<pt_3rd  << "   y3star : "  << y3star << " SumPt:    "  << SumPt << std::endl; 
      
      _h_750_JetMultiplicity->fill(njets_btw,weight);
      _h_750_Pt3rdJet->fill(pt_3rd,weight);
      _h_750_SumPt->fill(SumPt,weight);
      _h_750_y3star->fill(y3star,weight);
      
      //If Dijet mass is greater than 1000
      iter_=0;
      njets_btw =-999, pt_3rd =-999, SumPt =-999, y3star=-999;
 
      if ( dijet_mass > 1000.0*GeV) {
        if (iter_ ==0){njets_btw =0;}
        for( unsigned int i = 2; i < jet_list.size(); ++i) {
          double PtJet_btw = jet_list[i].momentum().pT();
          double EtaJet_btw =jet_list[i].momentum().eta();
          if(PtJet_btw > 15.0*GeV){
          if ((EtaJet_btw >eta_jmin) && (EtaJet_btw < eta_jmax) ){
            iter_++;
            if (iter_ ==1 ){ SumPt = 0.0;
              pt_3rd = jet_list[i].momentum().pT();
              y3star =  fabs(Z.rapidity() - 0.5*(jet_list[0].momentum().rapidity() + jet_list[1].momentum().rapidity()));

            }
            SumPt+=  jet_list[i].momentum().pT();
            njets_btw++;
          }
        }
       }
      }

      _h_1000_JetMultiplicity->fill(njets_btw,weight);
      _h_1000_Pt3rdJet->fill(pt_3rd,weight);
      _h_1000_SumPt->fill(SumPt,weight);
      _h_1000_y3star->fill(y3star,weight);

          
  
      //If Dijet mass is greater than 200
       iter_=0;
       njets_btw =-999, pt_3rd =-999, SumPt =-999, y3star=-999;  
      if ( dijet_mass > 200.0*GeV) {
          if (iter_ ==0){njets_btw =0;}
        for( unsigned int i = 2; i < jet_list.size(); ++i) {
          double PtJet_btw = jet_list[i].momentum().pT();
          double EtaJet_btw =jet_list[i].momentum().eta();
          if(PtJet_btw > 15.0*GeV){
          if ((EtaJet_btw >eta_jmin) && (EtaJet_btw < eta_jmax) ){
            iter_++;
            if (iter_ ==1 ){ SumPt = 0.0;
              pt_3rd = jet_list[i].momentum().pT();
              y3star =  fabs(Z.rapidity() - 0.5*(jet_list[0].momentum().rapidity() + jet_list[1].momentum().rapidity()));

            }
            SumPt+=  jet_list[i].momentum().pT();
            njets_btw++;
          }
        }
       }
      }

      _h_200_JetMultiplicity->fill(njets_btw,weight);
      _h_200_Pt3rdJet->fill(pt_3rd,weight);
      _h_200_SumPt->fill(SumPt,weight);
      _h_200_y3star->fill(y3star,weight);

      
       //If Dijet mass is greater than 1250
       iter_=0;
        njets_btw =-999, pt_3rd =-999, SumPt =-999, y3star=-999;  
      if ( dijet_mass > 1250.0*GeV) {
           if (iter_ ==0){njets_btw =0;}
        for( unsigned int i = 2; i < jet_list.size(); ++i) {
          double PtJet_btw = jet_list[i].momentum().pT();
          double EtaJet_btw =jet_list[i].momentum().eta();
          if(PtJet_btw > 15.0*GeV){
          if ((EtaJet_btw >eta_jmin) && (EtaJet_btw < eta_jmax) ){
            iter_++;
            if (iter_ ==1 ){ SumPt = 0.0;
              pt_3rd = jet_list[i].momentum().pT();
              y3star =fabs(Z.rapidity() - 0.5*(jet_list[0].momentum().rapidity() + jet_list[1].momentum().rapidity()));

            }
            SumPt+=  jet_list[i].momentum().pT();
            njets_btw++;
          }
        }
       }
      }

      _h_1250_JetMultiplicity->fill(njets_btw,weight);
      _h_1250_Pt3rdJet->fill(pt_3rd,weight);
      _h_1250_SumPt->fill(SumPt,weight);
      _h_1250_y3star->fill(y3star,weight);

      
      /// @todo Do the event by event analysis here
      
    }
    
    
    /// Normalise histograms etc., after the run
    void finalize() {
            double crossSection = 0.888;

            _h_Nevents->fill(sumOfWeights());
	    /// @todo Normalise, scale and otherwise manipulate histograms here
             cout <<"here comes the crossection :  " <<crossSection << "   and the Sum of Weights :   " <<sumOfWeights() <<endl;
/*
	     scale(_h_750_JetMultiplicity, crossSection/sumOfWeights()); // norm to cross section
             scale(_h_750_Pt3rdJet, crossSection/sumOfWeights());
             scale(_h_750_SumPt, crossSection/sumOfWeights());
             scale(_h_750_y3star, crossSection/sumOfWeights()); 
	     scale(_h_1250_JetMultiplicity, crossSection/sumOfWeights()); // norm to cross section
	     scale(_h_1250_Pt3rdJet, crossSection/sumOfWeights());
             scale(_h_1250_SumPt, crossSection/sumOfWeights());
             scale(_h_1250_y3star, crossSection/sumOfWeights());
             scale(_h_1000_JetMultiplicity, crossSection/sumOfWeights()); // norm to cross section
             scale(_h_1000_Pt3rdJet, crossSection/sumOfWeights());
             scale(_h_1000_SumPt, crossSection/sumOfWeights());
             scale(_h_1000_y3star, crossSection/sumOfWeights());
             scale(_h_200_JetMultiplicity, crossSection/sumOfWeights()); // norm to cross section
             scale(_h_200_Pt3rdJet, crossSection/sumOfWeights());
             scale(_h_200_SumPt, crossSection/sumOfWeights());
             scale(_h_200_y3star, crossSection/sumOfWeights());
             scale(_h_Mjj, crossSection/sumOfWeights());       
  */     
	    // normalize(_h_YYYY); # normalize to unity


/*
          normalize(_h_750_JetMultiplicity);
          normalize(_h_750_Pt3rdJet);
          normalize(_h_750_SumPt);
          normalize(_h_750_y3star);     
          normalize(_h_1250_JetMultiplicity);
          normalize(_h_1250_Pt3rdJet);
          normalize(_h_1250_SumPt);
          normalize(_h_1250_y3star);
          normalize(_h_1000_JetMultiplicity);
          normalize(_h_1000_Pt3rdJet);
          normalize(_h_1000_SumPt);
          normalize(_h_1000_y3star);      
*/

    }
	private:
    // Data members like post-cuts event weight counters go here
	private:
    AIDA::IHistogram1D *_h_750_JetMultiplicity;
    AIDA::IHistogram1D *_h_750_Pt3rdJet;
    AIDA::IHistogram1D *_h_750_SumPt;
    AIDA::IHistogram1D *_h_750_y3star;
    AIDA::IHistogram1D *_h_1000_JetMultiplicity;
    AIDA::IHistogram1D *_h_1000_Pt3rdJet;
    AIDA::IHistogram1D *_h_1000_SumPt;
    AIDA::IHistogram1D *_h_1000_y3star;
    AIDA::IHistogram1D *_h_1250_JetMultiplicity;
    AIDA::IHistogram1D *_h_1250_Pt3rdJet;
    AIDA::IHistogram1D *_h_1250_SumPt;
    AIDA::IHistogram1D *_h_1250_y3star;
    AIDA::IHistogram1D *_h_200_JetMultiplicity;
    AIDA::IHistogram1D *_h_200_Pt3rdJet;
    AIDA::IHistogram1D *_h_200_SumPt;
    AIDA::IHistogram1D *_h_200_y3star;
    AIDA::IHistogram1D *_h_Mjj;
    AIDA::IHistogram1D *_h_Nevents;
    AIDA::IHistogram1D *_h_Dphi;
    AIDA::IHistogram1D *_h_Deta;
    AIDA::IHistogram1D *_h_RPthard;
    AIDA::IHistogram1D *_h_zstar;


};



// The hook for the plugin system
DECLARE_RIVET_PLUGIN(CMS_FSQ_13_035_test);

}
