//------------------------------------------------------------------------------
// Only Data 
// (Temporal solution to avoid the problems in the Trees w/o GEN branches!)
//------------------------------------------------------------------------------
//#ifdef __ISMC 
//#include "/gpfs/csic_users/brochero/PAF_7TeV/GenInfo_Data.h"
//#endif 

#include "TopLegacy.h"

#include <fstream>
#include <iostream>
#include <math.h>

#if !defined(__CINT__)
ClassImp(TopLegacy);
#endif

TopLegacy::TopLegacy(TTree* tree) : PAFAnalysis(tree)
{
}

//------------------------------------------------------------------------------
// InsideLoop
//------------------------------------------------------------------------------
void TopLegacy::InsideLoop(){

  // Weight
  if (sampleName.Contains("Data")) 
     weight = fWeight;
  else{
    weight = fWeight * fPUWeight->GetWeight((float)T_Event_nTruePU);//True      
    //weight = fWeight * fPUWeight->GetWeight((float)T_Event_nPU);//Observed      
  }

  (*nEvents)++;

  //----------------------------------------------------------------------------
  // Reset data members
  //----------------------------------------------------------------------------
  JetNumber.clear();
  BtagJetNumber.clear();

  SysJetVar  .clear();
  Pt_clean   .clear();
  Eta_clean  .clear();
  Btag_clean .clear();
  PartonFlavour_clean.clear();

  S_Muon_Charge.clear();
  S_Elec_Charge.clear();
  S_Elec_SC_Et .clear();

  S_Muon.clear();
  S_Elec.clear();

  //----------------------------------------------------------------------------
  // Acceptance (GEN info)
  //----------------------------------------------------------------------------
  Gen_Muon_Charge.clear();
  Gen_Elec_Charge.clear();
  Gen_Muon.clear();
  Gen_Elec.clear();
  NGen_Jet.clear();
  NGen_b.clear();
  PtGen_Jet.clear();
  PtGen_b.clear();

  MET_Gen         =0.0;
  nSGenMuon       =  0;
  nSGenElec       =  0;

  HT              =0.0;
  missingEt       =  0;
  dileptonInvMass =  0;
  nJets           =  0;
  nJetsBtag       =  0;
  nGenMuon        =  0;
  nGenElec        =  0;
  nGenTau         =  0;
  nSelMuon        =  0;
  nSelElec        =  0;

  var             =0.0;

  //----------------------------------------------------------------------------
  // Acceptance: Only fot ttbar Samples
  //----------------------------------------------------------------------------

  if ((sampleName.Contains("TTbar") || sampleName.Contains("TTJets")) && !fileSuffix.Contains("Bkg")){
    
    MET_Gen=T_METgen_ET;
    
    //----------------------------------------------------------------------------
    // Select muons and electrons at Generation Level (Includes Taus)
    //----------------------------------------------------------------------------
    GetGenMuon();
    GetGenElec();
        
    //----------------------------------------------------------------------------
    // Jets and btag at GEN level (particle or parton?)
    //----------------------------------------------------------------------------
    for(unsigned int jjet=0; jjet<T_JetAKCHS_GenJet_Px->size(); jjet++){
      
      TLorentzVector GenJet(T_JetAKCHS_GenJet_Px    ->at(jjet),
			    T_JetAKCHS_GenJet_Py    ->at(jjet),
			    T_JetAKCHS_GenJet_Pz    ->at(jjet),
			    T_JetAKCHS_GenJet_Energy->at(jjet));
      
      if(GenJet.Pt()>30 && fabs(GenJet.Eta())<2.4){
	
	NGen_Jet.push_back(jjet);
	PtGen_Jet.push_back(GenJet.Pt());

	// b-jets
 	int GenFlavor=0;
	GenFlavor=T_JetAKCHS_Parton_Flavour->at(jjet);
	
	if(fabs(GenFlavor) == 5){ // ID for b quarks +/- 5
	  NGen_b.push_back(jjet);	
	  PtGen_b.push_back(GenJet.Pt());
	} // if(bjet)
	
      }// if(GenJet eta && pT)
    }// for(jjet)
    
    //----------------------------------------------------------------------------
    // Organize all jets by pT
    //----------------------------------------------------------------------------
    for (unsigned int ii=0;ii<NGen_Jet.size();ii++){
      for (unsigned int jj=ii;jj<NGen_Jet.size();jj++){
	
	if(PtGen_Jet[ii]<PtGen_Jet[jj]){
	  
	  float tempN;
	  tempN=NGen_Jet[ii];
	  NGen_Jet[ii]=NGen_Jet[jj];
	  NGen_Jet[jj]=tempN;
	  float tempPt;
	  tempPt=PtGen_Jet[ii];
	  PtGen_Jet[ii]=PtGen_Jet[jj];
	  PtGen_Jet[jj]=tempPt;
	}
      }
    }
    
    
    //----------------------------------------------------------------------------
    // Organize by pT: b-Jets
    //----------------------------------------------------------------------------    
    for (unsigned int ii=0;ii<NGen_b.size();ii++){
      for (unsigned int jj=ii;jj<NGen_b.size();jj++){
	
	if(PtGen_b[ii]<PtGen_b[jj]){
	  
	  float tempN;
	  tempN=NGen_b[ii];
	  NGen_b[ii]=NGen_b[jj];
	  NGen_b[jj]=tempN;
	  float tempPt;
	  tempPt=PtGen_b[ii];
	  PtGen_b[ii]=PtGen_b[jj];
	  PtGen_b[jj]=tempPt;
 	}
      }
    }

    
    
    TLorentzVector GenMuon1;
    TLorentzVector GenMuon2;
    
    TLorentzVector GenElec1;
    TLorentzVector GenElec2;
    
    //----------------------------------------------------------------------------
    // Filling of the TopTree with GEN info
    //----------------------------------------------------------------------------    
    if(nSGenMuon==2 && 
       (Gen_Muon_Charge[0]*Gen_Muon_Charge[1]<0.)) FillTree(Gen_Muon[0], Gen_Muon[1], 3);// mumu     
    if(nSGenElec==2 && 
       (Gen_Elec_Charge[0]*Gen_Elec_Charge[1]<0.)) FillTree(Gen_Elec[0], Gen_Elec[1], 4);// ee 
    if(nSGenMuon==1 && 
       nSGenElec==1 && 
       (Gen_Muon_Charge[0]*Gen_Elec_Charge[0]<0.)) FillTree(Gen_Elec[0], Gen_Muon[0], 5);// mue || emu 
    
  } // if (sampleName.Contains("TTbar")) 

  
  //----------------------------------------------------------------------------
  // MET Type I Corrected 13-004. Also change the variable in the FillTree
  //----------------------------------------------------------------------------
  // MC
  // missingEt = T_METType1CorrectedPFMetSmeared_ET;

  // metx      = T_METType1CorrectedPFMetSmeared_ET*TMath::Cos(T_METType1CorrectedPFMetSmeared_Phi);
  // mety      = T_METType1CorrectedPFMetSmeared_ET*TMath::Sin(T_METType1CorrectedPFMetSmeared_Phi);
  // metxSys   = metx;
  // metySys   = mety;

  // DATA
  missingEt = T_METPFTypeI_ET;

  metx      = T_METPFTypeI_ET*TMath::Cos(T_METPFTypeI_Phi);
  mety      = T_METPFTypeI_ET*TMath::Sin(T_METPFTypeI_Phi);
  metxSys   = metx;
  metySys   = mety;

  //----------------------------------------------------------------------------
  // Accept only events with a good vertex
  //----------------------------------------------------------------------------
  if (SelectedVertexIndex() < 0) return;

  //----------------------------------------------------------------------------
  // Select muons and electrons
  //----------------------------------------------------------------------------
  GetSelectedMuon();
  GetSelectedElec();

  //----------------------------------------------------------------------------
  // Get number of generated leptons 
  //----------------------------------------------------------------------------
  if(!sampleName.Contains("Data")) SelectedGenLepton();

  //----------------------------------------------------------------------------
  // Select the pair of leptons
  //----------------------------------------------------------------------------
  TLorentzVector Muon1;
  TLorentzVector Muon2;

  TLorentzVector Elec1;
  TLorentzVector Elec2;

  Bool_t mumuEvent = false;
  Bool_t eeEvent   = false;
  Bool_t emuEvent  = false;

  Int_t flav_seleclep1 = -1;
  Int_t flav_seleclep2 = -1;
  Int_t ind_seleclep1  = -1;
  Int_t ind_seleclep2  = -1;
  
  Double_t Wpair_max = -1;  // Max pt1 + pt2

  //----------------------------------------------------------------------------
  // Loop over the 1st lepton
  //----------------------------------------------------------------------------
  for (UInt_t i=0; i<nSelElec+nSelMuon; i++) {

    Double_t charge_lep1, charge_lep2, pt_lep1, pt_lep2;
    Int_t    flav_lep1, flav_lep2;
    
    if (i < nSelElec) {

      charge_lep1 = S_Elec_Charge[i];
      pt_lep1     = S_Elec[i].Pt();
      flav_lep1   = electronFlavor;

    } else {

      charge_lep1 = S_Muon_Charge[i-nSelElec];
      pt_lep1     = S_Muon[i-nSelElec].Pt();
      flav_lep1   = muonFlavor;    
    }

    //--------------------------------------------------------------------------
    // Loop over the 2nd lepton
    //--------------------------------------------------------------------------
    for (UInt_t j=i+1; j<nSelElec+nSelMuon; j++) {
      
      if  (j < nSelElec)  { 
	
	charge_lep2 = S_Elec_Charge[j];
	pt_lep2     = S_Elec[j].Pt();
	flav_lep2   = electronFlavor;
	
      } else { 

	charge_lep2 = S_Muon_Charge[j-nSelElec];
	pt_lep2     = S_Muon[j-nSelElec].Pt();
	flav_lep2   = muonFlavor;    	
      }

      //------------------------------------------------------------------------
      // Check for opposite charge leptons
      //------------------------------------------------------------------------
      if (charge_lep1*charge_lep2 < 0) {
	
	Double_t W_pair = pt_lep1 + pt_lep2;
	
	if (W_pair> Wpair_max) {
	  
	  Wpair_max      = W_pair;
	  ind_seleclep1  = i;
	  ind_seleclep2  = j;
	  flav_seleclep1 = flav_lep1;
	  flav_seleclep2 = flav_lep2;
	}
      }// Opposite charge lepton 
    }// loop over the 2nd lepton
  }// Loop over leptons

  //----------------------------------------------------------------------------
  // Select the lepton pair with highest weight
  //----------------------------------------------------------------------------
  if (Wpair_max > 0) {
    
    Int_t indfinal_seleclep1 = -1;
    Int_t indfinal_seleclep2 = -1;
    
    if (flav_seleclep1 == 1) indfinal_seleclep1 = ind_seleclep1;
    if (flav_seleclep1 == 2) indfinal_seleclep1 = ind_seleclep1 - nSelElec;
    
    if (flav_seleclep2 == 1) indfinal_seleclep2 = ind_seleclep2;
    if (flav_seleclep2 == 2) indfinal_seleclep2 = ind_seleclep2 - nSelElec;
    
    //--------------------------------------------------------------------------
    // Final classification
    //--------------------------------------------------------------------------
    if (flav_seleclep1 == 2 && flav_seleclep2 == 2) {

      mumuEvent = true;

      Muon1 = S_Muon[indfinal_seleclep1];
      Muon2 = S_Muon[indfinal_seleclep2];

    } else if (flav_seleclep1 == 1 && flav_seleclep2 == 1) {

      eeEvent = true;

      Elec1 = S_Elec[indfinal_seleclep1];
      Elec2 = S_Elec[indfinal_seleclep2];
      
    } else if (flav_seleclep1 == 1 && flav_seleclep2 == 2) {

      emuEvent = true;
      
      Elec1 = S_Elec[indfinal_seleclep1];
      Muon1 = S_Muon[indfinal_seleclep2];

    } else if (flav_seleclep1 == 2 && flav_seleclep2 == 1) {

      emuEvent = true;
      
      Muon1 = S_Muon[indfinal_seleclep1];
      Elec1 = S_Elec[indfinal_seleclep2];

    }
  }

  //----------------------------------------------------------------------------  
  // For the jet cleaning
  //----------------------------------------------------------------------------
  TVector3 p_particle_jc[2];

  if (mumuEvent){

    p_particle_jc[0].SetXYZ(Muon1.Px(), Muon1.Py(), Muon1.Pz());
    p_particle_jc[1].SetXYZ(Muon2.Px(), Muon2.Py(), Muon2.Pz());
    
  } 
  else if(eeEvent){
    
    p_particle_jc[0].SetXYZ(Elec1.Px(), Elec1.Py(), Elec1.Pz());
    p_particle_jc[1].SetXYZ(Elec2.Px(), Elec2.Py(), Elec2.Pz());
    
  } 
  else if (emuEvent){
    
    p_particle_jc[0].SetXYZ(Elec1.Px(), Elec1.Py(), Elec1.Pz());
    p_particle_jc[1].SetXYZ(Muon1.Px(), Muon1.Py(), Muon1.Pz());
  }

  //------------------------------------------------------------------------
  // Jet Selection
  //------------------------------------------------------------------------
  TVector3 vJet;
  TVector3 vGenJet;

  UInt_t   jSize   =    0;

  jSize = T_JetAKCHS_Px->size(); 

  if (mumuEvent || eeEvent || emuEvent) {// if(dilepton)

    for (UInt_t jt=0; jt<jSize; jt++) {

      vJet.SetXYZ(T_JetAKCHS_Px ->at(jt), T_JetAKCHS_Py ->at(jt), T_JetAKCHS_Pz ->at(jt));

      //------------------------------------------------------------------------
      // Systematic Variations: JES and JER
      //------------------------------------------------------------------------
      var = 0; 
      Double_t PtVar = 0; 

      if ( sys_source == sys_jes || sys_source == sys_jer ){

	Double_t jecCorr = 1.0; 
	// Only defined for MC (JEC) 
	// Also uncomment in the FillTree
	//Double_t jecCorr = T_JetAKCHS_Corr->at(jt); 
	
      //------------------------------------------------------------------------
      // JES
      //------------------------------------------------------------------------
      if ( sys_source == sys_jes ) {
	
	Double_t jecUnc  = 1.;
	// Only defined for MC (JEC) 
	// Also uncomment in the FillTree
	//Double_t jecUnc  = T_JetAKCHS_Uncertainty->at(jt); 
	
	Double_t unc = sqrt(jecUnc*jecUnc);
	if ( sys_direction == sys_Up)  var  = 1. + unc;
	else if ( sys_direction == sys_Down ) var = 1. - unc;

	PtVar   = var*vJet.Pt();
	
      }

      //------------------------------------------------------------------------
      // JER
      //------------------------------------------------------------------------
      if ( sys_source == sys_jer ) {
	
	Double_t factor_up=0, factor_nom=0, factor_down=0;
	
	//---------------------------------------------------------------------------------------------
	// Twiky Recipe (AN CMS-2014/416)
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
	//---------------------------------------------------------------------------------------------
	if (luminosity > 10000){// 2012 8TeV Data
	  if (fabs (vJet.Eta()) < 0.5 ) {    
	    factor_up   =1.105;
	    factor_nom  =1.079;
	    factor_down =1.053;
	  }
	  else if (fabs (vJet.Eta()) < 1.1 ) {    
	    factor_up   =1.127;
	    factor_nom  =1.099;
	    factor_down =1.071;
	  }
	  else if (fabs (vJet.Eta()) < 1.7 ) {    
	    factor_up   =1.150;
	    factor_nom  =1.121;
	    factor_down =1.092;
	  }
	  else if (fabs (vJet.Eta()) < 2.3 ) {    
	    factor_up   =1.254;
	    factor_nom  =1.208;
	    factor_down =1.162;
	  }
	  else if (fabs (vJet.Eta()) < 2.8 ) {    
	    factor_up   =1.316;
	    factor_nom  =1.254;
	    factor_down =1.192;
	  }
	  else if (fabs (vJet.Eta()) < 3.2 ) {    
	    factor_up   =1.458;
	    factor_nom  =1.395;
	    factor_down =1.332;
	  }
	  else if (fabs (vJet.Eta()) < 5.0 ) {    
	    factor_up   =1.247;
	    factor_nom  =1.056;
	    factor_down =0.865;
	  }
	}// 2012 8TeV Data

	else {// 2011 7TeV Data
	  if (fabs (vJet.Eta()) < 0.5 ) {    
	    factor_up   =1.115;
	    factor_nom  =1.052;
	    factor_down =0.990;
	  }
	  else if (fabs (vJet.Eta()) < 1.1 ) {    
	    factor_up   =1.114;
	    factor_nom  =1.057;
	    factor_down =1.001;
	  }	
	  else if (fabs (vJet.Eta()) < 1.7 ) {    
	    factor_up   =1.161;
	    factor_nom  =1.096;
	    factor_down =1.032;
	  }	
	  else if (fabs (vJet.Eta()) < 2.3 ) {    
	    factor_up   =1.228;
	    factor_nom  =1.134;
	    factor_down =1.042;
	  }	
	  else if (fabs (vJet.Eta()) < 5.0 ) {    
	    factor_up   =1.488;
	    factor_nom  =1.288;
	    factor_down =1.089;
	  }	
	}// 2011 7TeV Data
	
	// Jet Generation Level (JER)
	vGenJet.SetXYZ(T_JetAKCHS_GenJet_Px ->at(jt), T_JetAKCHS_GenJet_Py ->at(jt), T_JetAKCHS_GenJet_Pz ->at(jt));
	if ( vGenJet.Pt() < 15 ) continue; // Compare JER using Gen-Jets with pT>15GeV

	Double_t RecoGenMatch_DeltaR = vJet.DeltaR(vGenJet);// Matching GEN-RECO jets
	if(RecoGenMatch_DeltaR > 0.5) continue;

	Double_t deltaPt = (vJet.Pt() - vGenJet.Pt()); 
	if      ( sys_direction == sys_Up)   var = max(0.0, (vGenJet.Pt() + deltaPt*factor_up)/vJet.Pt());
	else if ( sys_direction == sys_Nom)  var = max(0.0, (vGenJet.Pt() + deltaPt*factor_nom)/vJet.Pt());
	else if ( sys_direction == sys_Down) var = max(0.0, (vGenJet.Pt() + deltaPt*factor_down)/vJet.Pt());
	
	PtVar   = var*vJet.Pt();
	
      } // if(JER)
      
      //------------------------------------------------------------------------
      // JES + MET OR JER + MET: Propagation to the MET
      //------------------------------------------------------------------------
      Double_t jetuncorrectedx = T_JetAKCHS_Px->at(jt)/jecCorr;
      Double_t jetuncorrectedy = T_JetAKCHS_Py->at(jt)/jecCorr;
      
      // Only propagate to jets outside of a 0.5 cone from the leptons:
      // jet cleanning, but not cut in pT, eta or ID
      TVector3 vJetpf;
      vJetpf.SetXYZ(T_JetAKCHS_Px->at(jt),T_JetAKCHS_Py->at(jt),T_JetAKCHS_Pz->at(jt));  

      Double_t dRminformet = 1.e9;
      for (UInt_t kkk=0; kkk<2; kkk++){
	Double_t dRformet = vJetpf.DeltaR(p_particle_jc[kkk]);
	if (dRformet < dRminformet) dRminformet = dRformet;
      }
      
      if ( dRminformet > 0.5 ) { 
	metxSys += jetuncorrectedx;
	metySys += jetuncorrectedy;
	
	metxSys -= jetuncorrectedx*var;
	metySys -= jetuncorrectedy*var; 
      }
      
      }//if(JES || JER)   
      
      Double_t JetPt = vJet.Pt(); // Nominal
      if ( sys_source == sys_jes || sys_source == sys_jer )  JetPt = PtVar; // Systematic variation
      
      //------------------------------------------------------------------------
      // Jet ID
      //------------------------------------------------------------------------
      Bool_t jetID = false;
      // do not include T_JetAKCHS_nDaughters->at(jt), was set to 0.
      jetID =  (T_JetAKCHS_NeutHadEnergyFrac->at(jt) < 0.99 && 
                T_JetAKCHS_NeutEmEnergyFrac->at(jt) < 0.99 &&
                (fabs(T_JetAKCHS_Eta->at(jt)) < 2.4 ? (T_JetAKCHS_CharEmEnergyFrac->at(jt) < 0.99 && 
						       T_JetAKCHS_CharHadEnergyFrac->at(jt) > 0. && T_JetAKCHS_ChargedMultiplicity->at(jt) > 0.) : true)); 
      
      if (JetPt > 30 && fabs(vJet.Eta()) < 2.4 && jetID) {
	
	//------------------------------------------------------------------------
	// Select jets outside of a 0.5 cone from the SELECTED leptons
	//------------------------------------------------------------------------
	Double_t dRmin = 999;
	
	for (UInt_t kk=0; kk<2; kk++) {	  
	  Double_t dR = vJet.DeltaR(p_particle_jc[kk]);
	  if (dR < dRmin) dRmin = dR;
	}
	if (dRmin < 0.5) continue; 
	
	//----------------------------------------------------------------------
	// Fill the jets
	//----------------------------------------------------------------------
	
	nJets++;
	
	JetNumber.push_back(jt);
	Eta_clean.push_back (vJet.Eta());
	Pt_clean.push_back (vJet.Pt()); //Systematic Variation in the jet pT must be applied in the TreeReader.C
	SysJetVar.push_back(var);	

	HT+= T_JetAKCHS_Et ->at(jt);

	Btag_clean.push_back(T_JetAKCHS_Tag_CombSVtx->at(jt)); 
	PartonFlavour_clean.push_back(T_JetAKCHS_Parton_Flavour->at(jt)); 
	
      }
    }// for(jt) All Jets
    
    //----------------------------------------------------------------------
    // Order the selected jets by pT. 
    //----------------------------------------------------------------------
    
    for (unsigned int i=0;i<nJets;i++){
      for (unsigned int j=i;j<nJets;j++){
	if(Pt_clean[i]<Pt_clean[j]){
	  
	  float tempPt;
	  tempPt=Pt_clean[i];
	  Pt_clean[i]=Pt_clean[j];
	  Pt_clean[j]=tempPt;
	  
	  float tempN;
	  tempN=JetNumber[i];
	  JetNumber[i]=JetNumber[j];
	  JetNumber[j]=tempN;
	  
	  float tempEta;
	  tempEta=Eta_clean[i];
	  Eta_clean[i]=Eta_clean[j];
	  Eta_clean[j]=tempEta;
	  
	  float tempVar;
	  tempVar=SysJetVar[i];
	  SysJetVar[i]=SysJetVar[j];
	  SysJetVar[j]=tempVar;
	  
	  float tempBtag;
	  tempBtag=Btag_clean[i];
	  Btag_clean[i]=Btag_clean[j];
	  Btag_clean[j]=tempBtag;
	  
	  float tempPF;
	  tempPF=PartonFlavour_clean[i];
	  PartonFlavour_clean[i]=PartonFlavour_clean[j];
	  PartonFlavour_clean[j]=tempPF;
	  
	}
      }
    }
    
    //--------------------------------------------------------------------------
    // b-tagging
    //--------------------------------------------------------------------------
    // official recipe (11/02/2015)
    // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation7TeVLegacy
    
    for (UInt_t l=0; l<nJets; l++) {

      if (sampleName.Contains("Data")){ 
	if (fBTagSF->IsTagged(Btag_clean[l], -999999, Pt_clean[l], Eta_clean[l]) ){ 
	  nJetsBtag++; 
	  BtagJetNumber.push_back(l);
	}
      } // if("Data")
      
      else if ( fBTagSF->IsTagged(Btag_clean[l], PartonFlavour_clean[l], Pt_clean[l], Eta_clean[l]) ){ 
	nJetsBtag++; 
	BtagJetNumber.push_back(l);
      }
    
    }// for(nJets) 

  } // if (mumuEvent || eeEvent || emuEvent)
  
  //--------------------------------------------------------------------------
  // LES: Propagation to the MET 
  //--------------------------------------------------------------------------
  if (sys_source == sys_les){
    for (UInt_t i=0; i<T_Elec_Energy->size(); i++){
      
      TLorentzVector Elec(T_Elec_Px    ->at(i),
                          T_Elec_Py    ->at(i),
                          T_Elec_Pz    ->at(i),
                          T_Elec_Energy->at(i));
      
      if (fabs(Elec.Eta()) < 1.5){ 
        if (sys_direction == sys_Up ) { 
	  metxSys += T_Elec_Px->at(i)*( 1 - (1+0.005) );
          metySys += T_Elec_Py->at(i)*( 1 - (1+0.005) ); 
	}
        else if (sys_direction == sys_Down ) { 
	  metxSys += T_Elec_Px->at(i)*( 1 - (1-0.005) );
          metySys += T_Elec_Py->at(i)*( 1 - (1-0.005) ); 
	}
      }
      else {
        if (sys_direction == sys_Up ) { 
	  metxSys += T_Elec_Px->at(i)*( 1 - (1+0.01) );
          metySys += T_Elec_Py->at(i)*( 1 - (1+0.01) ); 
	}
        else if (sys_direction == sys_Down ) { 
	  metxSys += T_Elec_Px->at(i)*( 1 - (1-0.01) );
          metySys += T_Elec_Py->at(i)*( 1 - (1-0.01) ); 
	}
      }
    }
  }//if(LES): MET Propagation

  //--------------------------------------------------------------------------
  // MET Propagation
  //--------------------------------------------------------------------------
  if (sys_source == sys_jes || sys_source == sys_jer || sys_source == sys_les){

    // comment if Pt only
    MET_Sys   = sqrt(metxSys*metxSys + metySys*metySys);
    missingEt = MET_Sys;
  }
  
  //--------------------------------------------------------------------------
  // Fill TopTrees - mumu channel
  //----------------------------------------------------------------------------
  if (PassTriggerMu() && mumuEvent) {
    
    // Fill Tree
    if (sampleName.Contains("TTbar") || sampleName.Contains("TTJets")){ 
      if( fileSuffix.Contains("Bkg") && nGenMuon!=2) FillTree(Muon1, Muon2, iMM);//TTbar Bkg
      if(!fileSuffix.Contains("Bkg") && nGenMuon==2) FillTree(Muon1, Muon2, iMM);//TTbar Signal
    }      
    else FillTree(Muon1, Muon2, iMM);
    
  } // if(mumu events)
  
  //----------------------------------------------------------------------------
  // Fill TopTrees - ee channel
  //----------------------------------------------------------------------------
  if (PassTriggerEG() && eeEvent) {
    
    //Tree
    if (sampleName.Contains("TTbar") || sampleName.Contains("TTJets")){ 
      if( fileSuffix.Contains("Bkg") && nGenElec!=2) FillTree(Elec1, Elec2, iEE);//TTbar Bkg
      if(!fileSuffix.Contains("Bkg") && nGenElec==2) FillTree(Elec1, Elec2, iEE);//TTbar Signal
    }      
    else FillTree(Elec1, Elec2, iEE);
  } // if(ee events)
  
  //----------------------------------------------------------------------------
  // Fill TopTrees - em channel
  //----------------------------------------------------------------------------
  if (PassTriggerEM() && emuEvent){    
    //Tree
    if (sampleName.Contains("TTbar") || sampleName.Contains("TTJets")){ 
      if( fileSuffix.Contains("Bkg") && !(nGenMuon==1 && nGenElec==1)) FillTree(Elec1, Muon1, iEM);//TTbar Bkg
      if(!fileSuffix.Contains("Bkg") &&  (nGenMuon==1 && nGenElec==1)) FillTree(Elec1, Muon1, iEM);//TTbar Signal
    }      
    
    else FillTree(Elec1, Muon1, iEM);
  } // if(mue event)
  
}// void(InsideLoop)

//------------------------------------------------------------------------------
// SetDataMembersAtTermination
//------------------------------------------------------------------------------
void TopLegacy::SetDataMembersAtTermination(){
  GetParameters();
}

//------------------------------------------------------------------------------
// Summary
//------------------------------------------------------------------------------
void TopLegacy::Summary(){
  
  AnalysisTree = ((TTree*) FindOutput("AnalysisTree"));
  
  cout << "---------------------------------------------------" <<endl;
  cout << "Number of entries in the Tree= " << AnalysisTree->GetEntries() <<endl;  
  cout << "---------------------------------------------------" <<endl;
  
}

//------------------------------------------------------------------------------
// GetParameters
//------------------------------------------------------------------------------
void TopLegacy::GetParameters()
{
  sampleName   = GetInputParameters()->TheNamedString("sampleName");
  fileSuffix   = GetInputParameters()->TheNamedString("fileSuffix");
  GetInputParameters()->TheNamedInt("sys_source", sys_source);
  GetInputParameters()->TheNamedInt("sys_direction", sys_direction);

  GetInputParameters()->TheNamedDouble("weight", fWeight); // cross section * luminosity / events in the sample
  GetInputParameters()->TheNamedDouble("luminosity", luminosity);

  //------------------------------------------------------------------------------
  // PU Reweight
  //------------------------------------------------------------------------------
  if(luminosity > 10000) fPUWeight = new PUWeight(luminosity,Summer12_53X,"2012"); // 19468.3pb-1
  else                   fPUWeight = new PUWeight(luminosity,Summer11_53X,"2011"); // 5100pb-1
  


      
  //--------------------------------------------------------------------------
  // New b-tagging code (11/02/2015)}
  //--------------------------------------------------------------------------
  int btagSysPar=0;
  if ( sys_source == sys_btag){
    if (sys_direction == sys_Up)   btagSysPar=1;
    if (sys_direction == sys_Down) btagSysPar=-1;
  }
  fBTagSF = new BTagSFUtil("CSV", "Medium", btagSysPar); 

}

//------------------------------------------------------------------------------
// PassTriggerMu
//------------------------------------------------------------------------------
Bool_t TopLegacy::PassTriggerMu()
{
  Bool_t pass = T_passTriggerDoubleMu; 
  if (T_Event_RunNumber==191090 ||  T_Event_RunNumber==193112 || T_Event_RunNumber==193116) pass=false; 
  return pass;
}

//------------------------------------------------------------------------------
// PassTriggerEG
//------------------------------------------------------------------------------
Bool_t TopLegacy::PassTriggerEG()
{
  Bool_t pass = T_passTriggerDoubleEl; 
  if (T_Event_RunNumber==191090 ||  T_Event_RunNumber==193112 || T_Event_RunNumber==193116) pass=false; 
  return pass;
}

//------------------------------------------------------------------------------
// PassTriggerEM
//------------------------------------------------------------------------------
Bool_t TopLegacy::PassTriggerEM()
{
  Bool_t pass = T_passTriggerElMu; 
  if (T_Event_RunNumber==191090 ||  T_Event_RunNumber==193112 || T_Event_RunNumber==193116) pass=false; 
  return pass;
}

//------------------------------------------------------------------------------
// SelectedVertexIndex
//------------------------------------------------------------------------------
Int_t TopLegacy::SelectedVertexIndex()
{
  Int_t goodVertexIndex = -999;
  nGoodVertex           =    0;
  
  for (UInt_t iVertex=0; iVertex<T_Vertex_z->size(); iVertex++) {
    
    if (fabs(T_Vertex_z ->at(iVertex)) < 24 &&
	T_Vertex_rho    ->at(iVertex)  <  2 &&
	T_Vertex_ndof   ->at(iVertex)  >  4 &&
	!T_Vertex_isFake->at(iVertex)) {
      
      nGoodVertex++;
      
      if (nGoodVertex == 1) goodVertexIndex = iVertex;
    }
  }
  
  return goodVertexIndex;
}

//------------------------------------------------------------------------------
// SelectedPfMuonPt
//------------------------------------------------------------------------------

// Tight Muon (TOP-13-004)
Double_t TopLegacy::SelectedTightPfMuonPt(UInt_t iMuon,
					  Int_t  iVertex)
{
  Double_t muonPt = -999;
  
  if (iVertex < 0) return muonPt;
  
  Bool_t pass = true;
  
  TLorentzVector Muon(T_Muon_Px    ->at(iMuon),
		      T_Muon_Py    ->at(iMuon),
		      T_Muon_Pz    ->at(iMuon),
		      T_Muon_Energy->at(iMuon));
    

  pass &=(T_Muon_IsGlobalMuon ->at(iMuon));
  pass &=(T_Muon_NormChi2GTrk->at(iMuon) < 10);
  pass &=(T_Muon_NValidHitsGTrk->at(iMuon) > 0);
  pass &=(T_Muon_NumOfMatchedStations->at(iMuon)>1);
  pass &=(T_Muon_IsTrackerMuonArbitrated->at(iMuon));
  pass &=(T_Muon_IPwrtAveBSInTrack->at(iMuon)< 0.2); //mm
  pass &=((T_Muon_vz->at(iMuon) - T_Vertex_z->at(iVertex))<0.5); //mm
  pass &=(T_Muon_NValidPixelHitsInTrk->at(iMuon)>0);
  pass &=(T_Muon_NLayers->at(iMuon)>5); 

  pass &= (Muon.Pt() > 20);
  pass &= (fabs(Muon.Eta()) < 2.4);
  
  Double_t relIso = (T_Muon_chargedHadronIsoR04->at(iMuon) + 
		     max(0.0 , 
			 T_Muon_neutralHadronIsoR04->at(iMuon) + 
			 T_Muon_photonIsoR04->at(iMuon) - 
			 0.5*T_Muon_sumPUPtR04->at(iMuon))
		     ) / Muon.Pt();
  pass &= (relIso < 0.12);

  if (pass) muonPt = Muon.Pt();

  return muonPt;
  
}

//------------------------------------------------------------------------------
// SelectedPfElecPt
//------------------------------------------------------------------------------
Double_t TopLegacy::SelectedPfElecPt(UInt_t iElec,
				     Int_t  iVertex)
{
  Double_t electronPt = -999;

  if (iVertex < 0) return electronPt;

  Bool_t pass = true;

  TLorentzVector Elec(T_Elec_Px    ->at(iElec),
		      T_Elec_Py    ->at(iElec),
		      T_Elec_Pz    ->at(iElec),
		      T_Elec_Energy->at(iElec));
  

  pass &= (Elec.Pt() > 20); // pT
  pass &= (fabs(Elec.Eta()) < 2.4); // eta
  pass &= (fabs(Elec.Eta()) < 1.4442 || fabs(Elec.Eta()) > 1.566); // eta

  // Effective Area Parametrization
  Double_t AEff03 = 0.;
  if      (fabs(Elec.Eta()) < 1.0)                              AEff03 = 0.13; // +/- 0.001
  else if (fabs(Elec.Eta()) >= 1.0 && fabs(Elec.Eta()) < 1.479) AEff03 = 0.14; // +/- 0.002
  else if (fabs(Elec.Eta()) >= 1.479 && fabs(Elec.Eta()) < 2.0) AEff03 = 0.07; // +/- 0.001
  else if (fabs(Elec.Eta()) >= 2.0 && fabs(Elec.Eta()) < 2.2)   AEff03 = 0.09; // +/- 0.001
  else if (fabs(Elec.Eta()) >= 2.2 && fabs(Elec.Eta()) < 2.3)   AEff03 = 0.11; // +/- 0.002
  else if (fabs(Elec.Eta()) >= 2.3 && fabs(Elec.Eta()) < 2.4)   AEff03 = 0.11; // +/- 0.003
  else if (fabs(Elec.Eta()) >= 2.4)                             AEff03 = 0.14; // +/- 0.004

  //////// dr 0.3 ////////
  Double_t relIso = (T_Elec_chargedHadronIso->at(iElec) + 
		     max(0.0, T_Elec_neutralHadronIso->at(iElec) + T_Elec_photonIso->at(iElec) - T_Event_RhoIso*AEff03)
		    )/Elec.Pt();

  pass &= (relIso < 0.1); 

  pass &= T_Elec_passConversionVeto->at(iElec);
  pass &= (T_Elec_nHits->at(iElec) < 1); 
  pass &= fabs(T_Elec_IPwrtPV->at(iElec)) < 0.04; // Impact parameter

  //----------------------------------------------------------------------------
  // Triggering MVA (pT>20GeV) 
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentification
  //----------------------------------------------------------------------------
  
  if (fabs(Elec.Eta())<0.8)        pass &= T_Elec_MVA->at(iElec) > 0.94; 
  else if (fabs(Elec.Eta())<1.479) pass &= T_Elec_MVA->at(iElec) > 0.85; 
  else                             pass &= T_Elec_MVA->at(iElec) > 0.62; 

  //----------------------------------------------------------------------------
  // Remove electrons close to GlobalMuons
  //----------------------------------------------------------------------------
  Double_t minDeltaR = 999;	

  TVector3 vElec(T_Elec_Px->at(iElec),
		 T_Elec_Py->at(iElec),
		 T_Elec_Pz->at(iElec)); 
      
  for (UInt_t j=0; j<T_Muon_Energy->size(); j++) {
   
    if ( T_Muon_IsGlobalMuon->at(j) ) { 
      
      TVector3 vMuon(T_Muon_Px->at(j),
		     T_Muon_Py->at(j),
		     T_Muon_Pz->at(j));
      
      Double_t deltaR = vElec.DeltaR(vMuon);
      
      if (deltaR < minDeltaR) minDeltaR = deltaR;
    }
  }  
  
  pass &= (minDeltaR > 0.1);
  
  if (pass) electronPt = Elec.Pt();

  return electronPt;
}

//------------------------------------------------------------------------------
// SelectedGenLepton
//------------------------------------------------------------------------------

void TopLegacy::SelectedGenLepton(){
  
  //----------------------------------------------------------------------------
  // Count generated muons and electrons
  //----------------------------------------------------------------------------
  nGenElec = T_Gen_ElecSt3_PID->size();
  nGenMuon = T_Gen_MuonSt3_PID->size();
  nGenTau  = T_Gen_TauSt3_PID->size();
  nGenLepton = nGenElec + nGenMuon + nGenTau;
  nTauElec = 0;
  nTauMuon = 0;

  for (UInt_t i=0; i<T_Gen_TauSt3_PID->size(); i++) {
    if (T_Gen_TauSt3_IsLepDec->at(i)) {
      if (abs(T_Gen_TauSt3_LepDec_PID->at(i)) == 11) nTauElec++;
      if (abs(T_Gen_TauSt3_LepDec_PID->at(i)) == 13) nTauMuon++;
    }
  }

  //----------------------------------------------------------------------------
  // Including taus as Signal
  //----------------------------------------------------------------------------
  if ( nGenLepton == 2 ) {
    if ( nGenElec == 2 || (nGenElec == 1 && nTauElec == 1) || nTauElec == 2 ) nGenElec = 2; 
    if ( nGenMuon == 2 || (nGenMuon == 1 && nTauMuon == 1) || nTauMuon == 2 ) nGenMuon = 2; 
    if ( (nGenElec == 1 && nGenMuon == 1) || (nGenElec == 1 && nTauMuon == 1) ||
	 (nGenMuon == 1 && nTauElec == 1) || (nTauElec == 1 && nTauMuon == 1) ) { nGenElec = 1; nGenMuon = 1; }
  } 
}

//------------------------------------------------------------------------------
// GetGenMuon
//------------------------------------------------------------------------------
void TopLegacy::GetGenMuon(){

  UInt_t muonSize = 0;
  
  muonSize = T_Gen_MuonSt3_Energy->size();
  
  for (UInt_t i=0; i<muonSize; i++) {
    
    TLorentzVector Muon(T_Gen_MuonSt3_Px    ->at(i),
			T_Gen_MuonSt3_Py    ->at(i),
			T_Gen_MuonSt3_Pz    ->at(i),
			T_Gen_MuonSt3_Energy->at(i));  
    
    if (Muon.Pt()>20 &&
	fabs(Muon.Eta()) < 2.4) { 
      
      if(T_Gen_MuonSt3_PID->at(i)== 13) Gen_Muon_Charge.push_back(-1.);
      if(T_Gen_MuonSt3_PID->at(i)==-13) Gen_Muon_Charge.push_back( 1.);
      
      Gen_Muon.push_back(Muon);
      
    }// if(GenMuon)
  }// for(muonSize)

//------------------------------------------------------------------------------
// GetGenMuon from Taus
//------------------------------------------------------------------------------
  for (UInt_t i=0; i<T_Gen_TauSt3_PID->size(); i++) {
    if (T_Gen_TauSt3_IsLepDec->at(i)) {
      if (abs(T_Gen_TauSt3_LepDec_PID->at(i)) == 13){
	
	
	TLorentzVector MuonTau(T_Gen_TauSt3_LepDec_Px    ->at(i),
			       T_Gen_TauSt3_LepDec_Py    ->at(i),
			       T_Gen_TauSt3_LepDec_Pz    ->at(i),
			       T_Gen_TauSt3_LepDec_Energy->at(i));  
	
	if (MuonTau.Pt()>20 &&
	    fabs(MuonTau.Eta()) < 2.4) {
	  
	  if(T_Gen_TauSt3_LepDec_PID->at(i)== 13) Gen_Muon_Charge.push_back(-1.);
	  if(T_Gen_TauSt3_LepDec_PID->at(i)==-13) Gen_Muon_Charge.push_back( 1.);
	  
	  Gen_Muon.push_back(MuonTau);
	  
	}// if(GenMuonTau)
	
	
      }// if(muon)
    }// if(LepDecay)
  }// for(Taus)
  
  nSGenMuon = Gen_Muon.size();
  
}

//------------------------------------------------------------------------------
// GetGenElec
//------------------------------------------------------------------------------
void TopLegacy::GetGenElec()
{
  UInt_t elecSize = 0;

  elecSize = T_Gen_ElecSt3_Energy->size();

  for (UInt_t i=0; i<elecSize; i++) {

    TLorentzVector Elec(T_Gen_ElecSt3_Px    ->at(i),
			T_Gen_ElecSt3_Py    ->at(i),
			T_Gen_ElecSt3_Pz    ->at(i),
			T_Gen_ElecSt3_Energy->at(i));
    
    if (Elec.Pt()>20 &&
	fabs(Elec.Eta()) < 2.4) { 
      
      if(T_Gen_ElecSt3_PID->at(i)== 11) Gen_Elec_Charge.push_back(-1.);
      if(T_Gen_ElecSt3_PID->at(i)==-11) Gen_Elec_Charge.push_back( 1.);
      
      Gen_Elec.push_back(Elec);
    }
  }
  
  //------------------------------------------------------------------------------
  // GetGenElec from Taus
  //------------------------------------------------------------------------------
  for (UInt_t i=0; i<T_Gen_TauSt3_PID->size(); i++) {
    if (T_Gen_TauSt3_IsLepDec->at(i)) {
      if (abs(T_Gen_TauSt3_LepDec_PID->at(i)) == 11){
	
	
	TLorentzVector ElecTau(T_Gen_TauSt3_LepDec_Px    ->at(i),
			       T_Gen_TauSt3_LepDec_Py    ->at(i),
			       T_Gen_TauSt3_LepDec_Pz    ->at(i),
			       T_Gen_TauSt3_LepDec_Energy->at(i));  
	
	if (ElecTau.Pt()>20 &&
	    fabs(ElecTau.Eta()) < 2.4) {
	  
	  if(T_Gen_TauSt3_LepDec_PID->at(i)== 11) Gen_Elec_Charge.push_back(-1.);
	  if(T_Gen_TauSt3_LepDec_PID->at(i)==-11) Gen_Elec_Charge.push_back( 1.);
	  
	  Gen_Elec.push_back(ElecTau);
	  
	}// if(GenMuonTau)
	
	
      }// if(muon)
    }// if(LepDecay)
  }// for(Taus)
  
  nSGenElec = Gen_Elec.size();
}


//------------------------------------------------------------------------------
// GetSelectedMuon
//------------------------------------------------------------------------------
void TopLegacy::GetSelectedMuon()
{
  UInt_t muonSize = 0;

  muonSize = T_Muon_Energy->size();

  for (UInt_t i=0; i<muonSize; i++) {
    
    if (SelectedTightPfMuonPt(i, SelectedVertexIndex()) > 0) { // TOP13-004
      
      S_Muon_Charge.push_back(T_Muon_Charge->at(i));
      
      float scale = 1.;
      if (sys_source == sys_les) { 
	if (sys_direction == sys_Up )   scale = 1.002; 
	if (sys_direction == sys_Down ) scale = 0.998; 
      }
      
      TLorentzVector Muon(scale*T_Muon_Px    ->at(i),
			  scale*T_Muon_Py    ->at(i),
			  scale*T_Muon_Pz    ->at(i),
			  scale*T_Muon_Energy->at(i));
      
      S_Muon.push_back(Muon);
      
    }
  }

  nSelMuon = S_Muon.size();
}

//------------------------------------------------------------------------------
// GetSelectedElec
//------------------------------------------------------------------------------
void TopLegacy::GetSelectedElec()
{
  UInt_t elecSize = 0;

  elecSize = T_Elec_Energy->size();

  for (UInt_t i=0; i<elecSize; i++) {

      if (SelectedPfElecPt(i, SelectedVertexIndex()) > 0) {
	
	S_Elec_Charge.push_back(T_Elec_Charge->at(i));
	S_Elec_SC_Et .push_back(T_Elec_SC_Et ->at(i));


        TLorentzVector ElecTmp(T_Elec_Px    ->at(i),
                               T_Elec_Py    ->at(i),
                               T_Elec_Pz    ->at(i),
                               T_Elec_Energy->at(i));
        float scale = 1.;
        if (sys_source == sys_les) {
	  if ( fabs(ElecTmp.Eta()) < 1.5 ) {
	    if (sys_direction == sys_Up ) scale = 1.005;
	    if (sys_direction == sys_Down ) scale = 0.995;
	  }
	  else {
	    if (sys_direction == sys_Up ) scale = 1.01;
	    if (sys_direction == sys_Down ) scale = 0.99;
	  }
        }
	
        TLorentzVector Elec(scale*T_Elec_Px    ->at(i),
                            scale*T_Elec_Py    ->at(i),
                            scale*T_Elec_Pz    ->at(i),
                            scale*T_Elec_Energy->at(i));
	
	S_Elec.push_back(Elec);
      }
  }
  
  nSelElec = S_Elec.size();
}

//------------------------------------------------------------------------------
// Initialise
//------------------------------------------------------------------------------

void TopLegacy::Initialise(){
  
  GetParameters();
  
  nEvents = InitCounterUI("nEvents","Number of analyzed events",0);
  
  //------------------------------------------------------------------------------
  // Tree Branches
  //------------------------------------------------------------------------------
  
  AnalysisTree = CreateTree("AnalysisTree","TopTree");
  
  AnalysisTree->Branch("TEvent",  &TEvent,  "TEvent/I");
  AnalysisTree->Branch("TWeight", &TWeight, "TWeight/F");
  AnalysisTree->Branch("TChannel",&TChannel,"TChannel/I");

  AnalysisTree->Branch("TNPV",&TNPV,"TNPV/I");
  AnalysisTree->Branch("TMET",&TMET,"TMET/F");
  AnalysisTree->Branch("TMET_Phi",&TMET_Phi,"TMET_Phi/F");
  AnalysisTree->Branch("TMETSig",&TMETSig,"TMETSig/F");

  AnalysisTree->Branch("TLep0Px",&TLep0Px,"TLep0Px/F");
  AnalysisTree->Branch("TLep0Py",&TLep0Py,"TLep0Py/F");
  AnalysisTree->Branch("TLep0Pz",&TLep0Pz,"TLep0Pz/F");
  AnalysisTree->Branch("TLep0E",&TLep0E,"TLep0E/F");
  AnalysisTree->Branch("TLep1Px",&TLep1Px,"TLep1Px/F");
  AnalysisTree->Branch("TLep1Py",&TLep1Py,"TLep1Py/F");
  AnalysisTree->Branch("TLep1Pz",&TLep1Pz,"TLep1Pz/F");
  AnalysisTree->Branch("TLep1E",&TLep1E,"TLep1E/F");

  AnalysisTree->Branch("TNJets",&TNJets,"TNJets/F");
  AnalysisTree->Branch("THT",&THT,"THT/F");
  AnalysisTree->Branch("TNJetsBtag",&TNJetsBtag,"TNJetsBtag/F");
  AnalysisTree->Branch("TBtagJet0",&TBtagJet0,"TBtagJet0/F");
  AnalysisTree->Branch("TBtagJet1",&TBtagJet1,"TBtagJet1/F");

  AnalysisTree->Branch("TJet0Px",&TJet0Px,"TJet0Px/F");
  AnalysisTree->Branch("TJet0Py",&TJet0Py,"TJet0Py/F");
  AnalysisTree->Branch("TJet0Pz",&TJet0Pz,"TJet0Pz/F");
  AnalysisTree->Branch("TJet0Et",&TJet0Et,"TJet0Et/F");
  AnalysisTree->Branch("TJet0E",&TJet0E,"TJet0E/F");
  AnalysisTree->Branch("TJet1Px",&TJet1Px,"TJet1Px/F");
  AnalysisTree->Branch("TJet1Py",&TJet1Py,"TJet1Py/F");
  AnalysisTree->Branch("TJet1Pz",&TJet1Pz,"TJet1Pz/F");
  AnalysisTree->Branch("TJet1Et",&TJet1Et,"TJet1Et/F");
  AnalysisTree->Branch("TJet1E",&TJet1E,"TJet1E/F");

  AnalysisTree->Branch("TBtagJet0Px",&TBtagJet0Px,"TBtagJet0Px/F");
  AnalysisTree->Branch("TBtagJet0Py",&TBtagJet0Py,"TBtagJet0Py/F");
  AnalysisTree->Branch("TBtagJet0Pz",&TBtagJet0Pz,"TBtagJet0Pz/F");
  AnalysisTree->Branch("TBtagJet0Et",&TBtagJet0Et,"TBtagJet0Et/F");
  AnalysisTree->Branch("TBtagJet0E",&TBtagJet0E,"TBtagJet0E/F");
  AnalysisTree->Branch("TBtagJet1Px",&TBtagJet1Px,"TBtagJet1Px/F");
  AnalysisTree->Branch("TBtagJet1Py",&TBtagJet1Py,"TBtagJet1Py/F");
  AnalysisTree->Branch("TBtagJet1Pz",&TBtagJet1Pz,"TBtagJet1Pz/F");
  AnalysisTree->Branch("TBtagJet1Et",&TBtagJet1Et,"TBtagJet1Et/F");
  AnalysisTree->Branch("TBtagJet1E",&TBtagJet1E,"TBtagJet1E/F");

  //-------------------------------------
  // Top pT reweight
  //-------------------------------------
  AnalysisTree->Branch("TtPx",&TtPx,"TtPx/F");
  AnalysisTree->Branch("TtPy",&TtPy,"TtPy/F");
  AnalysisTree->Branch("TtPz",&TtPz,"TtPz/F");
  AnalysisTree->Branch("TtE" ,&TtE ,"TtE/F");

  AnalysisTree->Branch("TtbarPx",&TtbarPx,"TtbarPx/F");
  AnalysisTree->Branch("TtbarPy",&TtbarPy,"TtbarPy/F");
  AnalysisTree->Branch("TtbarPz",&TtbarPz,"TtbarPz/F");
  AnalysisTree->Branch("TtbarE" ,&TtbarE ,"TtbarE/F");

  //-------------------------------------
  // JER and JES systematic variations
  //-------------------------------------
  AnalysisTree->Branch("TSysVarJet0",&TSysVarJet0,"TSysVarJet0/F");
  AnalysisTree->Branch("TSysVarJet1",&TSysVarJet1,"TSysVarJet1/F");
  AnalysisTree->Branch("TSysVarBtagJet0",&TSysVarBtagJet0,"TSysVarBtagJet0/F");
  AnalysisTree->Branch("TSysVarBtagJet1",&TSysVarBtagJet1,"TSysVarBtagJet1/F");
  
  AnalysisTree->Branch("TUncJet0",&TUncJet0,"TUncJet0/F");
  AnalysisTree->Branch("TUncJet1",&TUncJet1,"TUncJet1/F");
  AnalysisTree->Branch("TUncBtagJet0",&TUncBtagJet0,"TUncBtagJet0/F");
  AnalysisTree->Branch("TUncBtagJet1",&TUncBtagJet1,"TUncBtagJet1/F");
}

//------------------------------------------------------------------------------
// FillTrees with Dilepton Events
//------------------------------------------------------------------------------
void TopLegacy::FillTree(const TLorentzVector& lepton1, 
			 const TLorentzVector& lepton2,
			 UInt_t                iChannel){

  TEvent   = T_Event_EventNumber;
  TWeight  = weight;  
  
  //------------------------------------------------------------------------------
  //   ACCEPTANCE (Generation) ONLY ttbar Samples
  //------------------------------------------------------------------------------
  
  if(iChannel==3 ||  // mumu Gen
     iChannel==4 ||  // ee   Gen
     iChannel==5){   // emu or mue  Gen
    
    //-------------------------------------
    // Leptons
    //-------------------------------------
    if(lepton1.Pt()>lepton2.Pt()){
      TLep0Px = lepton1.Px();
      TLep0Py = lepton1.Py();
      TLep0Pz = lepton1.Pz();
      TLep0E  = lepton1.E();
      
      TLep1Px = lepton2.Px();
      TLep1Py = lepton2.Py();
      TLep1Pz = lepton2.Pz();
      TLep1E  = lepton2.E();
    }
    
    else{      
      TLep0Px = lepton2.Px();
      TLep0Py = lepton2.Py();
      TLep0Pz = lepton2.Pz();
      TLep0E  = lepton2.E();
      
      TLep1Px = lepton1.Px();
      TLep1Py = lepton1.Py();
      TLep1Pz = lepton1.Pz();
      TLep1E  = lepton1.E();
      
      if(iChannel==5) iChannel=-5; // emu=5 || mue=-5      
    }
    
    //-------------------------------------
    // MET
    //-------------------------------------
    TMET     = MET_Gen;
    TMET_Phi = T_METgen_Phi;
    TMETSig  = 0.0; // Is not relevant for Acc/Eff studies
    
    TNPV     = 0; // Is not relevant for Acc/Eff studies
    TChannel = iChannel;
    
    //-------------------------------------
    // Jets
    //-------------------------------------
    TNJets     = NGen_Jet.size();
    TNJetsBtag = NGen_b.size();
    THT        = 0; // Is not relevant for Acc/Eff studies
    
    // Jet initialization
    TJet0Px     = 0.0;
    TJet0Py     = 0.0;
    TJet0Pz     = 0.0;
    TJet0E      = 0.0;
    
    TJet1Px     = 0.0;
    TJet1Py     = 0.0;
    TJet1Pz     = 0.0;
    TJet1E      = 0.0;
    // b-jet initialization
    TBtagJet0Px = 0.0;
    TBtagJet0Py = 0.0;
    TBtagJet0Pz = 0.0;
    TBtagJet0E  = 0.0;
    
    TBtagJet1Px = 0.0;
    TBtagJet1Py = 0.0;
    TBtagJet1Pz = 0.0;
    TBtagJet1E  = 0.0;
    

    if(NGen_Jet.size()>0){
      TBtagJet0 = T_JetAKCHS_Parton_Flavour->at(NGen_Jet.at(0)); // Jet Flavour      
      
      TJet0Px   = T_JetAKCHS_GenJet_Px->at(NGen_Jet.at(0));
      TJet0Py   = T_JetAKCHS_GenJet_Py->at(NGen_Jet.at(0));
      TJet0Pz   = T_JetAKCHS_GenJet_Pz->at(NGen_Jet.at(0));
      TJet0E    = T_JetAKCHS_GenJet_Energy->at(NGen_Jet.at(0));      
    }
    
    if(NGen_Jet.size()>1){
      TBtagJet1 = T_JetAKCHS_Parton_Flavour->at(NGen_Jet.at(1)); // Jet Flavour      
      
      TJet1Px   = T_JetAKCHS_GenJet_Px->at(NGen_Jet.at(1));
      TJet1Py   = T_JetAKCHS_GenJet_Py->at(NGen_Jet.at(1));
      TJet1Pz   = T_JetAKCHS_GenJet_Pz->at(NGen_Jet.at(1));
      TJet1E    = T_JetAKCHS_GenJet_Energy->at(NGen_Jet.at(1));
    } 
    
    if(NGen_b.size()>0){
      TBtagJet0Px = T_JetAKCHS_GenJet_Px->at(NGen_b.at(0));
      TBtagJet0Py = T_JetAKCHS_GenJet_Py->at(NGen_b.at(0));
      TBtagJet0Pz = T_JetAKCHS_GenJet_Pz->at(NGen_b.at(0));
      TBtagJet0E  = T_JetAKCHS_GenJet_Energy->at(NGen_b.at(0));
      
    }
    
    if(NGen_b.size()>1){
      TBtagJet1Px = T_JetAKCHS_GenJet_Px->at(NGen_b.at(1));
      TBtagJet1Py = T_JetAKCHS_GenJet_Py->at(NGen_b.at(1));
      TBtagJet1Pz = T_JetAKCHS_GenJet_Pz->at(NGen_b.at(1));
      TBtagJet1E  = T_JetAKCHS_GenJet_Energy->at(NGen_b.at(1));
    }
    
    // Not relevant variables for Acc/Eff studies
    TJet0Et         = 0.0;
    TUncJet0        = 0.0;    
    TJet1Et         = 0.0;
    TUncJet1        = 0.0;
    
    TBtagJet0Et     = 0.0;
    TUncBtagJet0    = 0.0;
    TBtagJet1Et     = 0.0;
    TUncBtagJet1    = 0.0;
    
    TSysVarJet0     = 0.0;
    TSysVarJet1     = 0.0;
    TSysVarBtagJet0 = 0.0;
    TSysVarBtagJet1 = 0.0;    
    
    if(T_Gen_tSt3_PID->at(0) == 6){
      
      TtPx =T_Gen_tSt3_Px ->at(0);
      TtPy =T_Gen_tSt3_Py ->at(0);
      TtPz =T_Gen_tSt3_Pz ->at(0);
      TtE  =T_Gen_tSt3_Energy ->at(0);
      
      TtbarPx =T_Gen_tSt3_Px ->at(1);
      TtbarPy =T_Gen_tSt3_Py ->at(1);
      TtbarPz =T_Gen_tSt3_Pz ->at(1);
      TtbarE  =T_Gen_tSt3_Energy ->at(1);
      
    }
    
    else{
      
      TtPx =T_Gen_tSt3_Px ->at(1);
      TtPy =T_Gen_tSt3_Py ->at(1);
      TtPz =T_Gen_tSt3_Pz ->at(1);
      TtE  =T_Gen_tSt3_Energy ->at(1);
      
      TtbarPx =T_Gen_tSt3_Px ->at(0);
      TtbarPy =T_Gen_tSt3_Py ->at(0);
      TtbarPz =T_Gen_tSt3_Pz ->at(0);
      TtbarE  =T_Gen_tSt3_Energy ->at(0);
      
    }      
  } // if(ichannel == 3,4,5,-5)
  

  //------------------------------------------------------------------------------
  // RECO events
  //------------------------------------------------------------------------------
  
  else{

    TNPV     = nGoodVertex;

    //-------------------------------------
    // MET
    //-------------------------------------
    TMET     = missingEt;
    //TMET_Phi = T_METType1CorrectedPFMetSmeared_Phi; // MC
    TMET_Phi = T_METPFTypeI_Phi; // DATA
    TMETSig  = T_METPF_Sig;

    //-------------------------------------
    // pT Reweight for ttbar signal samples
    //-------------------------------------
    TtPx     = 0.0;
    TtPy     = 0.0;
    TtPz     = 0.0;
    TtE      = 0.0;
    
    TtbarPx  = 0.0;
    TtbarPy  = 0.0;
    TtbarPz  = 0.0;
    TtbarE   = 0.0;
    
    if((sampleName.Contains("MadSpin") || 
	sampleName.Contains("TTbar_Madgraph")) && 
       !fileSuffix.Contains("Bkg")){
      
      if(T_Gen_tSt3_PID->at(0) == 6){
	
      	TtPx =T_Gen_tSt3_Px ->at(0);
      	TtPy =T_Gen_tSt3_Py ->at(0);
      	TtPz =T_Gen_tSt3_Pz ->at(0);
      	TtE  =T_Gen_tSt3_Energy ->at(0);
	
      	TtbarPx =T_Gen_tSt3_Px ->at(1);
      	TtbarPy =T_Gen_tSt3_Py ->at(1);
      	TtbarPz =T_Gen_tSt3_Pz ->at(1);
      	TtbarE  =T_Gen_tSt3_Energy ->at(1);
	
      }
      
      else{
	
      	TtPx =T_Gen_tSt3_Px ->at(1);
      	TtPy =T_Gen_tSt3_Py ->at(1);
      	TtPz =T_Gen_tSt3_Pz ->at(1);
      	TtE  =T_Gen_tSt3_Energy ->at(1);
	
      	TtbarPx =T_Gen_tSt3_Px ->at(0);
      	TtbarPy =T_Gen_tSt3_Py ->at(0);
      	TtbarPz =T_Gen_tSt3_Pz ->at(0);
      	TtbarE  =T_Gen_tSt3_Energy ->at(0);
	
      }      
    } //if(top reweight)

    //-------------------------------------
    // Leptons
    //-------------------------------------
    if(lepton1.Pt()>lepton2.Pt()){
      
      TLep0Px = lepton1.Px();
      TLep0Py = lepton1.Py();
      TLep0Pz = lepton1.Pz();
      TLep0E  = lepton1.E();
      
      TLep1Px = lepton2.Px();
      TLep1Py = lepton2.Py();
      TLep1Pz = lepton2.Pz();
      TLep1E  = lepton2.E();
      
    }
    
    else{
      
      TLep0Px = lepton2.Px();
      TLep0Py = lepton2.Py();
      TLep0Pz = lepton2.Pz();
      TLep0E  = lepton2.E();
      
      TLep1Px = lepton1.Px();
      TLep1Py = lepton1.Py();
      TLep1Pz = lepton1.Pz();
      TLep1E  = lepton1.E();
      
      if(iChannel==2) iChannel=-2; // emu=2 || mue=-2
      
    }
    
    
    TChannel   = iChannel;
    
    //-------------------------------------
    // Jets
    //-------------------------------------
    TNJets     = nJets;
    THT        = HT;
    TNJetsBtag = nJetsBtag;
    
    // Jet initialization
    TJet0Px    = 0.0;
    TJet0Py    = 0.0;
    TJet0Pz    = 0.0;
    TJet0Et    = 0.0;
    TJet0E     = 0.0;
    TBtagJet0  = 0.0;
    TUncJet0   = 0.0;
 
    TJet1Px    = 0.0;
    TJet1Py    = 0.0;
    TJet1Pz    = 0.0;
    TJet1Et    = 0.0;
    TJet1E     = 0.0;
    TBtagJet1  = 0.0;
    TUncJet1   = 0.0;


    if(nJets>0){

      int Njet0=JetNumber[0];
      
      TJet0Px     = T_JetAKCHS_Px     ->at(Njet0);
      TJet0Py     = T_JetAKCHS_Py     ->at(Njet0);
      TJet0Pz     = T_JetAKCHS_Pz     ->at(Njet0);
      TJet0Et     = T_JetAKCHS_Et     ->at(Njet0);
      TJet0E      = T_JetAKCHS_Energy ->at(Njet0);
      
      TBtagJet0   = Btag_clean[0];
      TSysVarJet0 = SysJetVar[0];
      //TUncJet0    = T_JetAKCHS_Corr->at(Njet0);
      TUncJet0    = 0.0;

   }//if(nJets>0)
    
    if(nJets>1){
      int Njet1   = JetNumber[1];
      
      TJet1Px     = T_JetAKCHS_Px     ->at(Njet1);
      TJet1Py     = T_JetAKCHS_Py     ->at(Njet1);
      TJet1Pz     = T_JetAKCHS_Pz     ->at(Njet1);
      TJet1Et     = T_JetAKCHS_Et     ->at(Njet1);
      TJet1E      = T_JetAKCHS_Energy ->at(Njet1);

      TBtagJet1   = Btag_clean[1];
      TSysVarJet1 = SysJetVar[1];
      //TUncJet1    = T_JetAKCHS_Corr->at(Njet1);
      TUncJet1    = 0.0;
    }//if(nJets>1)

    TBtagJet0Px  = 0.0;
    TBtagJet0Py  = 0.0;
    TBtagJet0Pz  = 0.0;
    TBtagJet0Et  = 0.0;
    TBtagJet0E   = 0.0;
    TUncBtagJet0 = 0.0;

    TBtagJet1Px  = 0.0;
    TBtagJet1Py  = 0.0;
    TBtagJet1Pz  = 0.0;
    TBtagJet1Et  = 0.0;
    TBtagJet1E   = 0.0;
    TUncBtagJet1 = 0.0;

    if(nJetsBtag>0){
      int Nbtagjet0 = JetNumber[BtagJetNumber[0]];
      
      TBtagJet0Px   = T_JetAKCHS_Px     ->at(Nbtagjet0);
      TBtagJet0Py   = T_JetAKCHS_Py     ->at(Nbtagjet0);
      TBtagJet0Pz   = T_JetAKCHS_Pz     ->at(Nbtagjet0);
      TBtagJet0Et   = T_JetAKCHS_Et     ->at(Nbtagjet0);
      TBtagJet0E    = T_JetAKCHS_Energy ->at(Nbtagjet0);

      TSysVarBtagJet0 = SysJetVar[BtagJetNumber[0]];
      //TUncBtagJet0    = T_JetAKCHS_Corr->at(Nbtagjet0);
      TUncBtagJet0    = 0.0;
    }
    if(nJetsBtag>1){
      int Nbtagjet1=JetNumber[BtagJetNumber[1]];
      
      TBtagJet1Px = T_JetAKCHS_Px     ->at(Nbtagjet1);
      TBtagJet1Py = T_JetAKCHS_Py     ->at(Nbtagjet1);
      TBtagJet1Pz = T_JetAKCHS_Pz     ->at(Nbtagjet1);
      TBtagJet1Et = T_JetAKCHS_Et     ->at(Nbtagjet1);
      TBtagJet1E  = T_JetAKCHS_Energy ->at(Nbtagjet1);

      TSysVarBtagJet1 = SysJetVar[BtagJetNumber[1]];
      //TUncBtagJet1    = T_JetAKCHS_Corr->at(Nbtagjet1);
      TUncBtagJet1    = 0.0;
    }
        
    } //Else (Generation)

    AnalysisTree->Fill();
    
}



