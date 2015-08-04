////////////////////////////////////////////////////////////////////////////////
//
//    FILE: RunTopLegacy.C
// AUTHORS: Brochero, Top Group
//    DATE: 2014
//
// Before run, take into account:
// 
// DATA: 
// - Uncomment the first code lines of TopLegacy.C to include the   
//   definition of the GEN branches.
// - Use the correct MET variable. In data the available one is METPFTypeI.
// - Also change T_METPFTypeI_Phi in the FillTree
// 
// MC:
// - Comment the first code lines of TopLegacy.C.
// - Select the correct MET variable.
// - For ttbar, VV and tW, uncomment the lines asociated to the systematic 
//   uncertainties (if you are going to produce the syst. variations).
//
////////////////////////////////////////////////////////////////////////////////

TProof* proof = 0;

Double_t G_Event_Weight = 1;

void RunTopLegacy(int isample=0, // Temporal solution to fix the prblems with the PROOF sessions
		  TString  sampleName     = "TTbar_Powheg",
		  TString  selector       = "TopLegacy",
		  TString  fileSuffix     = "",
		  TString  sampleSys      = "", // for systematic samples, indicate which source is going to be studied 
		  Double_t lumiForPUdata  = 5100.0, // luminosity in http://www.hep.uniovi.es/jfernan/PUhistos
		  Int_t    sys_source     = -1,      // -1 nominal; enum sys_sources { sys_jes = 0, sys_jer, sys_les, sys_unenergy };
		  Int_t    sys_direction  = -1,      // -1 nominal; enum sys_directions { sys_Up = 0, sys_Nom, sys_Down};
		  Long64_t nEvents        = -1){
  
  gROOT->LoadMacro("$PAFPATH/PAF.C");

  Double_t G_Event_Lumi=0.0;
  if(lumiForPUdata>10000) G_Event_Lumi = 19468.3;   
  else G_Event_Lumi = 5100.0;

  cout << "Params: " << endl;
  cout << "sampleName      " << sampleName    << endl;
  cout << "selector        " << selector      << endl;
  cout << "fileSuffix      " << fileSuffix    << endl;
  cout << "sampleSys       " << sampleSys     << endl;
  cout << "lumiForNorm     " << G_Event_Lumi  << endl;
  cout << "lumiForPUdata   " << lumiForPUdata << endl;
  cout << "sys_source      " << sys_source    << endl;
  cout << "sys_direction   " << sys_direction << endl;
  cout << "nEvents         " << nEvents       << endl;

  //----------------------------------------------------------------------------
  // Base path to input files
  //----------------------------------------------------------------------------
  TString dataPath = "/gpfs/csic_projects/tier3data";     // IFCA   (gridui)
  
  // Read Data and MC information from the goBogle doc table
  gROOT->LoadMacro("/gpfs/csic_users/brochero/DatasetManager/DatasetManager.C+");

  
  if(G_Event_Lumi>10000) DatasetManager* dm = new DatasetManager("Legacy_Summer15_74X"); // 19468.3pb-1
  else                   DatasetManager* dm = new DatasetManager("Legacy_Summer11_53X"); // 5100pb-1
  
  // Use this if you know that the information on the google doc table has
  // changed and you need to update the information
  dm->RedownloadFiles();
  dm->LoadDataset(sampleName);  // Load information about a given dataset
  G_Event_Weight = dm->GetCrossSection() * G_Event_Lumi / dm->GetEventsInTheSample(); // Weight estimation
  
  cout << "\n\n" << endl;
  cout << "------------------------------------------------\n" << endl;
  cout << "      x-section = " << dm->GetCrossSection()      << endl;
  cout << "     luminosity = " << G_Event_Lumi               << endl;
  cout << "        nevents = " << dm->GetEventsInTheSample() << endl;
  cout << " base file name = " << dm->GetBaseFileName()      << endl;
  cout << " G_Event_Weight = " << G_Event_Weight             << endl;
  cout << "------------------------------------------------\n" << endl;
  cout << "\n\n" << endl;
  
 
  //----------------------------------------------------------------------------
  // Tree type
  //----------------------------------------------------------------------------
  gPAFOptions->SetTreeType(kMiniTrees);

  //----------------------------------------------------------------------------
  // PROOF mode
  //----------------------------------------------------------------------------
  //gPAFOptions->proofMode = kSequential;       // No PROOF

  gPAFOptions->proofMode = kLite;             // PROOF Lite

  //gPAFOptions->proofMode = kCluster;            // PROOF Cluster
  //gPAFOptions->NSlots = 10;                   // Number of slots

  // Optional parameters for PROOF Cluster mode
  // gPAFOptions->proofServer = "proof.ifca.es";  // PROOF server
  // gPAFOptions->proofServerPort = 1093;         // PROOF port
  // gPAFOptions->maxSlavesPerNode = 9999;        // Max number of slaves / node
 
  //----------------------------------------------------------------------------
  // PROOF start
  //----------------------------------------------------------------------------
  cout << ">> Starting PROOF..." << endl;
  proof = InitProof();
  if (!proof && gPAFOptions->proofMode != kSequential) {
    cerr << "ERROR: I could not initialise a PROOF session!" << endl;
    return;
  }
  if (gPAFOptions->proofMode != kSequential) gPAFOptions->proofSession->SetLogLevel(2, TProofDebug::kOutput); 

  //for(int isample=0; isample < dm->GetFiles().size(); isample++){

  //----------------------------------------------------------------------------
  // Loaded Sample
  //----------------------------------------------------------------------------

    cout << "\n\n" << endl;
    cout << "------------------------------------------------\n" << endl;
    cout << " Path file                      = " << dm->GetLocalBasePath()  << endl;
    cout << " Total number of files          = " << dm->GetFiles().size()   << endl;
    cout << " Files processing               = " << isample                 << endl;
    cout << "------------------------------------------------\n" << endl;
    cout << "\n\n" << endl;
    
    //----------------------------------------------------------------------------
    // Data Temporal Solution!!!!!!!
    //----------------------------------------------------------------------------
    // if(G_Event_Lumi>10000){// 8TeV
    //   if(sampleName=="DataMuA") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_DoubleMuA_876");
    //   if(sampleName=="DataMuB") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_DoubleMuParkedB_4412");
    //   if(sampleName=="DataMuC") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_DoubleMuParkedC_7016");
    //   if(sampleName=="DataMuD") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_DoubleMuParkedD_7360");

    //   if(sampleName=="DataMuEGA") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_MuEGA_876");    
    //   if(sampleName=="DataMuEGB") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_MuEGB_4412");    
    //   if(sampleName=="DataMuEGC") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_MuEGC_7016");    
    //   if(sampleName=="DataMuEGD") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_MuEGD_7360");    

    //   if(sampleName=="DataEGA") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_DoubleElectronA_876");    
    //   if(sampleName=="DataEGB") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_DoubleElectronB_4412");    
    //   if(sampleName=="DataEGC") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_DoubleElectronC_7016");    
    //   if(sampleName=="DataEGD") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer12_53X/Legacy","Tree_DoubleElectronD_7360");    

    // }
    // else{// 7TeV
    // if(sampleName=="DataMuA") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer11_53X","Tree_DoubleMuA");
    // if(sampleName=="DataMuB") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer11_53X","Tree_DoubleMuB");
    
    // if(sampleName=="DataEGA") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer11_53X","Tree_DoubleEleA");
    // if(sampleName=="DataEGB") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer11_53X","Tree_DoubleEleB");
    
    // if(sampleName=="DataMuEGA") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer11_53X","Tree_MuEGA");
    // if(sampleName=="DataMuEGB") 
    // 	vector<TString> realdata= DatasetManager::GetRealDataFiles("MC_Summer11_53X","Tree_MuEGB");
    // }

    // gPAFOptions->dataFiles.push_back(realdata[isample]);
    // G_Event_Weight=1.0;
    
    //----------------------------------------------------------------------------
    //Monte Carlo
    //----------------------------------------------------------------------------
    gPAFOptions->dataFiles.push_back((dm->GetFiles())[isample]);


  //----------------------------------------------------------------------------
  // Output file name
  //----------------------------------------------------------------------------
  TString outputDir = "/gpfs/csic_projects/cms/brochero/TopResults_Legacy_7TeV";
  gSystem->mkdir(outputDir, kTRUE);
  
  std::ostringstream oss;      
  oss << G_Event_Lumi;

  TString LumiString = oss.str();
  TString outputFile = outputDir
    + "/"
    + "TrLe_" 
    + "v2btag_"
    + LumiString
    + "pb-1_"
    + sampleName
    + sampleSys
    + fileSuffix;

  char filenumber[200];
  sprintf(filenumber,"%i",isample);

  if(dm->GetFiles().size()==1) outputFile = outputFile + ".root";
  else                         outputFile = outputFile + filenumber + ".root";  

  gPAFOptions->outputFile = outputFile;

  //----------------------------------------------------------------------------
  // Parameters for the analysis
  //----------------------------------------------------------------------------
  // See packages/InputParameters/InputParameters.h for information on how
  // to use this class.

  gPAFOptions->inputParameters = new InputParameters();

  gPAFOptions->inputParameters->SetNamedString("sampleName",     sampleName.Data());
  gPAFOptions->inputParameters->SetNamedDouble("weight",         G_Event_Weight);
  gPAFOptions->inputParameters->SetNamedInt   ("sys_source",     sys_source);
  gPAFOptions->inputParameters->SetNamedInt   ("sys_direction",  sys_direction);
  gPAFOptions->inputParameters->SetNamedDouble("luminosity",     lumiForPUdata);
  gPAFOptions->inputParameters->SetNamedString("fileSuffix",     fileSuffix.Data());


  //----------------------------------------------------------------------------
  // Definition of Gen variables in DATA
  //----------------------------------------------------------------------------
  if(sampleName.Contains("Data")) gSystem->AddIncludePath("-D__ISMC"); //Works?
  
  //----------------------------------------------------------------------------
  // Number of events (Long64_t)
  //----------------------------------------------------------------------------
  gPAFOptions->SetNEvents(nEvents);

  //----------------------------------------------------------------------------
  // First event (Long64_t)
  //----------------------------------------------------------------------------
  gPAFOptions->SetFirstEvent(0);
  
  //----------------------------------------------------------------------------
  // Name of analysis class
  //----------------------------------------------------------------------------
  // If 0 the default name schema will be used, i.e. depending on the value
  // of gPAFOptions->treeType: MyAnalysisTESCO or MyAnalsyisMiniTrees

  gPAFOptions->SetAnalysis(selector.Data());

  gPAFOptions->AddPackage("PUWeight");
  gPAFOptions->AddPackage("BTagSFUtil");

  //----------------------------------------------------------------------------
  // Control output and checks
  //----------------------------------------------------------------------------
  // + If true (default) the output file is reopened so the objects in the
  //   file can be interactively accessed. The object in the output are also
  //   listed

  gPAFOptions->reopenOutputFile = false;

  //----------------------------------------------------------------------------
  // Run the analysis
  //----------------------------------------------------------------------------
  if (!RunAnalysis())
    cerr << "ERROR: There was a problem running the analysis!" << endl;

  gPAFOptions->dataFiles.clear();

  //} // For over all samples
}
