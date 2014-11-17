#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <vector>
#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLegend.h"
#include <iomanip>

void reconstruction(){

  // Open the raw data file:
  unsigned int    events(-9999);
  unsigned int    run(-9999);
  unsigned int    tdcPhase(9999);

  std::vector<int> *hits_S0 = 0;
  std::vector<int> *hits_S1 = 0;
  std::vector<int> *hits_S2 = 0;
  std::vector<int> *hits_S3 = 0;

  TFile *dutFile = TFile::Open("../data/RAW_USC.00000478.0001.A.storageManager.00.0000.root");
  TTree *duttree = (TTree *) dutFile->Get("tree");
  duttree->SetBranchAddress("event", &events);
  duttree->SetBranchAddress("run", &run);
  duttree->SetBranchAddress("dut0_row", &hits_S0);
  duttree->SetBranchAddress("dut1_row", &hits_S1);
  duttree->SetBranchAddress("dut2_row", &hits_S2);
  duttree->SetBranchAddress("dut3_row", &hits_S3);
  duttree->SetBranchAddress("tdcPhase", &tdcPhase);


  // Book all histograms:
  TH1F *hitmap_S0 = new TH1F("hitmap_S0","hitmap_S0", 256, 0. , 256.);
  TH1F *hitmap_S1 = new TH1F("hitmap_S1","hitmap_S1", 256, 0. , 256.);
  TH1F *hitmap_S2 = new TH1F("hitmap_S2","hitmap_S2", 256, 0. , 256.);
  TH1F *hitmap_S3 = new TH1F("hitmap_S3","hitmap_S3", 256, 0. , 256.);

  TH1F *hitmap_S0_masked = new TH1F("hitmap_S0_masked","hitmap_S0_masked", 256, 0. , 256.);
  TH1F *hitmap_S1_masked = new TH1F("hitmap_S1_masked","hitmap_S1_masked", 256, 0. , 256.);
  TH1F *hitmap_S2_masked = new TH1F("hitmap_S2_masked","hitmap_S2_masked", 256, 0. , 256.);
  TH1F *hitmap_S3_masked = new TH1F("hitmap_S3_masked","hitmap_S3_masked", 256, 0. , 256.);

  TH1F *h_clusters_S2 = new TH1F("clusters_S2","clusters_S2", 256, 0. , 256.);
  TH1F *h_clusters_S3 = new TH1F("clusters_S3","clusters_S3", 256, 0. , 256.);


  // Prepare output Reco ROOT Tree:
  std::vector<int> clusterwidth_S0;
  std::vector<double> clusterposition_S0;
  std::vector<int> clusterwidth_S1;
  std::vector<double> clusterposition_S1;
  std::vector<int> clusterwidth_S2;
  std::vector<double> clusterposition_S2;
  std::vector<int> clusterwidth_S3;
  std::vector<double> clusterposition_S3;
  std::vector<double> stubs_M1;
  std::vector<double> stubs_M2;
  std::vector<int> hits_valid_S0;
  std::vector<int> hits_valid_S1;
  std::vector<int> hits_valid_S2;
  std::vector<int> hits_valid_S3;

  TFile *recofile = new TFile("../data/RECO_USC.00000478.0001.A.storageManager.00.0000.root","RECREATE");
  TTree *recotree = new TTree("recotree","reconstructed");
  recotree->Branch("events" , &events );
  recotree->Branch("run" , &run );
  recotree->Branch("row_S0", &hits_S0);
  recotree->Branch("row_S1", &hits_S1);
  recotree->Branch("row_S2", &hits_S2);
  recotree->Branch("row_S3", &hits_S3);
  recotree->Branch("row_S0_masked", &hits_valid_S0);
  recotree->Branch("row_S1_masked", &hits_valid_S1);
  recotree->Branch("row_S2_masked", &hits_valid_S2);
  recotree->Branch("row_S3_masked", &hits_valid_S3);
  recotree->Branch("clusters_S0", &clusterposition_S0);
  recotree->Branch("clusters_S1", &clusterposition_S1);
  recotree->Branch("clusters_S2", &clusterposition_S2);
  recotree->Branch("clusters_S3", &clusterposition_S3);
  recotree->Branch("clusterwidth_S0", &clusterwidth_S0);
  recotree->Branch("clusterwidth_S1", &clusterwidth_S1);
  recotree->Branch("clusterwidth_S2", &clusterwidth_S2);
  recotree->Branch("clusterwidth_S3", &clusterwidth_S3);
  recotree->Branch("stubs_module_1",&stubs_M1);
  recotree->Branch("stubs_module_2",&stubs_M2);
  recotree->Branch("tdc", &tdcPhase);


  // The main events loop:
  for (int i = 0; i < duttree->GetEntries(); i++) {
    // Some printout to see where we are:
    if(i%1000 == 0) std::cout << "Processing event " << i << std::endl;
    duttree->GetEntry(i);

    // Masking noisy strips
    // If a strip is "noisy" can be read from the hit map histograms
    // and results in a hardcoded cut on strip numbers.
    for(size_t hit = 0 ; hit < hits_S0->size(); hit++) hitmap_S0->Fill(hits_S0->at(hit));
    // ...

    // Hint: S0 and S1 have 254 strips each, S2, S3 have 256 strips each.

    // Store all strips which are fine in vectors for further analyses
    // hits_valid_S0;


    // Clustering
    // On each sensor, cluster all strips which belong together. 
    // The most simple solution is to just add all up which are next to each other
    // storing the first and last strip.

    std::vector<int> cluster_firstStrip_S0;
    std::vector<int> cluster_lastStrip_S0;

    std::vector<int> cluster_firstStrip_S1;
    std::vector<int> cluster_lastStrip_S1;

    std::vector<int> cluster_firstStrip_S2;
    std::vector<int> cluster_lastStrip_S2;

    std::vector<int> cluster_firstStripS_S3;
    std::vector<int> cluster_lastStrip_S3;

    // Run over all strips on every module to connect them:
    for(size_t k = 0; k < hits_valid_S0.size(); ++k){
      // ...
    }

    // Some extra cross check that we didn't mess things up:
    if(cluster_firstStrip_S0.size() != cluster_lastStrip_S0.size()){
      std::cerr << "Error in vector size cluster/S0 for first/last strip.\n";
      exit;
    }

    //Calculate the cluster width and position for all sensors:
    for(size_t strip = 0; strip < cluster_firstStrip_S0.size(); ++strip) {
      // ...
      // clusterwidth_S0
      // clusterposition_S0
    }
    // ...


    // Produce the stubs for the two modules
    // Stubs are defined between sensors S0 and S1 (module 1) and S2/S3 (module 2).
		
    // Stubs Module 1:
    // (restricting to one-cluster events)
    if(clusterposition_S0.size() == 1 && clusterposition_S1.size() == 1) {
      stubs_M1.push_back(clusterposition_S1.front() - clusterposition_S0.front());
    }
    // ...


    recotree->Fill();

    clusterposition_S0.clear();
    clusterposition_S1.clear();
    clusterposition_S2.clear();
    clusterposition_S3.clear();

    clusterwidth_S0.clear();
    clusterwidth_S1.clear();
    clusterwidth_S2.clear();
    clusterwidth_S3.clear();

    stubs_M1.clear();
    stubs_M2.clear();

    hits_valid_S0.clear();
    hits_valid_S1.clear();
    hits_valid_S2.clear();
    hits_valid_S3.clear();

  } //end of events loop


  recofile->cd();

  hitmap_S0->Write();
  hitmap_S1->Write();
  hitmap_S2->Write();
  hitmap_S3->Write();
  hitmap_S0_masked->Write();
  hitmap_S1_masked->Write();
  hitmap_S2_masked->Write();
  hitmap_S3_masked->Write();
  h_clusters_S2->Write();
  h_clusters_S3->Write();
  recotree->Write();

  recofile->Close();

  //h_position->Draw();
}

