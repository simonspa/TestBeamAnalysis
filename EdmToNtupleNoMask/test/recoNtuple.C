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

void recoNtuple(){

  // Open the raw data file:
  unsigned int    events(-9999);
  unsigned int    run(-9999);
  unsigned int    tdcPhase(9999);

  std::vector<int> *strip_S0 = 0;
  std::vector<int> *strip_S1 = 0;
  std::vector<int> *strip_S2 = 0;
  std::vector<int> *strip_S3 = 0;
  std::vector<int> *channel_S0 = 0;
  std::vector<int> *channel_S1 = 0;
  std::vector<int> *channel_S2 = 0;
  std::vector<int> *channel_S3 = 0;

  TFile *dutFile = TFile::Open("RAW_USC.00000478.0001.A.storageManager.00.0000.root");
  TTree *duttree = (TTree *) dutFile->Get("tree");
  duttree->SetBranchAddress("event", &events);
  duttree->SetBranchAddress("run", &run);
  duttree->SetBranchAddress("dut0_row", &strip_S0);
  duttree->SetBranchAddress("dut1_row", &strip_S1);
  duttree->SetBranchAddress("dut2_row", &strip_S2);
  duttree->SetBranchAddress("dut3_row", &strip_S3);
  duttree->SetBranchAddress("dut0_channel", &channel_S0);
  duttree->SetBranchAddress("dut1_channel", &channel_S1);
  duttree->SetBranchAddress("dut2_channel", &channel_S2);
  duttree->SetBranchAddress("dut3_channel", &channel_S3);
  duttree->SetBranchAddress("tdcPhase", &tdcPhase);
	
  int dutEntries = duttree->GetEntries();

  TH1F *h_row_S0 = new TH1F("row_S0","row_S0", 256, 0. , 256.);
  TH1F *h_row_S1 = new TH1F("row_S1","row_S1", 256, 0. , 256.);
  TH1F *h_row_S2 = new TH1F("row_S2","row_S2", 256, 0. , 256.);
  TH1F *h_row_S3 = new TH1F("row_S3","row_S3", 256, 0. , 256.);
  TH1F *h_row_S2_masked = new TH1F("row_S2_masked","row_S2_masked", 256, 0. , 256.);
  TH1F *h_row_S3_masked = new TH1F("row_S3_masked","row_S3_masked", 256, 0. , 256.);

  TH1F *h_clusters_S2 = new TH1F("clusters_S2","clusters_S2", 256, 0. , 256.);
  TH1F *h_clusters_S3 = new TH1F("clusters_S3","clusters_S3", 256, 0. , 256.);

  //creat recofile
  std::vector<double> stubs_M1;
  std::vector<double> stubs_M2;
  std::vector<int> *row_S2_masked = 0;
  std::vector<int> *row_S3_masked = 0;

  TFile *recofile = new TFile("RECO_USC.00000478.0001.A.storageManager.00.0000.root","RECREATE");
  TTree *recotree = new TTree("recotree","reconstructed");
  recotree->Branch("events" , &events );
  recotree->Branch("run" , &run );
  recotree->Branch("row_S0", &strip_S0);
  recotree->Branch("row_S1", &strip_S1);
  recotree->Branch("row_S2", &strip_S2);
  recotree->Branch("row_S3", &strip_S3);
  recotree->Branch("row_S2_masked", &row_S2_masked);
  recotree->Branch("row_S3_masked", &row_S3_masked);
  recotree->Branch("channel_S0", &channel_S0);
  recotree->Branch("channel_S1", &channel_S1);
  recotree->Branch("channel_S2", &channel_S2);
  recotree->Branch("channel_S3", &channel_S3);
  recotree->Branch("clusters_S0", &clusters_position_S0);
  recotree->Branch("clusters_S1", &clusters_position_S1);
  recotree->Branch("clusters_S2", &clusters_position_S2);
  recotree->Branch("clusters_S3", &clusters_position_S3);
  recotree->Branch("clusterwidth_S0", &clusterwidth_S0);
  recotree->Branch("clusterwidth_S1", &clusterwidth_S1);
  recotree->Branch("clusterwidth_S2", &clusterwidth_S2);
  recotree->Branch("clusterwidth_S3", &clusterwidth_S3);
  recotree->Branch("stubs_module_1",&stubs_module_1);
  recotree->Branch("stubs_module_2",&stubs_module_2);
  recotree->Branch("tdc", &tdcPhase);

  //events loop
  for (int i = 0; i < dutEntries; i++) {

    duttree->GetEntry(i);

    // Masking noisy strips
    // If a strip is "noisy" can be read from the hit map histograms
    // and results in a hardcoded cut on strip numbers.

    // Store all strips which are fine in vectors for further analyses:
    std::vector<int> S0_strips_masked;
    std::vector<int> S1_strips_masked;
    std::vector<int> S2_strips_masked;
    std::vector<int> S3_strips_masked;
	
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
    for(int k = 0; k < S0_strips_masked.size(); ++k){
      // ...
    }

    if(cluster_firstStrip_S0.size() != cluster_lastStrip_S0.size()){
      std::cerr << "Error in vector size cluster/S0 for first/last strip.\n";
      exit(583);
    }

    //Calculate the cluster width and position for all sensors:
    std::vector<int> clusterwidth_S0;
    std::vector<double> clusterposition_S0;

    std::vector<int> clusterwidth_S1;
    std::vector<double> clusterposition_S1;

    std::vector<int> clusterwidth_S2;
    std::vector<double> clusterposition_S2;

    std::vector<int> clusterwidth_S3;
    std::vector<double> clusterposition_S3;

    for(size_t strip = 0; strip < cluster_firstStrip_S0.size(); ++strip) {
      // ...
    }

		
    // Produce the stubs for the two modules
    // Stubs are defined between sensros S0 and S1 (module 1) and S2/S3 (module2).

		
    //stubs module 1
    if ( v_firstStripOfCluster.size() == 1 &&  v_firstStripOfCluster1.size() == 1){
      for(size_t iStrip = 0; iStrip < v_firstStripOfCluster.size() && 	v_firstStripOfCluster1.size(); ++iStrip){
	double pos_S0 = (double(v_firstStripOfCluster.at(0) + v_lastStripOfCluster.at(0))/2.);
	double pos_S1 = (double(v_firstStripOfCluster1.at(0) + v_lastStripOfCluster1.at(0))/2.);
	stubs_module_1.push_back(pos_S1 - pos_S0);
      }
    }
    //stubs module 2
    if ( v_firstStripOfCluster2.size() == 1 && v_firstStripOfCluster3.size() == 1){
      for(size_t iStrip = 0; iStrip < v_firstStripOfCluster2.size() && v_firstStripOfCluster3.size(); ++iStrip){
	double pos_S2 = (double(v_firstStripOfCluster2.at(0) + v_lastStripOfCluster2.at(0))/2.);
	double pos_S3 = (double(v_firstStripOfCluster3.at(0) + v_lastStripOfCluster3.at(0))/2.);
	stubs_module_2.push_back(pos_S3 - pos_S2);
      }
    }

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

    row_S2_masked->clear();
    row_S3_masked->clear();

	
  } //end of events loop


  recofile->cd();

  h_row_S0->Write();
  h_row_S1->Write();
  h_row_S2->Write();
  h_row_S3->Write();
  h_row_S2_masked->Write();
  h_row_S3_masked->Write();
  h_clusters_S2->Write();
  h_clusters_S3->Write();
  recotree->Write();

  recofile->Close();

  //h_position->Draw();
}

