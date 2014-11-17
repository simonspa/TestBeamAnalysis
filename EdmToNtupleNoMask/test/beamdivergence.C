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
#include "TStyle.h"
#include "stdio.h"

void beamdivergence(){

  // Book histogram:
  TH1F *h_div = new TH1F("","; #mum; Events", 100, -2000. , 2000.);

  // Prepare input data:
  std::vector<double> *clusters_position_S0= 0;
  std::vector<double> *clusters_position_S1= 0;
  std::vector<double> *clusters_position_S2= 0;
  std::vector<double> *clusters_position_S3= 0;

  // Use a zero-degree run for this analysis:
  TFile *recofile =  TFile::Open("MyRecoTTreeFile.root");
  TTree *recotree = (TTree *) recofile->Get("recotree");
  recotree->SetBranchAddress("clusters_S0", &clusters_position_S0);
  recotree->SetBranchAddress("clusters_S1", &clusters_position_S1);
  recotree->SetBranchAddress("clusters_S2", &clusters_position_S2);
  recotree->SetBranchAddress("clusters_S3", &clusters_position_S3);

  int dutEntries = recotree->GetEntries();

  // Run over all events in the tree:
  for (int i =0; i < dutEntries; i++){
    recotree->GetEntry(i);

    unsigned int cluster_size_S0 = clusters_position_S0->size();
    unsigned int cluster_size_S1 = clusters_position_S1->size();
    unsigned int cluster_size_S2 = clusters_position_S2->size();
    unsigned int cluster_size_S3 = clusters_position_S3->size();

    // Limit to event with only one cluster per sensor
    if(cluster_size_S0 == 1 && cluster_size_S1 == 1 && cluster_size_S2 == 1 && cluster_size_S3 == 1) {

      // By applying a cut to the stubs we get a cleaner signal:
      if(TMath::Abs(clusters_position_S0->front()-clusters_position_S1->front()) < 4. && TMath::Abs(clusters_position_S2->front()-clusters_position_S3->front()) < 4.) {
	// Fill histogram:
	// Average cluster positions on both M1 and M2
	// Subtract averaged cluster position M2 from M1
	// What do we have to do about the different pitches of the Modules?
	h_div->Fill(1);
      }
    }
  }

  gStyle->SetOptFit(1011);
  h_div->Draw();

}
