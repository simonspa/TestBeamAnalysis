#include <fstream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"

void stubs_angle(){

  // Get a canvas:
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->SetTickx();
  c1->SetTicky();

  double angle, MV;
  TGraph *g1 = new TGraph();

  // Read in text file:
  ifstream fp;
  fp.open("rad.txt");
  while(!fp.eof()) {
    for (int i = 0 ; i < 8; ++i) {
      if(fp >> angle >> MV )
	// Fill in the data points:
	g1->SetPoint(i,angle,MV);
    }
  }
  fp.close();

  // Fit the TGraph:
  double d = 0;
  TF1 *func = new TF1("fit","my-fit-function");
  func->SetParNames ("delta0","d","ang0");
  g1->Fit("fit");

  // Do some decoration:
  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(1);
  g1->Draw("ap");
  g1->SetMinimum(0);
  g1->SetMaximum(20);
  g1->GetXaxis()->SetTitle("angle [rad]");
  g1->GetYaxis()->SetTitle("stubs mean value");
  g1->SetTitle("rotated DUT");

  // Also add a legend:
  leg = new TLegend(0.3,0.8,0.6,0.9);
  leg->AddEntry(g1,"mv","p");
  leg->AddEntry(fit,"Fit","L");
  leg->Draw();
}
