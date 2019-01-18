#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<math.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

void read()
{
  gStyle->SetLineWidth(2);
  gStyle->SetOptFit(1);

  TFile *f = TFile::Open("Histograms.root"); //extract 2D scatter plots from the root file
  f->ls{};
  TH2D *hKpPim=(TH2D*)f->Get("K+Pi-");
  TH2D *hKmPip=(TH2D*)f->Get("K-Pi+");
  TH2D *hKpPip=(TH2D*)f->Get("K+Pi+");
  TH2D *hKmPim=(TH2D*)f->Get("K-Pi-");

  TF1 *fun[10];
  TF1 *funbkg[10];
  TF1 *funSig[10];
  //hKpPim->Draw();
  TH1D *hsig1[10]; //Histogram for signal 1
  TH1D *hsig2[10]; //Histogram for signal 2
  TH1D *hbkg1[10]; //Histogram for background 1
  TH1D *hbkg2[10]; //Histogram for background 2

  TH1D *hsigt[10];  //Histogram for total signal
  TH1D *hbkgt[10];  //Histogram for total background
  TH1D *hsigonly[10]; //Histogram of total signal minus background
  TCanvas *can[10];
  TCanvas *canbkg[10]; //Canvas for background

  Double_t Mass[10],Width[10],Yield[10];
  Double_t ErMass[10],ErWidth[10],ErYield[10];
  for(int i=0;i<10;i++)
    {
      can[i] = new TCanvas(Form("can%d",i),"",10,10,600,600);
      can[i]->SetLeftMargin(0.2);
      can[i]->SetRightMargin(0.05);
      can[i]->SetTopMargin(0.05);
      can[i]->SetBottomMargin(0.13);

      canbkg[i] = new TCanvas(Form("canbkg%d",i),"",10,10,600,600);
      canbkg[i]->SetLeftMargin(0.2);
      canbkg[i]->SetRightMargin(0.05);
      canbkg[i]->SetTopMargin(0.05);
      canbkg[i]->SetBottomMargin(0.13);


      hsigt[i]=new TH1D(Form("Signal unlike pair%d",i),"",90,0.6,1.5);
      hbkgt[i]=new TH1D(Form("bkgt%d",i),"",90,0.6,1.5);
      hsigonly[i]=new TH1D(Form("Final Signal %d",i),"",90,0.6,1.5);
      hsig1[i] = hKpPim->ProjectionX(Form("KpPim%d",i),i+1,i+1,"e");//Projection for each cos theta bin for K+Pi-
      hsig2[i] = hKmPip->ProjectionX(Form("KmPip%d",i),i+1,i+1,"e");//Projection for each cos theta bin for K-Pi+
      hbkg1[i] = hKpPip->ProjectionX(Form("KpPip%d",i),i+1,i+1,"e");//Projection for each cos theta bin for K+Pi+
      hbkg2[i] = hKmPim->ProjectionX(Form("KmPim%d",i),i+1,i+1,"e");//Projection for each cos theta bin for K-Pi-
      hsigt[i]->Add(hsig1[i],hsig2[i],1,1); //Total signal=signal1+signal2

      for(int bin=1;bin<hbkgt[i]->GetNbinsX();bin++)
	{ //Background combined together by like sign technique
	  double totbkg=2*TMath::Sqrt(hbkg1[i]->GetBinContent(bin)*hbkg2[i]->GetBinContent(bin));
	  double bkgerr;
	  bkgerr=0.01;
	  hbkgt[i]->SetBinContent(bin,totbkg);
	  hbkgt[i]->SetBinError(bin,bkgerr);
	}
      canbkg[i]->cd();
      hsigt[i]->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
      hsigt[i]->GetXaxis()->SetTitleOffset(1.2);
      hsigt[i]->GetYaxis()->SetTitleOffset(1.7);
      hsigt[i]->GetYaxis()->SetTitle("Counts");
      hsigt[i]->SetMarkerColor(1);
      hsigt[i]->SetMarkerStyle(20);
      hsigt[i]->SetMarkerSize(1.0);
      hsigt[i]->Draw();        //Signal Total=Signal1(K+Pi-) + Signal2(K-Pi+)

      hbkgt[i]->SetMarkerColor(2);
      hbkgt[i]->SetMarkerStyle(24);
      hbkgt[i]->SetMarkerSize(1.0);
      hbkgt[i]->Draw("same");


      canbkg[i]->SaveAs(Form("ScaledBkgInvariantMassDistInCosThetaStarBin%d.gif",i));
      canbkg[i]->SaveAs(Form("ScaledBkgInvariantMassDistInCosThetaStarBin%d.png",i));
      canbkg[i]->Print(Form("ScaledBkgInvariantMassDistInCosThetaStarBin%d.eps",i));
      canbkg[i]->cd();

      fun[i]= new TF1(Form("BWfunc%d",i),fitfunction,0.65,1.2,6); //Fit function to signalonly:BW+poly2

      fun[i]->SetParameters(1000,0.048,0.896,150,150,60);
      fun[i]->SetParName(0,"A");
      fun[i]->SetParName(1,"#Gamma");
      fun[i]->SetParName(2,"Mass");
      fun[i]->SetParName(3,"a_{0}"); //coefficient of degree2 in poly2
      fun[i]->SetParName(4,"a_{1}"); //coefficient of degree1 in poly2
      fun[i]->SetParName(5,"a_{2}"); //coefficient of degree0 in poly2
      //fun[i]->SetParLimits(0,0,10000000000);
      //fun[i]->SetParLimits(2,0.8,0.95);
      //Final Signal=Total Signal(due to K+Pi- & K-Pi+) - Total background (like sign K+Pi+ and K-Pi-)
      hsigonly[i]->Add(hsigt[i],hbkgt[i],1,-1);
      can[i]->cd();
      hsigonly[i]->GetXaxis()->SetTitle("M_{K#pi} (GeV/c^{2})");
      hsigonly[i]->GetXaxis()->SetTitleOffset(1.2);
      hsigonly[i]->GetYaxis()->SetTitleOffset(1.7);
      hsigonly[i]->GetYaxis()->SetTitle("Counts");
      hsigonly[i]->SetMarkerColor(1);
      hsigonly[i]->SetMarkerStyle(20);
      hsigonly[i]->SetMarkerSize(1.0);
      hsigonly[i]->Draw();
     // gPad->BuildLegend();
      hsigonly[i]->Fit(fun[i],"REM");
      funbkg[i]= new TF1(Form("Bkgfunc%d",i),"([0]*x*x+[1]*x+[2])",0.65,1.2);
      funbkg[i]->SetParameters(fun[i]->GetParameter(3),fun[i]->GetParameter(4),fun[i]->GetParameter(5));
      funbkg[i]->SetLineColor(3);
      funbkg[i]->Draw("same");
      funSig[i]= new TF1(Form("Sigfunc%d",i),"(0.01*([0]/2.*TMath::Pi())*([1])*(1./( pow( x-[2],2) + [1]*[1]*0.25)))",0.65,1.2);
      funSig[i]->SetParameters(fun[i]->GetParameter(0),fun[i]->GetParameter(1),fun[i]->GetParameter(2));
      funSig[i]->SetLineColor(4); //4 Blue
      funSig[i]->Draw("same"); //BreitWigner function+Poly2
//gPad->BuildLegend();
      can[i]->SaveAs(Form("InvariantMassDistInCosThetaStarBin%d.gif",i));
      can[i]->SaveAs(Form("InvariantMassDistInCosThetaStarBin%d.png",i));
      can[i]->Print(Form("InvariantMassDistInCosThetaStarBin%d.eps",i));
      can[i]->cd();
      Mass[i] = fun[i]->GetParameter(2); //Taking Values of Parameters from the fitted function
      Width[i] = fun[i]->GetParameter(1)*1000;
      Double_t y = fun[i]->GetParameter(0); //Get Yield
      ErMass[i] = fun[i]->GetParError(2);
      ErWidth[i] = fun[i]->GetParError(1)*1000;
      Double_t ery = fun[i]->GetParError(0); //error in  yield
      Yield[i] = y / 0.2; //normalized by bin width
      ErYield[i] = ery /0.2;
    }

  TLine *lm = new TLine(-1,0.896,1,0.896);
  lm->SetLineColor(4);
  lm->SetLineWidth(2);
  lm->SetLineStyle(2);
  TLine *lw = new TLine(-1,52,1,52);
  lw->SetLineColor(4);
  lw->SetLineWidth(2);
  lw->SetLineStyle(2);

  Double_t CosThetaStar[10] = {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};
  Double_t ErCosThetaStar[10] = {0,0,0,0,0,0,0,0,0,0};
  TGraphErrors *grM = new TGraphErrors(10,CosThetaStar,Mass,ErCosThetaStar,ErMass);
  TGraphErrors *grW = new TGraphErrors(10,CosThetaStar,Width,ErCosThetaStar,ErWidth); //in MeV
  TGraphErrors *grdNdCostTheta = new TGraphErrors(10,CosThetaStar,Yield,ErCosThetaStar,ErYield);

  TCanvas *c = new TCanvas("c","",10,10,600,600);
  c->SetLeftMargin(0.2);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.13);
  grM->SetTitle("Mass vs. cos#theta");
  grM->SetMaximum(0.905);
  grM->SetMinimum(0.89);
  grM->SetMarkerColor(1);
  grM->SetMarkerStyle(20);
  grM->SetMarkerSize(1.0);
  grM->GetXaxis()->SetTitle("cos#theta");
  grM->GetXaxis()->SetTitleSize(0.05);
  grM->GetXaxis()->SetTitleOffset(1.25);
  grM->GetXaxis()->SetTitleFont(42);
  grM->GetXaxis()->CenterTitle(true);
  grM->GetXaxis()->SetLabelSize(0.05);
  grM->GetXaxis()->SetLabelFont(42);
  grM->GetXaxis()->SetNdivisions(505);
  grM->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
  grM->GetYaxis()->SetTitleFont(42);
  grM->GetYaxis()->CenterTitle(true);
  grM->GetYaxis()->SetTitleSize(0.05);
  grM->GetYaxis()->SetTitleOffset(2.);
  grM->GetYaxis()->SetLabelSize(0.05);
  grM->GetYaxis()->SetLabelFont(42);
  grM->GetYaxis()->SetNdivisions(505);
  grM->Draw("AP");
  lm->Draw();
  c->SaveAs("MassVsCosThetaStar.gif");
  c->SaveAs("MassVsCosThetaStar.png");
  c->Print("MassVsCosThetaStar.eps");
  c->cd();
  TCanvas *c1 = new TCanvas("c1","",10,10,600,600);
  c1->SetLeftMargin(0.2);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.13);
  grW->SetTitle("Width vs. cos#theta");
  grW->SetMaximum(60);
  grW->SetMinimum(0);
  grW->SetMarkerColor(1);
  grW->SetMarkerStyle(20);
  grW->SetMarkerSize(1.0);
  grW->GetXaxis()->SetTitle("cos#theta");
  grW->GetXaxis()->SetTitleSize(0.05);
  grW->GetXaxis()->SetTitleOffset(1.25);
  grW->GetXaxis()->SetTitleFont(42);
  grW->GetXaxis()->CenterTitle(true);
  grW->GetXaxis()->SetLabelSize(0.05);
  grW->GetXaxis()->SetLabelFont(42);
  grW->GetXaxis()->SetNdivisions(505);
  grW->GetYaxis()->SetTitle("#Gamma (MeV/c^{2})");
  grW->GetYaxis()->SetTitleFont(42);
  grW->GetYaxis()->CenterTitle(true);
  grW->GetYaxis()->SetTitleSize(0.05);
  grW->GetYaxis()->SetTitleOffset(1.7);
  grW->GetYaxis()->SetLabelSize(0.05);
  grW->GetYaxis()->SetLabelFont(42);
  grW->GetYaxis()->SetNdivisions(505);
  grW->Draw("AP");
  lw->Draw();
  c1->SaveAs("WidthVsCosThetaStar.gif");
  c1->SaveAs("WidthVsCosThetaStar.png");
c1->Print("WidthVsCosThetaStar.eps");
  c1->cd();

  TF1 *funRho = new TF1("funRho",RhoZeroZero,-1,1.,2);
  funRho->SetParameters(1000,0.3);
  funRho->SetParNames("Const.","#rho_{00}");
  TCanvas *c2 = new TCanvas("c2","",10,10,600,600);
  c2->SetLeftMargin(0.2);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.05);
  c2->SetBottomMargin(0.13);
  grdNdCostTheta->SetMaximum(600000);
  grdNdCostTheta->SetMinimum(540000);
  grdNdCostTheta->SetTitle("dN/dcos#theta vs. cos#theta");
  grdNdCostTheta->SetMarkerColor(1);
  grdNdCostTheta->SetMarkerStyle(20);
  grdNdCostTheta->SetMarkerSize(1.0);
  grdNdCostTheta->GetXaxis()->SetTitle("cos#theta");
  grdNdCostTheta->GetXaxis()->SetTitleSize(0.05);
  grdNdCostTheta->GetXaxis()->SetTitleOffset(1.25);
  grdNdCostTheta->GetXaxis()->SetTitleFont(42);
  grdNdCostTheta->GetXaxis()->CenterTitle(true);
  grdNdCostTheta->GetXaxis()->SetLabelSize(0.05);
  grW->GetXaxis()->SetLabelFont(42);
  grdNdCostTheta->GetXaxis()->SetNdivisions(505);
  grdNdCostTheta->GetYaxis()->SetTitle("dN/dcos#theta");
  grdNdCostTheta->GetYaxis()->SetTitleFont(42);
  grdNdCostTheta->GetYaxis()->CenterTitle(true);
  grdNdCostTheta->GetYaxis()->SetTitleSize(0.05);
  grdNdCostTheta->GetYaxis()->SetTitleOffset(1.7);
  grdNdCostTheta->GetYaxis()->SetLabelSize(0.05);
  grdNdCostTheta->GetYaxis()->SetLabelFont(42);
  grdNdCostTheta->GetYaxis()->SetNdivisions(505);
  grdNdCostTheta->Draw("AP");
  grdNdCostTheta->Fit(funRho,"REM");
  c2->SaveAs("dNdcosthetaVsCosThetaStar.gif");
  c2->SaveAs("dNdcosthetaVsCosThetaStar.png");
  c2->Print("dNdcosthetaVsCosThetaStar.eps");
  c2->cd();
}

double fitfunction(double *x,double *par)
{
    return 0.01*(par[0]/2.*TMath::Pi())*(par[1])*(1./( pow( x[0]-par[2],2) + par[1]*par[1]*0.25)) + par[3]*x[0]*x[0] +par[4]*x[0] +par[5];
}
double RhoZeroZero(double *x,double *par)
{
  return par[0]*(3./4.)*((1.-par[1])+(3*par[1]-1)*TMath::Cos(x[0]));
}
