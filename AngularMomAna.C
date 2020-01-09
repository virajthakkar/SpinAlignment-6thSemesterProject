#define AngularMomAna_cxx
#include "AngularMomAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AngularMomAna::Loop()
{
  //Analysis Variables
  Float_t En = 200/2.0;  //in GeV units
  Float_t Mn = 0.938;  //in GeV units  //assuming Mp = Mn = 0.938 GeV
  Float_t Pn = TMath::Sqrt(En*En - Mn*Mn);  //in Gev Units //Momentum of colliding nucleons(along z axis)
  TFile *fout=new TFile("Angularmomvsb.root","recreate");
  TProfile *Lxb=new TProfile("Lxb","Lx vs b",60,0,15);
  TProfile *Lyb=new TProfile("Lyb","Ly vs b",60,0,15);
  TProfile *Lb=new TProfile("Lb","L vs b",60,0,15);
 TProfile *modLxb=new TProfile("|Lx|","|Lx| vs b",60,0,15);
 TProfile *modLyb=new TProfile("|Ly|","|Ly| vs b",60,0,15);
  //   In a ROOT session, you can do:
  //      Root > .L AngularMomAna.C
  //      Root > AngularMomAna t
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
ofstream fang;

fang.open("ang mom gold distribution root.txt");
fang.precision(10);
//file.setf(ios::scientific);
fang.setf(ios::showpoint);

ofstream avsb;

avsb.open("L vs b Au.txt");
avsb.precision(10);
//file.setf(ios::scientific);
avsb.setf(ios::showpoint);

  //Event Loop
  for (Long64_t jentry=0; jentry<nentries; jentry++)
    {

      //if(jentry > 10) break;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

   Float_t Lx = 0.;
    Float_t Ly = 0.;

    for(int i=0; i<394; i++)
	{
	  if(Stat[i] == 0) continue;
        if(i < 197)
	    {
	      Lx =  Lx + (Ny[i] * Pn)* (1.0/0.1975);
	      Ly =  Ly + ((-1.0)*Nx[i] * Pn) *(1.0/0.1975);
	    }
	  else
	    {
	      Lx =  Lx + (Ny[i] * (-1.0)*Pn) * (1.0/0.1975);
	      Ly =  Ly + ((-1.0)*Nx[i] * (-1.0)*Pn)*(1.0/0.1975);
	    }
	  //cout<<"Nxy: "<<Nx[i]<<", "<<Ny[i]<<endl;

	}

/*
Float_t L=0;
	for(int i=0; i<394; i++)
	{
	  if(Stat[i] == 0) continue;
        if(i < 197)
	    {
	      if(Nx[i]>0)
          {
              L =L+ sqrt(Nx[i]*Nx[i]+ Ny[i]*Ny[i])*(-Pn)/2;
          }
          else
            {
                L =L+ sqrt(Nx[i]*Nx[i]+ Ny[i]*Ny[i])*(Pn/2);
            }
        }
	  else
	    {
	     // Lx =  Lx + (Ny[i] * (-1.0)*Pn); // * (1.0/0.1975);
	    //  Ly =  Ly + ((-1.0)*Nx[i] * (-1.0)*Pn); // *(1.0/0.1975);
	    if(Nx[i]>0)
        {
            L =L + sqrt(Nx[i]*Nx[i]+ Ny[i]*Ny[i])*(Pn)/2;
        }
        else
        {
            L =L+ sqrt(Nx[i]*Nx[i]+ Ny[i]*Ny[i])*(-Pn/2);
        }
	    }
	  //cout<<"Nxy: "<<Nx[i]<<", "<<Ny[i]<<endl;

	}
*/

      /*
      //cout<<Pn;
      for(int i=0; i<197; i++)
	{
	  // cout<<"x"<<Nx[i]<<"y"<<Ny[i]<<endl;
	  if(Stat[i]==0) continue;

	  Lx +=  Ny[i] * Pn;
	  Ly += -Nx[i] * Pn;
	}

      double Lx2=0;
      double Ly2=0;

      for(int i=197; i<394; i++)
	{
	  if(Stat[i]==0) continue;
	  Lx2 += Ny[i]  * (-1.0*Pn);
	  Ly2 += -Nx[i] * (-1.0*Pn);
	}
      */

      //cout<<"#Event: "<<jentry<<"\tImp: "<<Imp<<" Lx  "<<Lx+Lx2<<"  Ly  "<<Ly+Ly2<<endl;
     double L= sqrt(Lx*Lx + Ly*Ly);
      Lb->Fill(Imp,abs(L));
     Lxb->Fill(Imp,Lx);
       modLxb->Fill(Imp,abs(Lx));
      Lyb->Fill(Imp,Ly);
      modLyb->Fill(Imp,abs(Ly));
      if(jentry % 1000==0)
	cout<<"#Event: "<<jentry<<"\tImp: "<<Imp<<" Lx  "<<Lx<<"  Ly  "<<Ly<<" L "<<L<<endl;
	fang<<"#Event: "<<jentry<<"\tImp: "<<Imp<<" Lx  "<<Lx<<"  Ly  "<<Ly<<" L "<<L<<endl;
avsb<<setw(25)<<Imp<<setw(25)<<L<<endl;
//cout<<"#Event: "<<jentry<<"\tImp: "<<Imp<<" L "<<L<<endl;
//fang<<"#Event: "<<jentry<<"\tImp: "<<Imp<<" L "<<L<<endl;
    }//end of event loop
  fout->Write();

}//end of program
