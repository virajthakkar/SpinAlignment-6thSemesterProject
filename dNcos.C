
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<math.h>


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//PID K+= 321
// PID Pi-=-211

void dNcos()
{
  ifstream fin;
  fin.open("PhojetPPCollision7TeV.txt",ios::in);
  if(!fin)
    {
      cout<<"Can not open the file  "<<endl;
      return -1;
    }
  Char_t eve1[6],eve2[6];
  Int_t nevt; //Event number
  Int_t nptls,process;//# particles in an event
  Int_t Pid; // PID of the particle
  Int_t slno; //Serial # of the particle
  Int_t prisec; //Primary (1)or secondary(2)
  Int_t charge; // charge of the particle
  Int_t totEvt; //Total no of events in the file
  Float_t Px,Py,Pz,Energy;  //four momentum of the particle

  //start
  fin>>eve1>>eve2>>totEvt;
  cout<<"Total events"<<totEvt<<endl;
  //Define 2D Histograms
    TH2D *KpPim = new TH2D("K+Pi-", "InvariantMass and CosTheta for K+Pi- pair", 90, 0.6, 1.5, 10, -1, 1);
    TH2D *KmPip = new TH2D("K-Pi+", "InvariantMass and CosTheta for K-Pi+ pair", 90, 0.6, 1.5, 10, -1, 1);
    TH2D *KpPip = new TH2D("K+Pi+", "InvariantMass and CosTheta for K+Pi+ pair", 90, 0.6, 1.5, 10, -1, 1);
    TH2D *KmPim = new TH2D("K-Pi-", "InvariantMass and CosTheta for K-Pi- pair", 90, 0.6, 1.5, 10, -1, 1);
        double PxKp[400];  //K+ four vector
        double PyKp[400];
        double PzKp[400];
        double EKp[400];

        double PxKm[400];  //K- four vector
        double PyKm[400];
        double PzKm[400];
        double EKm[400];


        double PxPim[400];  //Pi- four vector
        double PyPim[400];
        double PzPim[400];
        double EPim[400];

        double PxPip[400];  //Pi+ four vector
        double PyPip[400];
        double PzPip[400];
        double EPip[400];

//Definitions of vectors,variables for double "for loop" combinations
TLorentzVector vKp;
TLorentzVector vPim;
TLorentzVector vMother;
TVector3 b; //Boost
TVector3 dau; //daughter 3 momentum vector
double coss;
double M;

TLorentzVector vKm;
TLorentzVector vPip;
TLorentzVector vMother2;
TVector3 b2; //Boost
TVector3 dau2; //daughter 3 momentum vector
double coss2;
double M2;

TLorentzVector vMother3;
TVector3 b3; //Boost
TVector3 dau3; //daughter 3 momentum vector
double coss3;
double M3;
TLorentzVector vMother4;
TVector3 b4; //Boost
TVector3 dau4; //daughter 3 momentum vector
double coss4;
double M4;

//Define z axis ; z= P cross (0,0,1) direction of pp beam
TVector3 z;
//***********************************************************EVENT LOOP****************************************************
  for(Int_t ievt = 0; ievt < totEvt; ievt++)
    {
      fin>>nevt>>nptls>>process;
      cout<<"Processing Event #  =  "<<nevt<<" no. of Particles = "<<nptls<<endl;
       int c321=0;
       int c211m=0;
       int c321m=0;
       int c211=0;


      for(Int_t ipart =0;ipart < nptls; ipart++)
        {
          fin>>slno>>prisec>>Pid>>charge>>Px>>Py>>Pz>>Energy;

if(TMath::Sqrt(Px*Px + Py*Py )<0.001)
    continue;

          //selecting only primary particle
         // if(prisec != 1)continue;
	  if(Pid == 321) //K+
        {c321++;
       // cout<<setw(5)<<Pid<<setw(5)<<charge<<setw(15)<<Px<<setw(15)<<Py<<setw(15)<<Pz<<setw(15)<<Energy<<endl;
    PxKp[c321]=Px; PyKp[c321]=Py;PzKp[c321]=Pz;EKp[c321]=Energy;
//cout<<c321<<setw(5)<<Pid<<setw(5)<<charge<<setw(15)<<PxKp[c321]<<setw(15)<<PyKp[c321]<<setw(15)<<PzKp[c321]<<setw(15)<<EKp[c321]<<endl;
        }
	 else if(Pid ==-211) //Pi-
	  {c211m++;
// cout<<setw(5)<<Pid<<setw(5)<<charge<<setw(15)<<Px<<setw(15)<<Py<<setw(15)<<Pz<<setw(15)<<Energy<<endl;
	 PxPim[c211m]=Px; PyPim[c211m]=Py;PzPim[c211m]=Pz;EPim[c211m]=Energy;
	//cout<<c211m<<setw(5)<<Pid<<setw(5)<<charge<<setw(15)<<PxPim[c211m]<<setw(15)<<PyPim[c211m]<<setw(15)<<PzPim[c211m]<<setw(15)<<EPim[c211m]<<endl;
        }
   else if(Pid==-321)
    {c321m++;
PxKm[c321m]=Px; PyKm[c321m]=Py;PzKm[c321m]=Pz;EKm[c321m]=Energy;
}
    else if(Pid==211)
    {
        c211++;
PxPip[c211]=Px; PyPip[c211]=Py;PzPip[c211]=Pz;EPip[c211]=Energy;
}








        }//*****************************************************END OF PARTICLE LOOP*****************************************************************

//cout<<" c321   "<<c321<<"  c211m   "<<c211m<<" c321m "<<c321<<" c211 "<<c211<<endl;



//**********************               K+    and Pi -          *****************************

for(int i=1;i<=c321;i++) //Combinations of K+  and Pi-
{
    vKp.SetPxPyPzE(PxKp[i],PyKp[i],PzKp[i],EKp[i]);
     for(int j=1;j<=c211m;j++)
     {

vPim.SetPxPyPzE(PxPim[j],PyPim[j],PzPim[j],EPim[j]);

vMother=vKp+vPim;

//cout<<vKp[0]<<" "<<vKp[1]<<" "<<vKp[2]<<" "<<vKp[3]<<" "<<vKp.M()<<" "<<vPim[0]<<" "<<vPim[1]<<" "<<vPim[2]<<" "<<vPim[3]<<" "<<vPim.M()<<" "<<vMother[0]<<" "<<vMother[1]<<" "<<vMother[2]<<" "<<vMother[3]<<" "<<vMother.M()<<endl;

//cout<<" K+ "<<vKp.M2()<<" Pi- "<<vPim.M2()<<endl;

b.SetXYZ(-vMother.Px()/vMother.E(),-vMother.Py()/vMother.E(),-vMother.Pz()/vMother.E()); //boost



z.SetXYZ(vMother.Py()/sqrt( pow(vMother.Px(),2)+pow(vMother.Py(),2)  ),-vMother.Px()/sqrt( pow(vMother.Px(),2)+pow(vMother.Py(),2)  ),0); //z axis set

//vMother.Boost(b);
//vKp.Boost(b);
vPim.Boost(b);
//cout<<vKp[0]<<" "<<vKp[1]<<" "<<vKp[2]<<" "<<vKp[3]<<" "<<vKp.M()<<" "<<vPim[0]<<" "<<vPim[1]<<" "<<vPim[2]<<" "<<vPim[3]<<" "<<vPim.M()<<" "<<vMother[0]<<" "<<vMother[1]<<" "<<vMother[2]<<" "<<vMother[3]<<" "<<vMother.M()<<endl;


dau.SetXYZ(vPim.Px(),vPim.Py(),vPim.Pz());
//Calculating cos in Vector Meson Rest frame cos= V(daughter). z

coss=dau*z/dau.Mag();

 M=vMother.M();
//cout<<"cos "<<coss<<" Inv M " << vMother.M()<<endl;


KpPim -> Fill(M,coss);


}
}
//*******************************END of K+ and Pi-  ****************************************

//***************************K-   and    Pi+***************************************************
//

for(int i=1;i<=c321m;i++) //comb of K- and Pi+
{

vKm.SetPxPyPzE(PxKm[i],PyKm[i],PzKm[i],EKm[i]);

     for(int j=1;j<=c211;j++)
     {



vPip.SetPxPyPzE(PxPip[j],PyPip[j],PzPip[j],EPip[j]);


vMother2=vKm+vPip;

//cout<<vKp[0]<<" "<<vKp[1]<<" "<<vKp[2]<<" "<<vKp[3]<<" "<<vKp.M()<<" "<<vPim[0]<<" "<<vPim[1]<<" "<<vPim[2]<<" "<<vPim[3]<<" "<<vPim.M()<<" "<<vMother[0]<<" "<<vMother[1]<<" "<<vMother[2]<<" "<<vMother[3]<<" "<<vMother.M()<<endl;

//cout<<" K+ "<<vKp.M2()<<" Pi- "<<vPim.M2()<<endl;

b2.SetXYZ(-vMother2.Px()/vMother2.E(),-vMother2.Py()/vMother2.E(),-vMother2.Pz()/vMother2.E()); //Boost Vector



z.SetXYZ(vMother2.Py()/sqrt( pow(vMother2.Px(),2)+pow(vMother2.Py(),2)  ),-vMother2.Px()/sqrt( pow(vMother2.Px(),2)+pow(vMother2.Py(),2)  ),0);

//vMother.Boost(b);
//vKm.Boost(b2);
vPip.Boost(b2);
//cout<<vKp[0]<<" "<<vKp[1]<<" "<<vKp[2]<<" "<<vKp[3]<<" "<<vKp.M()<<" "<<vPim[0]<<" "<<vPim[1]<<" "<<vPim[2]<<" "<<vPim[3]<<" "<<vPim.M()<<" "<<vMother[0]<<" "<<vMother[1]<<" "<<vMother[2]<<" "<<vMother[3]<<" "<<vMother.M()<<endl;


dau2.SetXYZ(vPip.Px(),vPip.Py(),vPip.Pz()); //daughter
//Calculating cos in Vector Meson Rest frame cos= V(daughter). z

coss2=dau2*z/dau2.Mag();

M2=vMother2.M();
//cout<<"cos "<<coss2<<" Inv M " << M2<<endl;


KmPip -> Fill(M2,coss2);
}
}
//********************************END of K- and Pi +********************************************


//**********************                  K+    and Pi +          *****************************
for(int i=1;i<=c321;i++) //Combinations of K+  and Pi+
{
    vKp.SetPxPyPzE(PxKp[i],PyKp[i],PzKp[i],EKp[i]);
     for(int j=1;j<=c211;j++)
     {

vPip.SetPxPyPzE(PxPip[j],PyPip[j],PzPip[j],EPip[j]);

vMother3=vKp+vPip;

//cout<<vKp[0]<<" "<<vKp[1]<<" "<<vKp[2]<<" "<<vKp[3]<<" "<<vKp.M()<<" "<<vPim[0]<<" "<<vPim[1]<<" "<<vPim[2]<<" "<<vPim[3]<<" "<<vPim.M()<<" "<<vMother[0]<<" "<<vMother[1]<<" "<<vMother[2]<<" "<<vMother[3]<<" "<<vMother.M()<<endl;

//cout<<" K+ "<<vKp.M2()<<" Pi- "<<vPim.M2()<<endl;

b3.SetXYZ(-vMother3.Px()/vMother3.E(),-vMother3.Py()/vMother3.E(),-vMother3.Pz()/vMother3.E()); //boost



z.SetXYZ(vMother3.Py()/sqrt( pow(vMother3.Px(),2)+pow(vMother3.Py(),2)  ),-vMother3.Px()/sqrt( pow(vMother3.Px(),2)+pow(vMother3.Py(),2)  ),0); //z axis set

//vMother.Boost(b);
//vKp.Boost(b);
vPip.Boost(b3);
//cout<<vKp[0]<<" "<<vKp[1]<<" "<<vKp[2]<<" "<<vKp[3]<<" "<<vKp.M()<<" "<<vPim[0]<<" "<<vPim[1]<<" "<<vPim[2]<<" "<<vPim[3]<<" "<<vPim.M()<<" "<<vMother[0]<<" "<<vMother[1]<<" "<<vMother[2]<<" "<<vMother[3]<<" "<<vMother.M()<<endl;


dau3.SetXYZ(vPip.Px(),vPip.Py(),vPip.Pz());
//Calculating cos in Vector Meson Rest frame cos= V(daughter). z

coss3=dau3*z/dau3.Mag();

 M3=vMother3.M();
//cout<<"cos "<<coss3<<" Inv M " << vMother3.M()<<endl;


KpPip -> Fill(M3,coss3);


}
}
//*******************************END of K+ and Pi+  ****************************************


//**********************               K-   and Pi -          *****************************

for(int i=1;i<=c321m;i++) //Combinations of K- and Pi-
{
    vKm.SetPxPyPzE(PxKm[i],PyKm[i],PzKm[i],EKm[i]);
     for(int j=1;j<=c211m;j++)
     {

vPim.SetPxPyPzE(PxPim[j],PyPim[j],PzPim[j],EPim[j]);

vMother4=vKm+vPim;

//cout<<vKp[0]<<" "<<vKp[1]<<" "<<vKp[2]<<" "<<vKp[3]<<" "<<vKp.M()<<" "<<vPim[0]<<" "<<vPim[1]<<" "<<vPim[2]<<" "<<vPim[3]<<" "<<vPim.M()<<" "<<vMother[0]<<" "<<vMother[1]<<" "<<vMother[2]<<" "<<vMother[3]<<" "<<vMother.M()<<endl;

//cout<<" K+ "<<vKp.M2()<<" Pi- "<<vPim.M2()<<endl;

b4.SetXYZ(-vMother4.Px()/vMother4.E(),-vMother4.Py()/vMother4.E(),-vMother4.Pz()/vMother4.E()); //boost



z.SetXYZ(vMother4.Py()/sqrt( pow(vMother4.Px(),2)+pow(vMother4.Py(),2)  ),-vMother4.Px()/sqrt( pow(vMother4.Px(),2)+pow(vMother4.Py(),2)  ),0); //z axis set

//vMother.Boost(b);
//vKp.Boost(b);
vPim.Boost(b4);
//cout<<vKp[0]<<" "<<vKp[1]<<" "<<vKp[2]<<" "<<vKp[3]<<" "<<vKp.M()<<" "<<vPim[0]<<" "<<vPim[1]<<" "<<vPim[2]<<" "<<vPim[3]<<" "<<vPim.M()<<" "<<vMother[0]<<" "<<vMother[1]<<" "<<vMother[2]<<" "<<vMother[3]<<" "<<vMother.M()<<endl;


dau4.SetXYZ(vPim.Px(),vPim.Py(),vPim.Pz());
//Calculating cos in Vector Meson Rest frame cos= V(daughter). z

coss4=dau4*z/dau4.Mag();

 M4=vMother4.M();
//cout<<"cos "<<coss4<<" Inv M " << vMother4.M()<<endl;


KmPim -> Fill(M4,coss4);


}
}
//*******************************END of K- and Pi-  ****************************************







}//***************************end of event loop**********************************


//Fill Histograms to File
TFile *File = new TFile("Histograms.root", "recreate");
File->cd();
KpPim->Write();
KmPip->Write();
KpPip->Write();
KmPim->Write();


//double Nc=KpPim->GetBinContent(1,1);
//cout<<Nc;






  return 0;
}
