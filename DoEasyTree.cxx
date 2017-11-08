//g++ -g -O2 -std=c++0x DoEasyTree.cxx -L$SAD -lSharcAnalysis -I$SAD `grsi-config --cflags --all-libs --root` -o easy.out

// uses delcorr to apply energy loss correction to PID cuts
#include<Globals.h>
#include<TChannel.h>
#include<TTigress.h>
#include<TSharc.h>
#include<TTriFoil.h>
#include "TSharcAnalysis.h"
#include<TReaction.h>
#include<TRF.h>
#include<TSRIM.h>
#include<GValue.h>

#include<sstream>
#include<stdio.h>
#include<fstream>
#include<cmath>

#include<Rtypes.h>
#include<TSystem.h>
#include<TVector3.h>
#include<TFile.h>
#include<TChain.h>
#include<TStopwatch.h>
#include<TCutG.h>
#include<TMath.h>

struct InputOptions {
  std::string treename;
  int A;
  int nentries;
  TVector3 offset; //mm
  double TargetThick; // um
  std::string TargetMat; // um
  double EperU; // MeV
  UInt_t addmax; // size of addback container
  UInt_t piecewise; // Controls piecewise calibrations for detector 5
  double dbox_coeffs[2];
  double ubox_coeffs[2];
  double uqqq_coeffs[2];
  double pad5_coeffs1[2];
  double pad5_coeffs2[2];
  double pad_coeffs[4][2];
  double ecut;
  std::string name;
  std::string cutfile;
  std::string calfile;
  std::string logfile;
  std::string badstripsfile;
//  std::string typespec;
  std::string chaindata;
};

void InputOptionsPrint(InputOptions *runinfo){
  printf( DCYAN "Input file loaded: \n\t" RESET_COLOR "'%s' \n",runinfo->name.c_str());
	printf("\tstd::string treename   = %s\n",runinfo->treename.c_str());
	printf("\tint A                  = %i\n",runinfo->A);
	printf("\tDouble_t EperU         = %5.3f\n",runinfo->EperU);
	printf("\tint nentries           = %i\n",runinfo->nentries);
	printf("\tTVector3 offset        = %.02f\t%.02f\t%.02f\n",runinfo->offset.X(),runinfo->offset.Y(),runinfo->offset.Z());
	printf("\tDouble_t TargetThick   = %5.3f\n",runinfo->TargetThick);
	printf("\tDouble_t TargetMat     = %s\n",runinfo->TargetMat.c_str());
	printf("\tDouble_t ecut          = %5.3f\n",runinfo->ecut);		
	printf("\tUInt_t addmax          = %i\n",runinfo->addmax);
	printf("\tUInt_t piecewise       = %i\n",runinfo->piecewise);
	printf("\tdouble dbox_coeffs[]  = %.02f\t%.02f\n",runinfo->dbox_coeffs[0],runinfo->dbox_coeffs[1]);
	printf("\tdouble ubox_coeffs[]  = %.02f\t%.02f\n",runinfo->ubox_coeffs[0],runinfo->ubox_coeffs[1]);
	printf("\tdouble uqqq_coeffs[]  = %.02f\t%.02f\n",runinfo->uqqq_coeffs[0],runinfo->uqqq_coeffs[1]);
	printf("\tdouble pad5_coeffs1[]  = %.02f\t%.02f\n",runinfo->pad5_coeffs1[0],runinfo->pad5_coeffs1[1]);
	printf("\tdouble pad5_coeffs2[]  = %.02f\t%.02f\n",runinfo->pad5_coeffs2[0],runinfo->pad5_coeffs2[1]);
	printf("\tdouble pad_coeffs[0][]  = %.02f\t%.02f\n",runinfo->pad_coeffs[0][0],runinfo->pad_coeffs[0][1]);
	printf("\tdouble pad_coeffs[1][]  = %.02f\t%.02f\n",runinfo->pad_coeffs[1][0],runinfo->pad_coeffs[1][1]);
	printf("\tdouble pad_coeffs[2][]  = %.02f\t%.02f\n",runinfo->pad_coeffs[2][0],runinfo->pad_coeffs[2][1]);
	printf("\tdouble pad_coeffs[3][]  = %.02f\t%.02f\n",runinfo->pad_coeffs[3][0],runinfo->pad_coeffs[3][1]);
	printf("\tstd::string cutfile    = %s\n",runinfo->cutfile.c_str());
	printf("\tstd::string calfile    = %s\n",runinfo->calfile.c_str());
	printf("\tstd::string logfile    = %s\n",runinfo->logfile.c_str());
	printf("\tstd::string badstripsfile    = %s\n",runinfo->badstripsfile.c_str());
//	printf("\tstd::string typespec   = %s\n",runinfo->typespec.c_str());
	printf("\tstd::string chaindata  = %s\n",runinfo->chaindata.c_str());   //input analysis trees.

  return;
}

InputOptions *SetOptions(const char*);
void GainSwitch(UInt_t entry);
Bool_t InitVariables(InputOptions *);
Int_t CheckInputTree(const char*);
Int_t ProcessData(InputOptions*, const char *,const char *);

// system variables
Double_t pi = TMath::Pi();
Double_t d2r = TMath::DegToRad();
Double_t r2d = TMath::RadToDeg();
TReaction *reaction[5];
const char ion[] = {'p','p','d','c','t'};
TChannel *channel = 0;
TVector3 position;
TSRIM *srim[5];
TTigressHit *TH = 0;
Double_t thetaraw;
TVector3 rvec, tpos;
UInt_t typectot[]={0,0,0,0,0};
UInt_t typecrun[]={0,0,0,0,0};
UInt_t gamctot[]={0,0,0,0,0};
UInt_t gamcrun[]={0,0,0,0,0};
 
// cuts
TCutG *IsProtPid[4],*IsDeutPid[4],*IsTritPid[4],*IsCarbPid[4];
TCutG *IsProtElastic,*IsDeutElastic,*IsCarbElastic,*IsTritElastic;
TCutG *IsGoodCharge;
// for piecewise calibration
TCutG *IsProtPid5_1,*IsDeutPid5_1,*IsTritPid5_1,*IsCarbPid5_1;
TCutG *IsProtPid5_2,*IsDeutPid5_2,*IsTritPid5_2,*IsCarbPid5_2;

//////OUTTREE VARIABLES /////////
UShort_t type;
UShort_t det;
UShort_t front;
UShort_t back;
Double_t thetalab;
Double_t phi;
Double_t phicorr;
Double_t delcorr;

Double_t dcharge; //  front charge /25.0
Double_t dchargeb;//  back charge / 25.0
Double_t pcharge; //  /125.0 
Double_t dengalpha; // calfile applied (gain matched)
Double_t denergy;
Double_t penergy;
Double_t exc;
Double_t ekin;

Double_t thetacm;
Double_t phase; // RF phase used for sharc timing

Double_t ekinthry;
Double_t dengthry; // expected energy deposited in delta (elastics only)
Double_t pengthry; // expected energy deposited in pad (elastics only)

Double_t rekin;
Double_t rbeta;
Double_t rthetalab;

UShort_t addsize; // can't be larger than 10!
Double_t eadd[10];
Double_t tadd[10];
Double_t thetadd[10];
Double_t tigeng[10];
UShort_t tigdet[10];
UShort_t tigcry[10];
UShort_t tigseg[10];

UShort_t tsize; // can't be larger than 10!
Double_t tdop[10];
Double_t theth[10];
Double_t teng[10];
Double_t ttim[10];
UShort_t tdet[10];
UShort_t tcry[10];
UShort_t tseg[10];

Int_t    tbeam;
Bool_t   beam;

UInt_t raw_entry;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


InputOptions *runinfo;

int main(int argc, char **argv)	{ 
	
	if(argc!=4) {
	  printf( DRED "try: ./easy.out input_file analysis_tree output_tree  instead." RESET_COLOR "\n");
	  return 1;
  }
	if(CheckInputTree(argv[2])==-1){
    printf("\n\n Tree did not contain enough TGRSI Branches!!\n\n");
    return 1;
  }

 // GValue v("TIGZ",-5.0); // set tigress position offset
  runinfo = SetOptions(argv[1]); 
	InputOptionsPrint(runinfo); // duh
  if(!InitVariables(runinfo)) // makes reactions and histograms
    exit(1);

  Int_t s = ProcessData(runinfo,argv[2],argv[3]); // loop over data and fill histograms		
    //-> SetBranchVaribles();
    //-> SetupOutTree();
    //-> ClearTreeVaribles();
    //-> CloseOutTree();

	return 0;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void ClearTreeVariables(){

  type     = 10; // should be cleared to something that won't be mistaken as good data
  det      = 0;
  front    = 0;
  back     = 0;    
  thetalab = 0.00;
  phi      = 0.00;
  phicorr  = 0.00; // corrects increased detector thickness due to phi angle
  
  dcharge  = 0.00;
  dchargeb = 0.00;
  pcharge  = 0.00;
  dengalpha= 0.00;
  denergy  = 0.00;
  penergy  = 0.00;
  exc      = 0.00;
  ekin     = 0.00;
  
  rekin    = 0.0;
  rbeta    = 0.0;
  rthetalab= 0.0;

  ekinthry = 0.00;
  dengthry = 0.00;
  pengthry = 0.00;

  phase    = 0.00;
  thetacm  = 0.00;

  addsize  = 0; // can't be larger than 10!
  tsize    = 0; // can't be larger than 10!
  tbeam    = 0;
  beam     = false;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GainSwitch(){
  IsProtPid[0] = IsProtPid5_2;
  IsDeutPid[0] = IsDeutPid5_2;
  IsTritPid[0] = IsDeutPid5_2;
  IsCarbPid[0] = IsDeutPid5_2;
  
  runinfo->pad_coeffs[0][0] = runinfo->pad5_coeffs2[0];
  runinfo->pad_coeffs[0][1] = runinfo->pad5_coeffs2[1];
  printf( DRED "\n\tUsing pad 5 coefficients :   offs = %.2f\t gain = %.2f\n" RESET_COLOR "",runinfo->pad_coeffs[0][0],runinfo->pad_coeffs[0][1]);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// the charges here are actually alpha source gain matched charges, not raw pre-amplifier charges
UInt_t GetType(UInt_t entry, UInt_t d, Double_t delcorr, Double_t dchg, Double_t pchg, Double_t ekin, Double_t thetlab){

  if(d>=5 && d<=8){

     if(pchg<50.0){ // is no pad then see if it's elastic proton or deuteron
        if(IsProtElastic && IsProtElastic->IsInside(thetlab,dchg)) return 1;
        if(IsDeutElastic && IsDeutElastic->IsInside(thetlab,dchg)) return 2;
        if(IsCarbElastic && IsCarbElastic->IsInside(thetlab,dchg)) return 3;
        if(IsTritElastic && IsTritElastic->IsInside(thetlab,dchg)) return 4;
        return 10;
     }
     
     if(IsProtPid[d-5] && IsProtPid[d-5]->IsInside(pchg,dchg*delcorr)){
//	      return 0; // all proton PID hits are reconstructed as dp
       
       // kinematic energy of sr(p,p) + extra energy to prevent leaking across
       Double_t ekinpp = reaction[1]->GetTLab(thetlab*d2r,2)*1e3; 
       if(ekin>ekinpp+1500.0) 
         return 0;
       else
        return 1;
      }
     
     if(IsDeutPid[d-5] && IsDeutPid[d-5]->IsInside(pchg,dchg*delcorr))
       return 2;
     
     if(IsCarbPid[d-5] && IsCarbPid[d-5]->IsInside(pchg,dchg*delcorr))
       return 2;
     
     if(IsTritPid[d-5] && IsTritPid[d-5]->IsInside(pchg,dchg*delcorr))
       return 4;

  } else if(d>=9 && d<=16)// everything upstream is assumed to be dp
    return 0;
  
  return 10;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// manual pad calibration using input coefficients
Double_t CalibratePad(UShort_t d, Double_t padcharge){
  if(padcharge<50.0)
     return 0.0;
  return runinfo->pad_coeffs[d-5][0] + runinfo->pad_coeffs[d-5][1] * padcharge;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Double_t CalibrateDelta(UShort_t detector, Double_t delchgalpha){

  Double_t delenergy = 0.00;

  if(detector>=5 && detector<=8)
    delenergy = runinfo->dbox_coeffs[0] + runinfo->dbox_coeffs[1]*delchgalpha; 
  else if(detector>=9 && detector<=12)
    delenergy = runinfo->ubox_coeffs[0] + runinfo->ubox_coeffs[1]*delchgalpha; 
  else if(detector>=13)
    delenergy = runinfo->uqqq_coeffs[0] + runinfo->uqqq_coeffs[1]*delchgalpha; 

  return delenergy;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// including entry # as an argument makes it possible to do time dependent analysis/calibration (preamps & trifoil)
bool SetBranchVariables(UInt_t entry, TTigress *tigress, TSharc *sharc, TTriFoil *trifoil, TRF *rf) {

  ClearTreeVariables();

// don't mess around with multiplicity events
  if(sharc->GetMultiplicity()!=1)
		  return false;

  phase = rf->Phase(); // SHARC phase (time) set

  if(trifoil){
    tbeam = trifoil->TBeam();
    beam = trifoil->Beam();
  }

  TSharcHit *SH = sharc->GetSharcHit(0);

  det = SH->GetDetector();			
  front = SH->GetFrontStrip();
  back = SH->GetBackStrip();  
  
  // throw out bad strips at this level
  if(TSharcAnalysis::BadStrip(det,front,-1) || TSharcAnalysis::BadStrip(det,-1,back))
    return false;

  ////////////////////////////////////////
  dcharge = SH->GetFrontCharge();
  dchargeb = SH->GetBackCharge();
  if(!IsGoodCharge->IsInside(det,dcharge-dchargeb)) // front-back energy cut will minimize noise
     return false;

  // Get alpha calibrated (gain matched) delta energies out of analysis trees and apply low energy upstream cut
  dengalpha = SH->GetDeltaE(); // returns d_energy_front
  if(det>8 && dengalpha<runinfo->ecut)
     return false;
  
  denergy = CalibrateDelta(det,dengalpha);
  ////////////////////////////////////////
  
  thetaraw = TSharc::GetPosition(det,front,back).Theta()*r2d;
  
  position = TSharcAnalysis::GetPosition(det,front,back); // TSA has offset stored and applies it to position internally
  thetalab = TSharcAnalysis::RandomizeThetaLab(det,front,back);
  phi = position.Phi()*r2d;			
  phicorr = cos((phi + (8-det)*90.0)*d2r); 
  delcorr = phicorr*sin(thetalab*d2r); 
  ////////////////////////////////////////
  pcharge = SH->GetPadCharge();
  penergy = CalibratePad(det,pcharge);
  
  // get theta without any position offset applied, to match kinematic cuts
  type = GetType(entry,det,delcorr,dengalpha,pcharge,denergy+penergy,thetaraw);
  if(type==10) // at this point throw away everything I can't identify  
    return false;
	
  typectot[type]++;
  typecrun[type]++;

  ekin = TSharcAnalysis::GetReconstructedEnergy(position,det,denergy,penergy,ion[type]);
  exc = reaction[type]->GetExcEnergy(ekin*1e-3,thetalab*d2r,2);
  reaction[type]->SetExc(exc); // set CM excitation energy
  thetacm = reaction[type]->ConvertThetaLabToCm(thetalab*d2r,2);
  
  ekinthry = reaction[type]->GetTLab(thetalab*d2r,2)*1e3;
  dengthry = TSharcAnalysis::GetMeasuredEnergy(position,det,ekinthry,ion[type]).at(0); // use target thickness & default deadlayers to get expected denergy
  pengthry = TSharcAnalysis::GetMeasuredEnergy(position,det,ekinthry,ion[type],"",denergy).at(1); // use target, deadlayers & MEASURED denergy to get expected penergy
  if(pengthry<10.0)
    pengthry = 0.0;
  
 // if(type==4) printf("\n%6i.  type=%i   ekinthry=%5.0f   dengthry=%5.0f   pengthry=%5.0f   ekin=%5.0f    thetalab=%3.0f   exc=%.1f",entry,type,ekinthry,dengthry,pengthry,ekin,thetalab,exc*1e3);
  rthetalab = reaction[type]->ConvertThetaCmToLab(thetacm,3);
  
  //rvec = new TVector3();
  //rvec->SetMagThetaPhi(1.0,rthetalab,(180.0+phi)*d2r);
  rvec.Clear();
  rvec.SetMagThetaPhi(1.0,rthetalab,(180.0+phi)*d2r);

  rekin = reaction[type]->GetTLabFromThetaCm(pi-thetacm,3)*1e3;
  rekin = srim[type]->GetEnergy(rekin,TSharcAnalysis::GetTargetThickness(rthetalab)); // get recoil energy loss as it leaves the target
  rbeta = reaction[type]->AnalysisBeta(rekin*1e-3,3);
  /*
  if(type==0 && TMath::Abs(exc)<0.1){
     printf("\n\n%i.\t srim[500 MeV], x=[0.1 um] E=%.2f, x=[1 um] E=%.2f, x=[3 um] E=%.2f, x=[10 um] E=%.2f",entry,srim[type]->GetEnergy(500000.0,0.1),srim[type]->GetEnergy(500e3,1.0),srim[type]->GetEnergy(500e3,3.0),srim[type]->GetEnergy(500e3,10.0));
     printf("\n%i.\tthetalab=%4.2f\tthetacm=%4.2f\texc=%4.1f\trekin1=%4.2f\trekin=%4.2f\trbeta=%.3f\trthetalab=%4.3f\tthickness=%.3f",entry,thetalab,thetacm*r2d,exc*1e3,reaction[type]->GetTLabFromThetaCm(pi-thetacm,3)*1e3,rekin,rbeta,rthetalab*r2d,TSharcAnalysis::GetTargetThickness(rthetalab));
}*/
  // fill gammas 
  for(int t=0; t<tigress->GetAddbackMultiplicity(); t++){ 
    if(t==runinfo->addmax)    // THIS HAS BEEN HARD CODED AS 10 IN BRANCH VARIABLE DECLARATIONS
      break;
    addsize++;
    gamctot[type]++;
    gamcrun[type]++;
    
    TH = tigress->GetAddbackHit(t);
    tigdet[t] = TH->GetDetector();
    tigcry[t] = TH->GetCrystal();
    tigseg[t] = TH->GetSegment();
    tigeng[t] = TH->GetEnergy();
    //eadd[t] = TH->GetDoppler(rbeta,rvec); // this is the improved doppler correction i can also add a tvector!!
    tadd[t] = TH->GetTimeToTrigger();
    //thetadd[t] = TH->GetPosition().Theta()*r2d;
    tpos = TTigress::GetPosition(tigdet[t],tigcry[t],tigseg[t]);
    thetadd[t] = tpos.Theta()*r2d;
    eadd[t] = tigeng[t]/sqrt(1-rbeta*rbeta)*(1-rbeta*cos(tpos.Angle(rvec)));
    //printf("\n%i. %i. tigdet=%i\ttigcry=%i\ttigseg=%i\tTH_theta = %.3f\tT_theta = %.3f",entry,t,tigdet[t],tigcry[t],tigseg[t],TH->GetPosition().Theta()*r2d,thetadd[t]);
  }
  for(int t=0; t<tigress->Size(); t++){ 
    if(t==10)    // THIS HAS BEEN HARD CODED AS 10 IN BRANCH VARIABLE DECLARATIONS
      break;
    tsize++;
    
    TH = tigress->GetTigressHit(t);
    tdet[t] = TH->GetDetector();
    tcry[t] = TH->GetCrystal();
    tseg[t] = TH->GetSegment();
    teng[t] = TH->GetEnergy();
   // tdop[t] = TH->GetDoppler(rbeta,rvec); // this is the improved doppler correction i can also add a tvector!!
    ttim[t] = TH->GetTimeToTrigger();
    tpos = TTigress::GetPosition(tdet[t],tcry[t],tseg[t]);
    theth[t] = tpos.Theta()*r2d;
    tdop[t] = teng[t]/sqrt(1-rbeta*rbeta)*(1-rbeta*cos(tpos.Angle(rvec)));
  }
  exc*=1e3;
  thetacm*=r2d;

  return true;
}; 



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TTree *SetupOutTree(InputOptions *runinfo, TTree *tree, const char *outfile) {
  std::string filename = tree->GetCurrentFile()->GetName();
  size_t x = filename.find_last_of("analysis");
  std::string runnumber = (filename.substr(x+1,9));
//    printf("filename is %s, x is %i, runnumber is %s\n",filename.c_str(),x,runnumber.c_str());
//  TFile *file = new TFile(Form("%s_%s_%s.root",runinfo->treename.c_str(),runinfo->name.c_str(),runnumber.c_str()),"recreate");
  TFile *file = new TFile(outfile,"RECREATE");
  printf("Creating output file:  %s \n",file->GetName());
  file->cd();
  
  TTree *outtree = new TTree("EasyTree","EasyTree");
  // SHARC
  outtree->Branch("type",&type,"type/s");
  outtree->Branch("det",&det,"det/s");
  outtree->Branch("front",&front,"front/s");
  outtree->Branch("back",&back,"back/s");      
  outtree->Branch("thetalab",&thetalab,"thetalab/D");
  outtree->Branch("phi",&phi,"phi/D");
  outtree->Branch("phicorr",&phicorr,"phicorr/D");
  outtree->Branch("delcorr",&delcorr,"delcorr/D");
  
  outtree->Branch("dcharge",&dcharge,"dcharge/D");
  outtree->Branch("dchargeb",&dchargeb,"dchargeb/D");
  outtree->Branch("pcharge",&pcharge,"pcharge/D");
  outtree->Branch("dengalpha",&dengalpha,"dengalpha/D");
  outtree->Branch("denergy",&denergy,"denergy/D");
  outtree->Branch("penergy",&penergy,"penergy/D");
  outtree->Branch("exc",&exc,"exc/D");
  outtree->Branch("ekin",&ekin,"ekin/D");
  outtree->Branch("thetacm",&thetacm,"thetacm/D");
  
  outtree->Branch("rekin",&rekin,"rekin/D");
  outtree->Branch("rbeta",&rbeta,"rbeta/D");
  outtree->Branch("rthetalab",&rthetalab,"rthetalab/D");

  //RF
  outtree->Branch("phase",&phase,"phase/D");
 
   //Theoryz
  outtree->Branch("ekinthry",&ekinthry,"ekinthry/D"); 
  outtree->Branch("dengthry",&dengthry,"dengthry/D");
  outtree->Branch("pengthry",&pengthry,"pengthry/D");
  // TIGRESS
  outtree->Branch("addsize",&addsize,"addsize/s");
  outtree->Branch("tigdet",&tigdet,"tigdet[addsize]/s");
  outtree->Branch("tigcry",&tigcry,"tigcry[addsize]/s");
  outtree->Branch("tigseg",&tigseg,"tigseg[addsize]/s");

  outtree->Branch("tigeng",&tigeng,"tigeng[addsize]/D");  // addback, non doppler
  outtree->Branch("eadd",&eadd,"eadd[addsize]/D");        // addback; dopplerified.
  outtree->Branch("tadd",&tadd,"tadd[addsize]/D");
  outtree->Branch("thetadd",&thetadd,"thetadd[addsize]/D"); 
  
  // non addback tigress.
  outtree->Branch("tsize",&tsize,"tsize/s");
  outtree->Branch("tdet",&tdet,"tdet[tsize]/s");
  outtree->Branch("tcry",&tcry,"tcry[tsize]/s");
  outtree->Branch("tseg",&tseg,"tseg[tsize]/s");

  outtree->Branch("teng",&teng,"teng[tsize]/D");  // non-addback, non doppler
  outtree->Branch("tdop",&tdop,"tdop[tsize]/D");  // non-addback; dopplerified
  outtree->Branch("ttim",&ttim,"ttim[tsize]/D");
  outtree->Branch("theth",&theth,"theth[tsize]/D");


//  outtree->Branch("thetacry",&thetacry,"thetacry[addsize]/D");
  // TRIFOIL
  outtree->Branch("tbeam",&tbeam,"tbeam/I");
  outtree->Branch("beam",&beam,"beam/B");

  
  return outtree;
}

void CloseOutTree(TTree *outtree){

  TFile *outfile = outtree->GetCurrentFile();
  printf("Saving and closing output file: %s\n",outfile->GetName());
  outfile->cd();
  outtree->Write();
  outfile->Close();

  return;
}

Int_t CheckInputTree(const char *infile){

  TChain tree("AnalysisTree");
  tree.Add(infile);
  
  if(!tree.FindBranch("TSharc"))
    return -1;
  if(!tree.FindBranch("TTigress"))
    return -1;
  //  if(!tree->FindBranch("TTriFoil"))
  //    return;
  if(!tree.FindBranch("TRF"))
    return -1;
    
  return 1;
}

Int_t ProcessData(InputOptions *runinfo, const char *infile, const char *outfile){

  TChain *chain = new TChain("AnalysisTree");
  chain->Add(infile);
  //chain->Add(runinfo->chaindata.c_str());

	if(runinfo->nentries==0 ||  runinfo->nentries > chain->GetEntries())
	  runinfo->nentries = chain->GetEntries();
  
  int ntrees = chain->GetNtrees();
  int nChainEntries = chain->GetEntries();
  int treeNumber, lastTreeNumber = -1;
  int numTreeEvents,treeentries;
  bool switched = false;

	printf( DGREEN "\n\tUsing %i entries [starting at %i]... \n\n" RESET_COLOR "",runinfo->nentries,0);
	TStopwatch w;

  for( int i=0 ; i< runinfo->nentries; i++){
    chain->LoadTree(i);
    treeNumber = chain->GetTreeNumber();

    if (treeNumber != lastTreeNumber) {
      printf( DCYAN "\nChanging to tree number %d from %d at chain entry number %d / %d." RESET_COLOR "\n",treeNumber, lastTreeNumber, i,runinfo->nentries);
      lastTreeNumber = treeNumber;
    } else {
      continue;
    }

    TTree *tree = chain->GetTree();
    treeentries = tree->GetEntries();

    TTree *outtree = SetupOutTree(runinfo,tree,outfile); // done
    TSharc *sharc = 0;	
    TTigress *tigress = 0;
    TTriFoil *trifoil = 0;
    TRF *rf = 0;
    
    tree->SetBranchAddress("TSharc", &sharc);
    tree->SetBranchAddress("TTigress", &tigress);
    tree->SetBranchAddress("TTriFoil", &trifoil);
    tree->SetBranchAddress("TRF", &rf);


    for(int i=0; i<5; i++){
       typecrun[i] = 0; // reset type counter each run
       gamcrun[i] = 0; // reset gamma type counter each run
    }

    for(int j=0;j<treeentries;j++) {
      if(i+j > runinfo->nentries)
        break;
      if(sharc)   sharc->Clear();
      if(tigress) tigress->Clear();
      if(trifoil) trifoil->Clear();
      if(rf)      rf->Clear();
      
      if(!switched && i+j>=runinfo->piecewise){
        GainSwitch();
        switched = true;
      }

      tree->GetEntry(j);  
      if(SetBranchVariables(i+j,tigress,sharc,trifoil,rf)) 
        outtree->Fill();

      if(j%50000==0){
        printf( DYELLOW "\t on entry %i / %i            \r" RESET_COLOR "",j,treeentries);
        fflush(stdout);
      }
    }
    printf( DGREEN "\t on entry %i / %i            \n" RESET_COLOR "",treeentries,treeentries);
    fflush(stdout);
    CloseOutTree(outtree);
    i += treeentries -5; // what the hell??
    
    // printf off some run stats and write to log file is requested
    printf( DCYAN "\n\t Particles [type]:\t%8i\t%8i\t%8i\t%8i",0,1,2,3);         
    printf( DYELLOW "\n\t\tRUN \t\t");         
    for(int i=0; i<4; i++)
        printf("%8i\t",typecrun[i]);
    printf( DGREEN "\n\t\tTOT \t\t");         
    for(int i=0; i<4; i++)
        printf("%8i\t",typectot[i]);
    printf( DCYAN "\n\t Gammas    [type]:\t%8i\t%8i\t%8i\t%8i",0,1,2,3);         
    printf( DYELLOW "\n\t\tRUN \t\t");         
    for(int i=0; i<4; i++)
        printf("%8i\t",gamcrun[i]);
    printf( DGREEN "\n\t\tTOT \t\t");         
    for(int i=0; i<4; i++)
        printf("%8i\t",gamctot[i]);
    printf( RESET_COLOR "\n\n");

    if(runinfo->logfile.length()){
      static std::ofstream logofile(runinfo->logfile.c_str());
      std::string strfile = chain->GetFile()->GetName();
      logofile << strfile.substr(strfile.find("analysis"),17).c_str() << "\t";
      logofile << typecrun[0] << "\t" << typecrun[1] << "\t" << typecrun[2] << "\t" << typecrun[3] << "\t";
      logofile << gamcrun[0] << "\t" << gamcrun[1] << "\t" << gamcrun[2] << "\t" << gamcrun[3] << "\n";
      logofile << std::flush;
      if(treeNumber==chain->GetNtrees()-1)
        logofile.close();
    }
  }

  printf( DGREEN "\n\n...Aaaaand we're done!" RESET_COLOR "\n\n");
  return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t InitVariables(InputOptions *runinfo){
 
  if(!runinfo->A || !runinfo->TargetThick || !runinfo->EperU || runinfo->cutfile.length()==0)
     return false;

  const char *calfilename = runinfo->calfile.c_str();
  printf( DCYAN "Cal file loaded: \n\t" RESET_COLOR "'%s' \n",calfilename);
  TChannel::ReadCalFile(calfilename);
  if(TChannel::GetNumberOfChannels()==0){
    printf("Couldn't find any channels... exiting.\n");
    return false;
  }
	printf("\t%i channels available\n",TChannel::GetNumberOfChannels());
	
// reaction stuff ///////////////////////////////////////////////////////////////////////////////////////////////
 
  TSharcAnalysis::SetTarget(runinfo->offset.X(),runinfo->offset.Y(),runinfo->offset.Z(),runinfo->TargetThick,runinfo->TargetMat.c_str()); // pos(mm), thickness(um)
  char beamnuc[6], reconuc[6], lightnuc[6];
  sprintf(beamnuc,"sr%i",runinfo->A);
  sprintf(reconuc,"sr%i",runinfo->A+1);
  sprintf(lightnuc,"sr%i",runinfo->A-1);
	
  TSharcAnalysis::SetBadStrips(runinfo->badstripsfile.c_str());
  double beame = TSharcAnalysis::GetBeamEnergyInTarget(beamnuc,runinfo->A*runinfo->EperU*1e3)*1e-3;
	reaction[0] = new TReaction(beamnuc,"d","p",reconuc,beame,0.0,true); // (d,p)
	reaction[0]->Print();
  reaction[1] = new TReaction(beamnuc,"p","p",beamnuc,beame,0.0,true); // (p,p)
	reaction[1]->Print();
	reaction[2] = new TReaction(beamnuc,"d","d",beamnuc,beame,0.0,true); // (d,d)
	reaction[2]->Print();
	reaction[3] = new TReaction(beamnuc,"c12","c12",beamnuc,beame,0.0,true); // (c,c)
	reaction[3]->Print();
	reaction[4] = new TReaction(beamnuc,"d","t",lightnuc,beame,0.0,true); // (d,t)
	reaction[4]->Print();
  
  if(runinfo->A==96){
    srim[0] = new TSRIM(Form("%s_in_%s.txt",beamnuc,runinfo->TargetMat.c_str())); // (d,p)
    srim[1] = new TSRIM(Form("%s_in_%s.txt",beamnuc,runinfo->TargetMat.c_str())); // (p,p)
    srim[4] = new TSRIM(Form("%s_in_%s.txt",lightnuc,runinfo->TargetMat.c_str()));// (d,t)
  }else if(runinfo->A==95){
    srim[0] = new TSRIM(Form("%s_in_%s.txt",reconuc,runinfo->TargetMat.c_str())); // (d,p)
    srim[1] = new TSRIM(Form("%s_in_%s.txt",beamnuc,runinfo->TargetMat.c_str())); // (p,p)
    srim[4] = new TSRIM(Form("%s_in_%s.txt",lightnuc,runinfo->TargetMat.c_str()));// (d,t)
  }else if(runinfo->A==94){
    srim[0] = new TSRIM(Form("%s_in_%s.txt",reconuc,runinfo->TargetMat.c_str())); // (d,p)
    srim[1] = new TSRIM(Form("%s_in_%s.txt",beamnuc,runinfo->TargetMat.c_str())); // (p,p)
    srim[4] = new TSRIM(Form("%s_in_%s.txt",beamnuc,runinfo->TargetMat.c_str())); // (d,t)
  }
  srim[2] = srim[1]; // (d,d)
  srim[3] = srim[1]; // (c,c)

// cut stuff ////////////////////////////////////////////////////////////////////////////////////////////////////

 	TFile *cutfile = new TFile(runinfo->cutfile.c_str(),"READ");
  if(!cutfile->IsOpen())
    return false;
  
  for(int d=5; d<=8; d++){
     if(d==5){
        if(runinfo->A>94){ // takes care of piecewise calibration
           IsProtPid5_1 = (TCutG*)cutfile->Get("IsProtPid5_1");  
           IsProtPid5_2 = (TCutG*)cutfile->Get("IsProtPid5_2");  
           IsDeutPid5_1 = (TCutG*)cutfile->Get("IsDeutPid5_1"); 
           IsDeutPid5_2 = (TCutG*)cutfile->Get("IsDeutPid5_2"); 
           IsTritPid5_1 = (TCutG*)cutfile->Get("IsTritPid5_1"); 
           IsTritPid5_2 = (TCutG*)cutfile->Get("IsTritPid5_2"); 
           IsCarbPid5_1 = (TCutG*)cutfile->Get("IsCarbPid5_1"); 
           IsCarbPid5_2 = (TCutG*)cutfile->Get("IsCarbPid5_2"); 
           printf( DBLUE "Piecewise Cuts loaded:" RESET_COLOR "\n\tIsProtPid5_1\t= 0x%08x\tIsDeutPid5_1\t= 0x%08x\tIsTritPid5_1\t= 0x%08x\tIsCarbPid5_1\t= 0x%08x\n",
                IsProtPid5_1,IsDeutPid5_1,IsTritPid5_1,IsCarbPid5_1);
           printf( DBLUE "Piecewise Cuts loaded:" RESET_COLOR "\n\tIsProtPid5_2\t= 0x%08x\tIsDeutPid5_2\t= 0x%08x\tIsTritPid5_2\t= 0x%08x\tIsCarbPid5_2\t= 0x%08x\n",
                IsProtPid5_2,IsDeutPid5_2,IsTritPid5_2,IsCarbPid5_2);
        } else{
           IsProtPid5_1 = (TCutG*)cutfile->Get("IsProtPid5");  
           IsProtPid5_2 = (TCutG*)cutfile->Get("IsProtPid5");  
           IsDeutPid5_1 = (TCutG*)cutfile->Get("IsDeutPid5"); 
           IsDeutPid5_2 = (TCutG*)cutfile->Get("IsDeutPid5"); 
           IsTritPid5_1 = (TCutG*)cutfile->Get("IsTritPid5"); 
           IsTritPid5_2 = (TCutG*)cutfile->Get("IsTritPid5"); 
           IsCarbPid5_1 = (TCutG*)cutfile->Get("IsCarbPid5"); 
           IsCarbPid5_2 = (TCutG*)cutfile->Get("IsCarbPid5"); 
        }
        // cut is pointing to first set of cuts in piecewise calibration
        IsProtPid[d-5] = IsProtPid5_1;
        IsDeutPid[d-5] = IsDeutPid5_1;
        IsTritPid[d-5] = IsTritPid5_1;
        IsCarbPid[d-5] = IsCarbPid5_1;
    } else {
        IsProtPid[d-5] = (TCutG*)cutfile->Get(Form("IsProtPid%i",d));  
        IsDeutPid[d-5] = (TCutG*)cutfile->Get(Form("IsDeutPid%i",d)); 
        IsTritPid[d-5] = (TCutG*)cutfile->Get(Form("IsTritPid%i",d)); 
        IsCarbPid[d-5] = (TCutG*)cutfile->Get(Form("IsCarbPid%i",d)); 
    }
    printf( DCYAN "Cuts loaded:" RESET_COLOR "\n\tIsProtPid[%i]\t= 0x%08x\tIsDeutPid[%i]\t= 0x%08x\tIsTritPid[%i]\t= 0x%08x\tIsCarbPid[%i]\t= 0x%08x\n",
          d-5,IsProtPid[d-5],d-5,IsDeutPid[d-5],d-5,IsTritPid[d-5],d-5,IsCarbPid[d-5]);
  }
  // make all these cuts overlap so that they can be combined in a fit and draw nicely      
  IsProtElastic = (TCutG*)cutfile->Get("IsProtKin");  
  IsDeutElastic = (TCutG*)cutfile->Get("IsDeutKin"); // this cut extends all the way to the proton cut so that there is no gap between them (double peak fitting works better)  
  IsCarbElastic = (TCutG*)cutfile->Get("IsCarbKin"); // this cut is for carboon  
  IsTritElastic = (TCutG*)cutfile->Get("IsTritKin"); // this cut is for tritoon  
  printf( DCYAN "Cuts loaded:" RESET_COLOR "\n\tIsProtElastic\t= 0x%08x\n\tIsDeutElastic\t= 0x%08x\n\tIsCarbElastic\t= 0x%08x\n\tIsTritElastic\t= 0x%08x\n",IsProtElastic,IsDeutElastic,IsCarbElastic,IsTritElastic);
 
  //IsGoodCharge = (TCutG*)cutfile->Get("FrontBackCharge");
  IsGoodCharge = (TCutG*)cutfile->Get("GoodFrontBack"); // slightly larger. More background, more data.
  printf( DCYAN "Cut loaded:" RESET_COLOR "\n\tIsGoodCharge\t= 0x%08x\n\t\n\n",IsGoodCharge);

  cutfile->Close();   

  if(!IsProtPid[0] || !IsDeutPid[0]){
    printf("Couldn't find enough cuts... exiting.\n");
    return false;
  }
  
  return true;
}


void trim(std::string * line, const std::string trimChars = " \f\n\r\t\v") {
//Removes the the string "trimCars" from  the string 'line'
   if (line->length() == 0)
      return;
   std::size_t found = line->find_first_not_of(trimChars);
   if (found != std::string::npos)
      *line = line->substr(found, line->length());
   found = line->find_last_not_of(trimChars);
   if (found != std::string::npos)
      *line = line->substr(0, found + 1);
   return;
}



InputOptions *SetOptions(const char *filename) {
   //Makes TChannels from a cal file.
   std::string infilename;
   infilename.append(filename);

   InputOptions *runinfo = new InputOptions;

   runinfo->treename = "";
   runinfo->A = 0;
   runinfo->nentries = -1;
   runinfo->ecut = 0;
   runinfo->addmax = 0;
   runinfo->piecewise = 0;
   runinfo->EperU = 0.00;
   runinfo->dbox_coeffs[0] = 0.0;
   runinfo->dbox_coeffs[1] = 1.0;
   runinfo->ubox_coeffs[0] = 0.0;
   runinfo->ubox_coeffs[1] = 1.0;
   runinfo->uqqq_coeffs[0] = 0.0;
   runinfo->uqqq_coeffs[1] = 1.0;
   runinfo->pad5_coeffs1[0] = 0.0;
   runinfo->pad5_coeffs1[1] = 1.0;
   runinfo->pad5_coeffs2[0] = 0.0;
   runinfo->pad5_coeffs2[1] = 1.0;
   runinfo->pad_coeffs[0][0] = 0.0;
   runinfo->pad_coeffs[0][1] = 1.0;
   runinfo->pad_coeffs[1][0] = 0.0;
   runinfo->pad_coeffs[1][1] = 1.0;
   runinfo->pad_coeffs[2][0] = 0.0;
   runinfo->pad_coeffs[2][1] = 1.0;
   runinfo->pad_coeffs[3][0] = 0.0;
   runinfo->pad_coeffs[3][1] = 1.0;
   runinfo->TargetThick= 0.00;
   runinfo->TargetMat= "";
   runinfo->cutfile = "";
   runinfo->calfile = "";
   runinfo->logfile = "";
   runinfo->offset.SetXYZ(0.00,0.00,0.00);
   runinfo->badstripsfile = "";
//   runinfo->typespec = "";
   runinfo->chaindata = "";
   runinfo->name = "";

   ifstream infile;
   infile.open(infilename.c_str());
   if(!infile || infilename.length()==0) {
      printf( DRED "\n\n\t Error :  Could not open file." RESET_COLOR "\n\n");
      return runinfo;
   }

   std::string line;
   int linenumber = 0;

   bool brace_open = false;
   int detector = 0;
   std::string name;

   //std::pair < int, int >pixel = std::make_pair(0, 0);

   while (std::getline(infile, line)) {
      linenumber++;
      trim(&line);
      int comment = line.find("//");
      if (comment != std::string::npos) {
         line = line.substr(0, comment);
      }
      if (!line.length())
         continue;
      int openbrace = line.find("{");
      int closebrace = line.find("}");
		int colon = line.find(":");

		if(openbrace  == std::string::npos &&
			closebrace == std::string::npos &&
			colon  == std::string::npos)
  		continue;
		//printf("line : %s\n",line.c_str());

      //*************************************//
      if (closebrace != std::string::npos) {

         brace_open = false;
         return runinfo; 
      }
      //*************************************//
      if (openbrace != std::string::npos) {
         brace_open = true;
         name = line.substr(0, openbrace).c_str();
         runinfo->name = name;
      }
      //*************************************//
      if (brace_open) {
         int ntype = line.find(":");
         if (ntype != std::string::npos) {
            std::string type = line.substr(0, ntype);
            line = line.substr(ntype + 1, line.length());
            trim(&line);
            std::istringstream ss(line);
            int j = 0;
            while (type[j]) {
               char c = *(type.c_str() + j);
               c = toupper(c);
               type[j++] = c;
            }
            //printf("type = %s\n",type.c_str());
            if(type.compare("TREENAME")==0) {
               runinfo->treename = line;
            } else if(type.compare("A")==0) {
               int tempint; ss>>tempint;
               runinfo->A = tempint;
            } else if(type.compare("NENTRIES")==0) {
               int tempint; ss>>tempint;
               runinfo->nentries = tempint;
            } else if(type.compare("ECUT")==0) {
              double tempdub; ss>>tempdub;
               runinfo->ecut = tempdub;     
            } else if(type.compare("EPERU")==0) {
              double tempdub; ss>>tempdub;
               runinfo->EperU = tempdub;                              
            } else if(type.compare("TARGETTHICK")==0) {
               double tempdub; ss>>tempdub;
               runinfo->TargetThick = tempdub;                              
            } else if(type.compare("TARGETMAT")==0) {
               runinfo->TargetMat = line;                              
            } else if(type.compare("OFFSET")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==3) {
                  runinfo->offset.SetXYZ(tempvec.at(0),tempvec.at(1),tempvec.at(2));
               }
            } else if(type.compare("DBOX_COEFFS")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==2) {
                  runinfo->dbox_coeffs[0] = tempvec.at(0);
                  runinfo->dbox_coeffs[1] = tempvec.at(1);
               }
            } else if(type.compare("UBOX_COEFFS")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==2) {
                  runinfo->ubox_coeffs[0] = tempvec.at(0);
                  runinfo->ubox_coeffs[1] = tempvec.at(1);
               }
            } else if(type.compare("UQQQ_COEFFS")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==2) {
                  runinfo->uqqq_coeffs[0] = tempvec.at(0);
                  runinfo->uqqq_coeffs[1] = tempvec.at(1);
               }
            } else if(type.compare("PAD5_COEFFS1")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==2) {
                  runinfo->pad5_coeffs1[0] = tempvec.at(0);
                  runinfo->pad5_coeffs1[1] = tempvec.at(1);
                  // set coefficients to begin with this value
                  runinfo->pad_coeffs[0][0] = tempvec.at(0);
                  runinfo->pad_coeffs[0][1] = tempvec.at(1);
               }
            } else if(type.compare("PAD5_COEFFS2")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==2) {
                  runinfo->pad5_coeffs2[0] = tempvec.at(0);
                  runinfo->pad5_coeffs2[1] = tempvec.at(1);
               }
            } else if(type.compare("PAD6_COEFFS")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==2) {
                  runinfo->pad_coeffs[1][0] = tempvec.at(0);
                  runinfo->pad_coeffs[1][1] = tempvec.at(1);
               }
            } else if(type.compare("PAD7_COEFFS")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==2) {
                  runinfo->pad_coeffs[2][0] = tempvec.at(0);
                  runinfo->pad_coeffs[2][1] = tempvec.at(1);
               }
            } else if(type.compare("PAD8_COEFFS")==0) {
               double value;
               std::vector<double> tempvec; 
               while (ss >> value) {   tempvec.push_back(value); }
               if(tempvec.size()==2) {
                  runinfo->pad_coeffs[3][0] = tempvec.at(0);
                  runinfo->pad_coeffs[3][1] = tempvec.at(1);
               }
            } else if(type.compare("NAME")==0) {
               runinfo->name = line;
            } else if(type.compare("CHAINDATA")==0) {
               runinfo->chaindata = line;
            } else if(type.compare("CUTFILE")==0) {
               runinfo->cutfile = line;
            } else if(type.compare("CALFILE")==0) {
               runinfo->calfile = line;
            } else if(type.compare("LOGFILE")==0) {
               runinfo->logfile = line;
            } else if(type.compare("BADSTRIPSFILE")==0) {
               runinfo->badstripsfile = line;
            } else if(type.compare("ADDMULTIMAX")==0) {
               double val; ss>> val;
               runinfo->addmax = val;
            } else if(type.compare("PIECEWISECAL")==0) {
               double val; ss>> val;
               runinfo->piecewise = val;
            }
         }
      }

   }

   return runinfo;
}

