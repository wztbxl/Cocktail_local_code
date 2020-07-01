#include "/star/u/syang/Macro/headers.h"
#include "/star/u/syang/Macro/function.C"

void plotCocktail(TString cenName = "080", Bool_t withRho = kFALSE){

	Int_t mCenIdx;
	if (cenName.EqualTo("080")) {
		mCenIdx = 0;
	}
	else if (cenName.EqualTo("010")) {
		mCenIdx = 1;
	}
	else if (cenName.EqualTo("1040")) {
		mCenIdx = 2;
	}
	else if (cenName.EqualTo("4080")) {
		mCenIdx = 3;
	}
	else if (cenName.EqualTo("4060")) {
		mCenIdx = 4;
	}
	else if (cenName.EqualTo("6080")) {
		mCenIdx = 5;
	}
	else if (cenName.EqualTo("6070")) {
		mCenIdx = 6;
	}
	else if (cenName.EqualTo("7080")) {
		mCenIdx = 7;
	}
	else{
		cout<<"Centrality name is wrong !"<<endl;
		return;
	}

	gStyle->SetOptFit(1111);

	const Double_t dY = 2.4;  // input parent particle rapidity range [-1.2, 1.2]
	const Double_t dM = 0.005;

	/*
	 * Centrality Sequence:
	 * 0 - 0-80%
	 * 1 - 0-10%
	 * 2 - 10-40%
	 * 3 - 40-80%
	 * 4 - 40-60%
	 * 5 - 60-80%
	 * 6 - 60-70%
	 * 7 - 70-80%
	 */
	const Int_t    nCenBins = 8;
	const Int_t    CentralityLow[nCenBins] = {0,  0,  10, 40, 40, 60, 60, 70};  // %
	const Int_t    CentralityHi[nCenBins]  = {80, 10, 40, 80, 60, 80, 70, 80};

	//const Double_t Nbinary[nCenBins] = {350.08, 1146.12, 467.44, 63.05, 103.61, 22.48, 31.71, 13.26}; // UU@193 GeV collisions
	const Double_t Nbinary[nCenBins] = {256.52, 1146.12, 467.44, 63.05, 103.61, 22.48, 31.71, 13.26}; // Au+Au 54 GeV collisions
	//const Double_t Nbinary[nCenBins] = {252.90, 1146.12, 467.44, 63.05, 103.61, 22.48, 31.71, 13.26}; // Au+Au 62 GeV collisions
	const Double_t dNdy_pi0_UU[nCenBins] = {115.756, 358.484, 165.895, 33.3489, 51.6463, 15.7225, 21.053, 10.4558};
	const Double_t dNdy_pi0_54[nCenBins] = {67.6539};//Qian's slides mean of pi+/pi-
	const Double_t dNdy_pi0_62[nCenBins] = {77.2};//62.4 

	//const Double_t dNdy_eta_AuAu_scale[nCenBins] = {1.34818, 1.17589, 1.29988, 1.59568, 1.60484, 1.73864}; // scale to PHENIX eta/pi0 ratio
	const Double_t dNdy_eta_AuAu_scale[nCenBins] = {1.34818, 1.17589, 1.29988, 1.59568, 1.60484, 1.73864, 1.73864, 1.73864}; // scale to PHENIX eta/pi0 ratio, use 60-80% factor for 60-70%, 70-80%
	//const Double_t dNdy_pi0_AuAu[nCenBins] = {98.5, 288, 135, 31.4, 47.5, 16}; // from Chi's direct photon analysis
    //const Double_t dNdy_eta_AuAu_raw[nCenBins] = {7.86, 22.9815, 10.7726, 2.50562, 3.79036, 1.27675}; // from Chi's direct photon analysis

	//const Double_t dNtotcc = 1.90e-2*Nbinary[mCenIdx]; // 0.8mb/42mb
	//const Double_t dNtotcc = 3.58e-3*252.90;//62GeV
  	const Double_t dNtotcc = 2.86e-3*Nbinary[mCenIdx];//54GeV how to get the total cross section 0.1mb/35mb
	const Double_t dNtotbb = 8.8e-5*Nbinary[mCenIdx]; // 3.7ub/42mb
	const Double_t dNtotdy = 1.e-6*Nbinary[mCenIdx]; // 42nb/42mb
	const Double_t BR_dy = 0.033575; //this brach ratio already included in the cross-section 

	Double_t dNdy_pi0,dNdy_eta,dNdy_omega,dNdy_phi,dNdy_etaprim,dNdy_jpsi,dNdy_rho,dNdy_psi;
	ifstream indata("./dndy.dat"); //this dndy is AuAu@200GeV minibias(0-80%) results
	if(indata.is_open()){
		indata>>dNdy_pi0>>dNdy_eta>>dNdy_rho>>dNdy_omega>>dNdy_phi>>dNdy_etaprim>>dNdy_jpsi>>dNdy_psi;
		cout<<"AuAu@200 GeV 0-80% particle yields:"<<endl;
		cout<<"dNdy_pi0: "<<dNdy_pi0<<endl;
		cout<<"dNdy_eta: "<<dNdy_eta<<endl;
		cout<<"dNdy_rho: "<<dNdy_rho<<endl;
		cout<<"dNdy_omega: "<<dNdy_omega<<endl;
		cout<<"dNdy_etaprim: "<<dNdy_etaprim<<endl;
		cout<<"dNdy_phi: "<<dNdy_phi<<endl;
		cout<<"dNdy_jpsi: "<<dNdy_jpsi<<endl;
		cout<<"dNdy_psi: "<<dNdy_psi<<endl;
	}else{
		cout<<"Failed to load dNdy parameters !"<<endl;
		return;
	}
	indata.close();

	Double_t dNdy_pi0Err = 0.08;
	Double_t dNdy_etaErr = (0.543-0.48)/0.48; // old numeber: 0.3
	Double_t dNdy_omegaErr = 0.33;
	Double_t dNdy_phiErr = 0.1;
	Double_t dNdy_etaprimErr = 1.;
	Double_t dNdy_jpsiErr = 0.15;
	Double_t dNdy_rhoErr = 0.42;
	Double_t dNdy_psiErr = 0.27;
	Double_t dNdy_ccErr = 0.15;
	Double_t dNdy_bbErr = 0.3;
	Double_t dNdy_dyErr = 0.3;

	//Double_t scaleFactor = dNdy_pi0_UU[mCenIdx]/dNdy_pi0;
	//dNdy_pi0     *= scaleFactor;
	//dNdy_eta     *= dNdy_eta_AuAu_scale[mCenIdx]*scaleFactor;
	//dNdy_rho     *= scaleFactor; // *** special *** 
	//dNdy_omega   *= scaleFactor;
	//dNdy_phi     *= scaleFactor;
	//dNdy_etaprim *= scaleFactor;
	//dNdy_jpsi    *= Nbinary[mCenIdx]/291.90; // 291.90 is AuAu minibias Nbinary
	//dNdy_psi     *= Nbinary[mCenIdx]/291.90; // 291.90 is AuAu minibias Nbinary

	//cout<<"****************"<<endl;
	//cout<<"dNdy_eta@Chi: "<<dNdy_eta_AuAu_scale[mCenIdx]*dNdy_eta_AuAu_raw[mCenIdx]*dNdy_pi0_UU[mCenIdx]/dNdy_pi0_AuAu[mCenIdx]<<"      dNdy_eta@Shuai: "<<dNdy_eta<<endl;
	//cout<<"****************"<<endl;
	//return;

	cout<<endl;
	//cout<<"AuAu->UU Scale Factor: "<<scaleFactor<<endl;
	cout<<"dNdy_pi0: "<<dNdy_pi0<<endl;
	cout<<"dNdy_eta: "<<dNdy_eta<<endl;
	cout<<"dNdy_rho: "<<dNdy_rho<<endl;
	cout<<"dNdy_omega: "<<dNdy_omega<<endl;
	cout<<"dNdy_etaprim: "<<dNdy_etaprim<<endl;
	cout<<"dNdy_phi: "<<dNdy_phi<<endl;
	cout<<"dNdy_jpsi: "<<dNdy_jpsi<<endl;
	cout<<"dNdy_psi: "<<dNdy_psi<<endl;

	const Double_t BR_pi0dal=1.174e-2;
	const Double_t BR_eta2ee = 2.7e-5;
	//const Double_t BR_etadal=7.e-3;
	const Double_t BR_etadal=6.9e-3;
	const Double_t BR_rho2ee=4.72e-5;
	const Double_t BR_omegadal=7.7e-4;
	const Double_t BR_omega2ee=7.28e-5;
	//const Double_t BR_etaprimdal=9.e-4;
	const Double_t BR_etaprimdal=4.7e-4; // calculate by Landsberg's formula
	const Double_t BR_phidal=1.15e-4;
	const Double_t BR_phi2ee=2.954e-4;
	const Double_t BR_jpsi2ee=5.94e-2;
	const Double_t BR_psi2ee=7.72e-3;

	const Double_t PtMin = 0.;
	const Double_t PtMax = 10.;

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
	TPDF *ps; 
	TString dir;
	system(Form("mkdir -p output/Cen%d%d", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	dir = Form("output/Cen%d%d", CentralityLow[mCenIdx], CentralityHi[mCenIdx]);
	if(withRho){
		ps = new TPDF(Form("%s/cen%d%d_cocktail_withRho.pdf", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]), 111);
	}else{
		ps = new TPDF(Form("%s/cen%d%d_cocktail_withoutRho.pdf", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]), 111);
	}
	ps->Off();

	//dalitz decay
	TFile *fpi0dal = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/pi0dalitzAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fpi0dal->IsOpen())){
		cout<<"Failed to open the pi0 dalitz decay file !"<<endl;
		return;
	}
	TFile *fetadal = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/etadalitzAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fetadal->IsOpen())){
		cout<<"Failed to open the eta dalitz decay file !"<<endl;
		return;
	}
	TFile *fomegadal = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/omegadalitzAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fomegadal->IsOpen())){
		cout<<"Failed to open the omega dalitz decay file !"<<endl;
		return;
	}
	TFile *fetaprimdal = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/etaprimdalitzAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fetaprimdal->IsOpen())){
		cout<<"Failed to open the etaprim dalitz decay file !"<<endl;
		return;
	}
	TFile *fphidal = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/phidalitzAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fphidal->IsOpen())){
		cout<<"Failed to open the phi dalitz decay file !"<<endl;
		return;
	}

	//twobody decay
	/*TFile *frho2ee = new TFile("/star/u/chiyang/bnl/dilepton/auau200/cocktail/run11_tpc20eta36phi/rho2ee.root");
	if(!(frho2ee->IsOpen())){
		cout<<"Failed to open the rho two body decay file !"<<endl;
		return;
	}*/
	TFile *fomega2ee = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/omega2eeAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fomega2ee->IsOpen())){
		cout<<"Failed to open the omega two body decay file !"<<endl;
		return;
	}
	TFile *fphi2ee = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/phi2eeAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fphi2ee->IsOpen())){
		cout<<"Failed to open the phi two body decay file !"<<endl;
		return;
	}
	TFile *fjpsi2ee = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/jpsi2eeAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fjpsi2ee->IsOpen())){
		cout<<"Failed to open the jpsi two body decay file !"<<endl;
		return;
	}
	/*TFile *fpsi2ee = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/Cen%d%d/psi/psi2eeAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fpsi2ee->IsOpen())){
		cout<<"Failed to open the psi two body decay file !"<<endl;
		return;
	}*/

	//heavy flavor decay cc,bb,drell-yan
	TFile *fcc = new TFile(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/input/54Result/Cen%d%d/cc2eeAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fcc->IsOpen())){
		cout<<"Failed to open the ccbar decay file !"<<endl;
		return;
	}
	/*TFile *fbb = new TFile(Form("/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/heavyflavor/bb2ee/output/Cen%d%d/bb2eeAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fbb->IsOpen())){
		cout<<"Failed to open the bbbar decay file !"<<endl;
		return;
	}
	TFile *fdy = new TFile(Form("/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/heavyflavor/dy2ee/output/Cen%d%d/dy2eeAll.root", CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	if(!(fdy->IsOpen())){
		cout<<"Failed to open the drell-yan process file !"<<endl;
		return;
	}*/

	//MC: -1<= Y_{mother particle (except heavy flavors which cover full rapidity range)} <=1;
	//MCAcc0: -1<= Y_{ee} <=1;
	//MCAcc1: -1<= Y_{ee} <=1 && -1<= eta_{e} <=1 && pT_{e}>=0.2
	//RCAcc1: -1<= Y_{ee} <=1 && -1<= eta_{e} <=1 && pT_{e}>=0.2 && weight the single track efficiency

	//dalitz decay
	TH2D *hMCMvsPtpi0dalitz = (TH2D*)fpi0dal->Get("hMCMvsPtpi0dalitz");
	hMCMvsPtpi0dalitz->Sumw2();
	TH2D *hMCAcc0MvsPtpi0dalitz = (TH2D*)fpi0dal->Get("hMCAcc0MvsPtpi0dalitz");
	hMCAcc0MvsPtpi0dalitz->Sumw2();
	TH2D *hMCAcc1MvsPtpi0dalitz = (TH2D*)fpi0dal->Get("hMCAcc1MvsPtpi0dalitz");
	hMCAcc1MvsPtpi0dalitz->Sumw2();
	TH2D *hRCAcc1MvsPt3Dpi0dalitz = (TH2D*)fpi0dal->Get("hRCAcc1MvsPt3Dpi0dalitz");
	hRCAcc1MvsPt3Dpi0dalitz->Sumw2();
	TH2D *hUpMCAcc1MvsPtpi0dalitz = (TH2D *)hMCAcc1MvsPtpi0dalitz->Clone("hUpMCAcc1MvsPtpi0dalitz");
	TH2D *hLowMCAcc1MvsPtpi0dalitz = (TH2D *)hMCAcc1MvsPtpi0dalitz->Clone("hLowMCAcc1MvsPtpi0dalitz");

	TH2D *hMCMvsPtetadalitz = (TH2D*)fetadal->Get("hMCMvsPtetadalitz");
	hMCMvsPtetadalitz->Sumw2();
	TH2D *hMCAcc0MvsPtetadalitz = (TH2D*)fetadal->Get("hMCAcc0MvsPtetadalitz");
	hMCAcc0MvsPtetadalitz->Sumw2();
	TH2D *hMCAcc1MvsPtetadalitz = (TH2D*)fetadal->Get("hMCAcc1MvsPtetadalitz");
	hMCAcc1MvsPtetadalitz->Sumw2();
	TH2D *hRCAcc1MvsPt3Detadalitz = (TH2D*)fetadal->Get("hRCAcc1MvsPt3Detadalitz");
	hRCAcc1MvsPt3Detadalitz->Sumw2();
	TH2D *hUpMCAcc1MvsPtetadalitz = (TH2D *)hMCAcc1MvsPtetadalitz->Clone("hUpMCAcc1MvsPtetadalitz");
	TH2D *hLowMCAcc1MvsPtetadalitz = (TH2D *)hMCAcc1MvsPtetadalitz->Clone("hLowMCAcc1MvsPtetadalitz");

	TH2D *hMCMvsPtomegadalitz = (TH2D*)fomegadal->Get("hMCMvsPtomegadalitz");
	hMCMvsPtomegadalitz->Sumw2();
	TH2D *hMCAcc0MvsPtomegadalitz = (TH2D*)fomegadal->Get("hMCAcc0MvsPtomegadalitz");
	hMCAcc0MvsPtomegadalitz->Sumw2();
	TH2D *hMCAcc1MvsPtomegadalitz = (TH2D*)fomegadal->Get("hMCAcc1MvsPtomegadalitz");
	hMCAcc1MvsPtomegadalitz->Sumw2();
	TH2D *hRCAcc1MvsPt3Domegadalitz = (TH2D*)fomegadal->Get("hRCAcc1MvsPt3Domegadalitz");
	hRCAcc1MvsPt3Domegadalitz->Sumw2();
	TH2D *hUpMCAcc1MvsPtomegadalitz = (TH2D *)hMCAcc1MvsPtomegadalitz->Clone("hUpMCAcc1MvsPtomegadalitz");
	TH2D *hLowMCAcc1MvsPtomegadalitz = (TH2D *)hMCAcc1MvsPtomegadalitz->Clone("hLowMCAcc1MvsPtomegadalitz");

	TH2D *hMCMvsPtetaprimdalitz = (TH2D*)fetaprimdal->Get("hMCMvsPtetaprimdalitz");
	hMCMvsPtetaprimdalitz->Sumw2();
	TH2D *hMCAcc0MvsPtetaprimdalitz = (TH2D*)fetaprimdal->Get("hMCAcc0MvsPtetaprimdalitz");
	hMCAcc0MvsPtetaprimdalitz->Sumw2();
	TH2D *hMCAcc1MvsPtetaprimdalitz = (TH2D*)fetaprimdal->Get("hMCAcc1MvsPtetaprimdalitz");
	hMCAcc1MvsPtetaprimdalitz->Sumw2();
	TH2D *hRCAcc1MvsPt3Detaprimdalitz = (TH2D*)fetaprimdal->Get("hRCAcc1MvsPt3Detaprimdalitz");
	hRCAcc1MvsPt3Detaprimdalitz->Sumw2();
	TH2D *hUpMCAcc1MvsPtetaprimdalitz = (TH2D *)hMCAcc1MvsPtetaprimdalitz->Clone("hUpMCAcc1MvsPtetaprimdalitz");
	TH2D *hLowMCAcc1MvsPtetaprimdalitz = (TH2D *)hMCAcc1MvsPtetaprimdalitz->Clone("hLowMCAcc1MvsPtetaprimdalitz");

	TH2D *hMCMvsPtphidalitz = (TH2D*)fphidal->Get("hMCMvsPtphidalitz");
	hMCMvsPtphidalitz->Sumw2();
	TH2D *hMCAcc0MvsPtphidalitz = (TH2D*)fphidal->Get("hMCAcc0MvsPtphidalitz");
	hMCAcc0MvsPtphidalitz->Sumw2();
  	cout<<2<<endl;
	TH2D *hMCAcc1MvsPtphidalitz = (TH2D*)fphidal->Get("hMCAcc1MvsPtphidalitz");
	hMCAcc1MvsPtphidalitz->Sumw2();
	TH2D *hRCAcc1MvsPt3Dphidalitz = (TH2D*)fphidal->Get("hRCAcc1MvsPt3Dphidalitz");
	hRCAcc1MvsPt3Dphidalitz->Sumw2();
	TH2D *hUpMCAcc1MvsPtphidalitz = (TH2D *)hMCAcc1MvsPtphidalitz->Clone("hUpMCAcc1MvsPtphidalitz");
	TH2D *hLowMCAcc1MvsPtphidalitz = (TH2D *)hMCAcc1MvsPtphidalitz->Clone("hLowMCAcc1MvsPtphidalitz");
	cout<<1<<endl;

	//two body decay
	//TH2D *hMCMvsPtrho2ee = (TH2D*)frho2ee->Get("hMCrhoeehist");
	//hMCMvsPtrho2ee->Sumw2();
	//TH2D *hMCAcc0MvsPtrho2ee = (TH2D*)frho2ee->Get("hMCAcc0rhoeehist");
	//hMCAcc0MvsPtrho2ee->Sumw2();
	//TH2D *hMCAcc1MvsPtrho2ee = (TH2D*)frho2ee->Get("hMCAcc1rhoeehist");
	//hMCAcc1MvsPtrho2ee->Sumw2();
	//TH2D *hRCAcc1MvsPt3Drho2ee = (TH2D*)frho2ee->Get("hRCrhoeehist");
	//hRCAcc1MvsPt3Drho2ee->Sumw2();
	//TH2D *hUpMCAcc1MvsPtrho2ee = (TH2D *)hMCAcc1MvsPtrho2ee->Clone("hUpMCAcc1MvsPtrho2ee");
	//TH2D *hLowMCAcc1MvsPtrho2ee = (TH2D *)hMCAcc1MvsPtrho2ee->Clone("hLowMCAcc1MvsPtrho2ee");

	TH2D *hMCMvsPtomega2ee = (TH2D*)fomega2ee->Get("hMCMvsPtomega2ee");
	hMCMvsPtomega2ee->Sumw2();
	TH2D *hMCAcc0MvsPtomega2ee = (TH2D*)fomega2ee->Get("hMCAcc0MvsPtomega2ee");
	hMCAcc0MvsPtomega2ee->Sumw2();
	TH2D *hMCAcc1MvsPtomega2ee = (TH2D*)fomega2ee->Get("hMCAcc1MvsPtomega2ee");
	hMCAcc1MvsPtomega2ee->Sumw2();
	TH2D *hRCAcc1MvsPt3Domega2ee = (TH2D*)fomega2ee->Get("hRCAcc1MvsPt3Domega2ee");
	hRCAcc1MvsPt3Domega2ee->Sumw2();
	TH2D *hUpMCAcc1MvsPtomega2ee = (TH2D *)hMCAcc1MvsPtomega2ee->Clone("hUpMCAcc1MvsPtomega2ee");
	TH2D *hLowMCAcc1MvsPtomega2ee = (TH2D *)hMCAcc1MvsPtomega2ee->Clone("hLowMCAcc1MvsPtomega2ee");

	TH2D *hMCMvsPtphi2ee = (TH2D*)fphi2ee->Get("hMCMvsPtphi2ee");
	hMCMvsPtphi2ee->Sumw2();
	TH2D *hMCAcc0MvsPtphi2ee = (TH2D*)fphi2ee->Get("hMCAcc0MvsPtphi2ee");
	hMCAcc0MvsPtphi2ee->Sumw2();
	TH2D *hMCAcc1MvsPtphi2ee = (TH2D*)fphi2ee->Get("hMCAcc1MvsPtphi2ee");
	hMCAcc1MvsPtphi2ee->Sumw2();
	TH2D *hRCAcc1MvsPt3Dphi2ee = (TH2D*)fphi2ee->Get("hRCAcc1MvsPt3Dphi2ee");
	hRCAcc1MvsPt3Dphi2ee->Sumw2();
	TH2D *hUpMCAcc1MvsPtphi2ee = (TH2D *)hMCAcc1MvsPtphi2ee->Clone("hUpMCAcc1MvsPtphi2ee");
	TH2D *hLowMCAcc1MvsPtphi2ee = (TH2D *)hMCAcc1MvsPtphi2ee->Clone("hLowMCAcc1MvsPtphi2ee");

	TH2D *hMCMvsPtjpsi2ee = (TH2D*)fjpsi2ee->Get("hMCMvsPtjpsi2ee");
	hMCMvsPtjpsi2ee->Sumw2();
	TH2D *hMCAcc0MvsPtjpsi2ee = (TH2D*)fjpsi2ee->Get("hMCAcc0MvsPtjpsi2ee");
	hMCAcc0MvsPtjpsi2ee->Sumw2();
	TH2D *hMCAcc1MvsPtjpsi2ee = (TH2D*)fjpsi2ee->Get("hMCAcc1MvsPtjpsi2ee");
	hMCAcc1MvsPtjpsi2ee->Sumw2();
	TH2D *hRCAcc1MvsPt3Djpsi2ee = (TH2D*)fjpsi2ee->Get("hRCAcc1MvsPt3Djpsi2ee");
	hRCAcc1MvsPt3Djpsi2ee->Sumw2();
	TH2D *hUpMCAcc1MvsPtjpsi2ee = (TH2D *)hMCAcc1MvsPtjpsi2ee->Clone("hUpMCAcc1MvsPtjpsi2ee");
	TH2D *hLowMCAcc1MvsPtjpsi2ee = (TH2D *)hMCAcc1MvsPtjpsi2ee->Clone("hLowMCAcc1MvsPtjpsi2ee");

	//TH2D *hMCMvsPtpsi2ee = (TH2D*)fpsi2ee->Get("hMCMvsPtpsi2ee");
	//hMCMvsPtpsi2ee->Sumw2();
	//TH2D *hMCAcc0MvsPtpsi2ee = (TH2D*)fpsi2ee->Get("hMCAcc0MvsPtpsi2ee");
	//hMCAcc0MvsPtpsi2ee->Sumw2();
	//TH2D *hMCAcc1MvsPtpsi2ee = (TH2D*)fpsi2ee->Get("hMCAcc1MvsPtpsi2ee");
	//hMCAcc1MvsPtpsi2ee->Sumw2();
	//TH2D *hRCAcc1MvsPt3Dpsi2ee = (TH2D*)fpsi2ee->Get("hRCAcc1MvsPt3Dpsi2ee");
	//hRCAcc1MvsPt3Dpsi2ee->Sumw2();
	//TH2D *hUpMCAcc1MvsPtpsi2ee = (TH2D *)hMCAcc1MvsPtpsi2ee->Clone("hUpMCAcc1MvsPtpsi2ee");
	//TH2D *hLowMCAcc1MvsPtpsi2ee = (TH2D *)hMCAcc1MvsPtpsi2ee->Clone("hLowMCAcc1MvsPtpsi2ee");

	//heavy flavor decay
	TH1D *hMCnCcEvts = (TH1D *)fcc->Get("hnEvts"); // 1-inclusive c statistics; 2-two charm string statistics
	TH2D *hMCFullYMvsPtcc2ee = (TH2D *)fcc->Get("hMCFullYMvsPtcc2ee");
	hMCFullYMvsPtcc2ee->SetNameTitle("hMCFullYMvsPtcc2ee","hMCFullYMvsPtcc2ee");
	hMCFullYMvsPtcc2ee->Sumw2();
	TH2D *hMCAcc0MvsPtcc2ee = (TH2D*)fcc->Get("hMCAcc0MvsPtcc2ee");
	hMCAcc0MvsPtcc2ee->SetNameTitle("hMCAcc0MvsPtcc2ee","hMCAcc0MvsPtcc2ee");
	hMCAcc0MvsPtcc2ee->Sumw2();
	TH2D *hMCAcc1MvsPtcc2ee = (TH2D*)fcc->Get("hMCAcc1MvsPtcc2ee");
	hMCAcc1MvsPtcc2ee->SetNameTitle("hMCAcc1MvsPtcc2ee","hMCAcc1MvsPtcc2ee");
	hMCAcc1MvsPtcc2ee->Sumw2();
	TH2D *hRCAcc1MvsPt3Dcc2ee = (TH2D*)fcc->Get("hRCAcc1MvsPt3Dcc2ee");
	hRCAcc1MvsPt3Dcc2ee->SetNameTitle("hRCAcc1MvsPt3Dcc2ee","hRCAcc1MvsPt3Dcc2ee");
	hRCAcc1MvsPt3Dcc2ee->Sumw2();
	TH2D *hUpMCAcc1MvsPtcc2ee = (TH2D *)hMCAcc1MvsPtcc2ee->Clone("hUpMCAcc1MvsPtcc2ee");
	TH2D *hLowMCAcc1MvsPtcc2ee = (TH2D *)hMCAcc1MvsPtcc2ee->Clone("hLowMCAcc1MvsPtcc2ee");

	//TH1D *hMCMvsPtbbnParent = (TH1D*)fbb->Get("nParent");
	//TH2D *hMCFullYMvsPtbb2ee = (TH2D *)fbb->Get("hMCFullYMvsPtbb2ee");
	//hMCFullYMvsPtbb2ee->SetNameTitle("hMCFullYMvsPtbb2ee","hMCFullYMvsPtbb2ee");
	//hMCFullYMvsPtbb2ee->Sumw2();
	//TH2D *hMCAcc0MvsPtbb2ee = (TH2D*)fbb->Get("hMCAcc0MvsPtbb2ee");
	//hMCAcc0MvsPtbb2ee->SetNameTitle("hMCAcc0MvsPtbb2ee","hMCAcc0MvsPtbb2ee");
	//hMCAcc0MvsPtbb2ee->Sumw2();
	//TH2D *hMCAcc1MvsPtbb2ee = (TH2D*)fbb->Get("hMCAcc1MvsPtbb2ee");
	//hMCAcc1MvsPtbb2ee->SetNameTitle("hMCAcc1MvsPtbb2ee","hMCAcc1MvsPtbb2ee");
	//hMCAcc1MvsPtbb2ee->Sumw2();
	//TH2D *hRCAcc1MvsPt3Dbb2ee = (TH2D*)fbb->Get("hRCAcc1MvsPt3Dbb2ee");
	//hRCAcc1MvsPt3Dbb2ee->SetNameTitle("hRCAcc1MvsPt3Dbb2ee","hRCAcc1MvsPt3Dbb2ee");
	//hRCAcc1MvsPt3Dbb2ee->Sumw2();
	//TH2D *hUpMCAcc1MvsPtbb2ee = (TH2D *)hMCAcc1MvsPtbb2ee->Clone("hUpMCAcc1MvsPtbb2ee");
	//TH2D *hLowMCAcc1MvsPtbb2ee = (TH2D *)hMCAcc1MvsPtbb2ee->Clone("hLowMCAcc1MvsPtbb2ee");

	//TH1D *hMCMvsPtdynParent = (TH1D*)fdy->Get("nParent");
	//TH2D *hMCFullYMvsPtdy2ee = (TH2D *)fdy->Get("hMCFullYMvsPtdy2ee");
	//hMCFullYMvsPtdy2ee->SetNameTitle("hMCFullYMvsPtdy2ee","hMCFullYMvsPtdy2ee");
	//hMCFullYMvsPtdy2ee->Sumw2();
	//TH2D *hMCAcc0MvsPtdy2ee = (TH2D*)fdy->Get("hMCAcc0MvsPtdy2ee");
	//hMCAcc0MvsPtdy2ee->SetNameTitle("hMCAcc0MvsPtdy2ee","hMCAcc0MvsPtdy2ee");
	//hMCAcc0MvsPtdy2ee->Sumw2();
	//TH2D *hMCAcc1MvsPtdy2ee = (TH2D*)fdy->Get("hMCAcc1MvsPtdy2ee");
	//hMCAcc1MvsPtdy2ee->SetNameTitle("hMCAcc1MvsPtdy2ee","hMCAcc1MvsPtdy2ee");
	//hMCAcc1MvsPtdy2ee->Sumw2();
	//TH2D *hRCAcc1MvsPt3Ddy2ee = (TH2D*)fdy->Get("hRCAcc1MvsPt3Ddy2ee");
	//hRCAcc1MvsPt3Ddy2ee->SetNameTitle("hRCAcc1MvsPt3Ddy2ee","hRCAcc1MvsPt3Ddy2ee");
	//hRCAcc1MvsPt3Ddy2ee->Sumw2();
	//TH2D *hUpMCAcc1MvsPtdy2ee = (TH2D *)hMCAcc1MvsPtdy2ee->Clone("hUpMCAcc1MvsPtdy2ee");
	//TH2D *hLowMCAcc1MvsPtdy2ee = (TH2D *)hMCAcc1MvsPtdy2ee->Clone("hLowMCAcc1MvsPtdy2ee");

	const Double_t npi0dal = hMCMvsPtpi0dalitz->GetEntries();
	const Double_t netadal = hMCMvsPtetadalitz->GetEntries();
	const Double_t nomegadal = hMCMvsPtomegadalitz->GetEntries();
	const Double_t netaprimdal = hMCMvsPtetaprimdalitz->GetEntries();
	const Double_t nphidal = hMCMvsPtphidalitz->GetEntries();

	//const Double_t nrho2ee = hMCMvsPtrho2ee->GetEntries();
	const Double_t nomega2ee = hMCMvsPtomega2ee->GetEntries();
	const Double_t nphi2ee = hMCMvsPtphi2ee->GetEntries();
	const Double_t njpsi2ee = hMCMvsPtjpsi2ee->GetEntries();
	//const Double_t npsi2ee = hMCMvsPtpsi2ee->GetEntries();

	const Double_t nInclusiveCharm = hMCnCcEvts->GetBinContent(1); // inclusive charm statistics
	//const Double_t ncc = hMCnCcEvts->GetBinContent(2); // two charm string statistics
	//const Double_t nbb = hMCMvsPtbbnParent->GetBinContent(1);
	//const Double_t nbb = hMCFullYMvsPtbb2ee->GetEntries(); //here just keep consistency with AuAu@200GeV PRC paper
	//const Double_t ndy = hMCMvsPtdynParent->GetBinContent(1);


	//Input ee pair |Y|<1
	hMCAcc0MvsPtpi0dalitz->Scale(dNdy_pi0*BR_pi0dal/npi0dal*dY/dM);
	hMCAcc0MvsPtetadalitz->Scale(dNdy_eta*BR_etadal/netadal*dY/dM);
	hMCAcc0MvsPtomegadalitz->Scale(dNdy_omega*BR_omegadal/nomegadal*dY/dM);
	hMCAcc0MvsPtetaprimdalitz->Scale(dNdy_etaprim*BR_etaprimdal/netaprimdal*dY/dM);
	hMCAcc0MvsPtphidalitz->Scale(dNdy_phi*BR_phidal/nphidal*dY/dM);
	//hMCAcc0MvsPtrho2ee->Scale(dNdy_rho*BR_rho2ee/nrho2ee*dY/dM);
	hMCAcc0MvsPtomega2ee->Scale(dNdy_omega*BR_omega2ee/nomega2ee*dY/dM);
	hMCAcc0MvsPtphi2ee->Scale(dNdy_phi*BR_phi2ee/nphi2ee*dY/dM);
	hMCAcc0MvsPtjpsi2ee->Scale(dNdy_jpsi*BR_jpsi2ee/njpsi2ee*dY/dM);
	//hMCAcc0MvsPtpsi2ee->Scale(dNdy_psi*BR_psi2ee/npsi2ee*dY/dM);
	hMCAcc0MvsPtcc2ee->Scale(dNtotcc/nInclusiveCharm/dM);
	//hMCAcc0MvsPtcc2ee->Scale(dNtotcc/ncc/dM/1.43); //1.43 is NInclusiveCharmEvents/NTwoCharmStringEvents
	//hMCAcc0MvsPtbb2ee->Scale(dNtotbb/nbb/dM);//note here is no cc->ee or bb->ee branch ratio
	//hMCAcc0MvsPtdy2ee->Scale(dNtotdy/ndy/dM);

	TH2D *hMCAcc0MvsPt = (TH2D *)hMCAcc0MvsPtpi0dalitz->Clone("hMCAcc0MvsPt");
	hMCAcc0MvsPt->SetTitle("hMCAcc0MvsPt");
	hMCAcc0MvsPt->Add(hMCAcc0MvsPtetadalitz);
	hMCAcc0MvsPt->Add(hMCAcc0MvsPtomegadalitz);
	hMCAcc0MvsPt->Add(hMCAcc0MvsPtetaprimdalitz);
	hMCAcc0MvsPt->Add(hMCAcc0MvsPtphidalitz);
	//if(withRho) hMCAcc0MvsPt->Add(hMCAcc0MvsPtrho2ee);
	hMCAcc0MvsPt->Add(hMCAcc0MvsPtomega2ee);
	hMCAcc0MvsPt->Add(hMCAcc0MvsPtphi2ee);
	hMCAcc0MvsPt->Add(hMCAcc0MvsPtjpsi2ee);
	//hMCAcc0MvsPt->Add(hMCAcc0MvsPtpsi2ee);
	hMCAcc0MvsPt->Add(hMCAcc0MvsPtcc2ee);
	//hMCAcc0MvsPt->Add(hMCAcc0MvsPtbb2ee);
	//hMCAcc0MvsPt->Add(hMCAcc0MvsPtdy2ee);

	TH2D *hMCAcc0MvsPtomega = (TH2D *)hMCAcc0MvsPtomegadalitz->Clone("hMCAcc0MvsPtomega");
	hMCAcc0MvsPtomega->SetTitle("hMCAcc0MvsPtomega");
	hMCAcc0MvsPtomega->Add(hMCAcc0MvsPtomega2ee);

	TH2D *hMCAcc0MvsPtphi = (TH2D *)hMCAcc0MvsPtphidalitz->Clone("hMCAcc0MvsPtphi");
	hMCAcc0MvsPtphi->SetTitle("hMCAcc0MvsPtphi");
	hMCAcc0MvsPtphi->Add(hMCAcc0MvsPtphi2ee);

	//Input ee pair |Y|<1 && pT_{e}>0.2 && |eta_{e}|<1
	hMCAcc1MvsPtpi0dalitz->Scale(dNdy_pi0*BR_pi0dal/npi0dal*dY/dM);
	hMCAcc1MvsPtetadalitz->Scale(dNdy_eta*BR_etadal/netadal*dY/dM);
	hMCAcc1MvsPtomegadalitz->Scale(dNdy_omega*BR_omegadal/nomegadal*dY/dM);
	hMCAcc1MvsPtetaprimdalitz->Scale(dNdy_etaprim*BR_etaprimdal/netaprimdal*dY/dM);
	hMCAcc1MvsPtphidalitz->Scale(dNdy_phi*BR_phidal/nphidal*dY/dM);
	//hMCAcc1MvsPtrho2ee->Scale(dNdy_rho*BR_rho2ee/nrho2ee*dY/dM);
	hMCAcc1MvsPtomega2ee->Scale(dNdy_omega*BR_omega2ee/nomega2ee*dY/dM);
	hMCAcc1MvsPtphi2ee->Scale(dNdy_phi*BR_phi2ee/nphi2ee*dY/dM);
	hMCAcc1MvsPtjpsi2ee->Scale(dNdy_jpsi*BR_jpsi2ee/njpsi2ee*dY/dM);
	//hMCAcc1MvsPtpsi2ee->Scale(dNdy_psi*BR_psi2ee/npsi2ee*dY/dM);
	hMCAcc1MvsPtcc2ee->Scale(dNtotcc/nInclusiveCharm/dM);
	//hMCAcc1MvsPtcc2ee->Scale(dNtotcc/ncc/dM/1.43); //1.43 is NInclusiveCharmEvents/NTwoCharmStringEvents
	//hMCAcc1MvsPtbb2ee->Scale(dNtotbb/nbb/dM);//note here is no cc->ee or bb->ee branch ratio
	//hMCAcc1MvsPtdy2ee->Scale(dNtotdy/ndy/dM);

	TH2D *hMCAcc1MvsPt = (TH2D *)hMCAcc1MvsPtpi0dalitz->Clone("hMCAcc1MvsPt");
	hMCAcc1MvsPt->SetTitle("hMCAcc1MvsPt");
	hMCAcc1MvsPt->Add(hMCAcc1MvsPtetadalitz);
	hMCAcc1MvsPt->Add(hMCAcc1MvsPtomegadalitz);
	hMCAcc1MvsPt->Add(hMCAcc1MvsPtetaprimdalitz);
	hMCAcc1MvsPt->Add(hMCAcc1MvsPtphidalitz);
	//if(withRho) hMCAcc1MvsPt->Add(hMCAcc1MvsPtrho2ee);
	hMCAcc1MvsPt->Add(hMCAcc1MvsPtomega2ee);
	hMCAcc1MvsPt->Add(hMCAcc1MvsPtphi2ee);
	hMCAcc1MvsPt->Add(hMCAcc1MvsPtjpsi2ee);
	//hMCAcc1MvsPt->Add(hMCAcc1MvsPtpsi2ee);
	hMCAcc1MvsPt->Add(hMCAcc1MvsPtcc2ee);
	//hMCAcc1MvsPt->Add(hMCAcc1MvsPtbb2ee);
	//hMCAcc1MvsPt->Add(hMCAcc1MvsPtdy2ee);

	TH2D *hMCAcc1MvsPtomega = (TH2D *)hMCAcc1MvsPtomegadalitz->Clone("hMCAcc1MvsPtomega");
	hMCAcc1MvsPtomega->SetTitle("hMCAcc1MvsPtomega");
	hMCAcc1MvsPtomega->Add(hMCAcc1MvsPtomega2ee);

	TH2D *hMCAcc1MvsPtphi = (TH2D *)hMCAcc1MvsPtphidalitz->Clone("hMCAcc1MvsPtphi");
	hMCAcc1MvsPtphi->SetTitle("hMCAcc1MvsPtphi");
	hMCAcc1MvsPtphi->Add(hMCAcc1MvsPtphi2ee);

	//Cocktail up limit
	hUpMCAcc1MvsPtpi0dalitz->Scale((1+dNdy_pi0Err)*dNdy_pi0*BR_pi0dal/npi0dal*dY/dM);
	hUpMCAcc1MvsPtetadalitz->Scale((1+dNdy_etaErr)*dNdy_eta*BR_etadal/netadal*dY/dM);
	hUpMCAcc1MvsPtomegadalitz->Scale((1+dNdy_omegaErr)*dNdy_omega*BR_omegadal/nomegadal*dY/dM);
	hUpMCAcc1MvsPtetaprimdalitz->Scale((1+dNdy_etaprimErr)*dNdy_etaprim*BR_etaprimdal/netaprimdal*dY/dM);
	hUpMCAcc1MvsPtphidalitz->Scale((1+dNdy_phiErr)*dNdy_phi*BR_phidal/nphidal*dY/dM);
	//hUpMCAcc1MvsPtrho2ee->Scale((1+dNdy_rhoErr)*dNdy_rho*BR_rho2ee/nrho2ee*dY/dM);
	hUpMCAcc1MvsPtomega2ee->Scale((1+dNdy_omegaErr)*dNdy_omega*BR_omega2ee/nomega2ee*dY/dM);
	hUpMCAcc1MvsPtphi2ee->Scale((1+dNdy_phiErr)*dNdy_phi*BR_phi2ee/nphi2ee*dY/dM);
	hUpMCAcc1MvsPtjpsi2ee->Scale((1+dNdy_jpsiErr)*dNdy_jpsi*BR_jpsi2ee/njpsi2ee*dY/dM);
	//hUpMCAcc1MvsPtpsi2ee->Scale((1+dNdy_psiErr)*dNdy_psi*BR_psi2ee/npsi2ee*dY/dM);
	hUpMCAcc1MvsPtcc2ee->Scale((1+dNdy_ccErr)*dNtotcc/nInclusiveCharm/dM);
	//hUpMCAcc1MvsPtcc2ee->Scale((1+dNdy_ccErr)*dNtotcc/ncc/dM/1.43); //1.43 is NInclusiveCharmEvents/NTwoCharmStringEvents
	//hUpMCAcc1MvsPtbb2ee->Scale((1+dNdy_bbErr)*dNtotbb/nbb/dM);//note here is no cc->ee or bb->ee branch ratio
	//hUpMCAcc1MvsPtdy2ee->Scale((1+dNdy_dyErr)*dNtotdy/ndy/dM);

	TH2D *hUpMCAcc1MvsPt = (TH2D *)hUpMCAcc1MvsPtpi0dalitz->Clone("hUpMCAcc1MvsPt");
	hUpMCAcc1MvsPt->SetTitle("hUpMCAcc1MvsPt");
	hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtetadalitz);
	hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtomegadalitz);
	hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtetaprimdalitz);
	hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtphidalitz);
	//if(withRho) hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtrho2ee);
	hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtomega2ee);
	hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtphi2ee);
	hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtjpsi2ee);
	//hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtpsi2ee);
	hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtcc2ee);
	//hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtbb2ee);
	//hUpMCAcc1MvsPt->Add(hUpMCAcc1MvsPtdy2ee);

	//Cocktail low limit
	hLowMCAcc1MvsPtpi0dalitz->Scale((1-dNdy_pi0Err)*dNdy_pi0*BR_pi0dal/npi0dal*dY/dM);
	hLowMCAcc1MvsPtetadalitz->Scale((1-dNdy_etaErr)*dNdy_eta*BR_etadal/netadal*dY/dM);
	hLowMCAcc1MvsPtomegadalitz->Scale((1-dNdy_omegaErr)*dNdy_omega*BR_omegadal/nomegadal*dY/dM);
	hLowMCAcc1MvsPtetaprimdalitz->Scale((1-dNdy_etaprimErr)*dNdy_etaprim*BR_etaprimdal/netaprimdal*dY/dM);
	hLowMCAcc1MvsPtphidalitz->Scale((1-dNdy_phiErr)*dNdy_phi*BR_phidal/nphidal*dY/dM);
	//hLowMCAcc1MvsPtrho2ee->Scale((1-dNdy_rhoErr)*dNdy_rho*BR_rho2ee/nrho2ee*dY/dM);
	hLowMCAcc1MvsPtomega2ee->Scale((1-dNdy_omegaErr)*dNdy_omega*BR_omega2ee/nomega2ee*dY/dM);
	hLowMCAcc1MvsPtphi2ee->Scale((1-dNdy_phiErr)*dNdy_phi*BR_phi2ee/nphi2ee*dY/dM);
	hLowMCAcc1MvsPtjpsi2ee->Scale((1-dNdy_jpsiErr)*dNdy_jpsi*BR_jpsi2ee/njpsi2ee*dY/dM);
	//hLowMCAcc1MvsPtpsi2ee->Scale((1-dNdy_psiErr)*dNdy_psi*BR_psi2ee/npsi2ee*dY/dM);
	hLowMCAcc1MvsPtcc2ee->Scale((1-dNdy_ccErr)*dNtotcc/nInclusiveCharm/dM);
	//hLowMCAcc1MvsPtcc2ee->Scale((1-dNdy_ccErr)*dNtotcc/ncc/dM/1.43); //1.43 is NInclusiveCharmEvents/NTwoCharmStringEvents
	//hLowMCAcc1MvsPtbb2ee->Scale((1-dNdy_bbErr)*dNtotbb/nbb/dM);//note here is no cc->ee or bb->ee branch ratio
	//hLowMCAcc1MvsPtdy2ee->Scale((1-dNdy_dyErr)*dNtotdy/ndy/dM);

	TH2D *hLowMCAcc1MvsPt = (TH2D *)hLowMCAcc1MvsPtpi0dalitz->Clone("hLowMCAcc1MvsPt");
	hLowMCAcc1MvsPt->SetTitle("hLowMCAcc1MvsPt");
	hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtetadalitz);
	hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtomegadalitz);
	hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtetaprimdalitz);
	hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtphidalitz);
	//if(withRho) hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtrho2ee);
	hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtomega2ee);
	hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtphi2ee);
	hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtjpsi2ee);
	//hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtpsi2ee);
	hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtcc2ee);
	//hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtbb2ee);
	//hLowMCAcc1MvsPt->Add(hLowMCAcc1MvsPtdy2ee);

	//Reconstruct ee pair
	hRCAcc1MvsPt3Dpi0dalitz->Scale(dNdy_pi0*BR_pi0dal/npi0dal*dY/dM);
	hRCAcc1MvsPt3Detadalitz->Scale(dNdy_eta*BR_etadal/netadal*dY/dM);
	hRCAcc1MvsPt3Domegadalitz->Scale(dNdy_omega*BR_omegadal/nomegadal*dY/dM);
	hRCAcc1MvsPt3Detaprimdalitz->Scale(dNdy_etaprim*BR_etaprimdal/netaprimdal*dY/dM);
	hRCAcc1MvsPt3Dphidalitz->Scale(dNdy_phi*BR_phidal/nphidal*dY/dM);
	//hRCAcc1MvsPt3Drho2ee->Scale(dNdy_rho*BR_rho2ee/nrho2ee*dY/dM);
	hRCAcc1MvsPt3Domega2ee->Scale(dNdy_omega*BR_omega2ee/nomega2ee*dY/dM);
	hRCAcc1MvsPt3Dphi2ee->Scale(dNdy_phi*BR_phi2ee/nphi2ee*dY/dM);
	hRCAcc1MvsPt3Djpsi2ee->Scale(dNdy_jpsi*BR_jpsi2ee/njpsi2ee*dY/dM);
	//hRCAcc1MvsPt3Dpsi2ee->Scale(dNdy_psi*BR_psi2ee/npsi2ee*dY/dM);
	hRCAcc1MvsPt3Dcc2ee->Scale(dNtotcc/nInclusiveCharm/dM);
	//hRCAcc1MvsPt3Dcc2ee->Scale(dNtotcc/ncc/dM/1.43); //1.43 is NInclusiveCharmEvents/NTwoCharmStringEvents
	//hRCAcc1MvsPt3Dbb2ee->Scale(dNtotbb/nbb/dM);//note here is no cc->ee or bb->ee branch ratio
	//hRCAcc1MvsPt3Ddy2ee->Scale(dNtotdy/ndy/dM);

	TH2D *hRCAcc1MvsPt3D = (TH2D *)hRCAcc1MvsPt3Dpi0dalitz->Clone("hRCAcc1MvsPt3D");
	hRCAcc1MvsPt3D->SetTitle("hRCAcc1MvsPt3D");
	hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Detadalitz);
	hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Domegadalitz);
	hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Detaprimdalitz);
	hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Dphidalitz);
	//if(withRho) hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Drho2ee);
	hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Domega2ee);
	hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Dphi2ee);
	hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Djpsi2ee);
	//hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Dpsi2ee);
	hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Dcc2ee);
	//hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Dbb2ee);
	//hRCAcc1MvsPt3D->Add(hRCAcc1MvsPt3Ddy2ee);

	TH2D *hRCAcc1MvsPt3Domega = (TH2D *)hRCAcc1MvsPt3Domegadalitz->Clone("hRCAcc1MvsPt3Domega");
	hRCAcc1MvsPt3Domega->SetTitle("hRCAcc1MvsPt3Domega");
	hRCAcc1MvsPt3Domega->Add(hRCAcc1MvsPt3Domega2ee);

	TH2D *hRCAcc1MvsPt3Dphi = (TH2D *)hRCAcc1MvsPt3Dphidalitz->Clone("hRCAcc1MvsPt3Dphi");
	hRCAcc1MvsPt3Dphi->SetTitle("hRCAcc1MvsPt3Dphi");
	hRCAcc1MvsPt3Dphi->Add(hRCAcc1MvsPt3Dphi2ee);


	TLegend *leg = new TLegend(0.2,0.65,0.8,0.85);
	leg->SetBorderSize(0);
	leg->SetFillColor(10);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);

	Int_t binLow,binHi;

	hMCAcc0MvsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hMCAcc0MvsPt->GetYaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	binLow = hMCAcc0MvsPt->GetXaxis()->FindBin(PtMin+1.e-6);
	binHi = hMCAcc0MvsPt->GetXaxis()->FindBin(PtMax-1.e-6);
	TH1D *hMCAcc0M = (TH1D *)hMCAcc0MvsPt->ProjectionY("hMCAcc0M",binLow,binHi);
	hMCAcc0M->SetTitle(Form("hMCAcc0M_Pt_%3.1f_%3.1f",PtMin,PtMax));
	hMCAcc0M->GetYaxis()->SetTitle("dN/dM");

	hMCAcc1MvsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hMCAcc1MvsPt->GetYaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	TH1D *hMCAcc1M = (TH1D *)hMCAcc1MvsPt->ProjectionY("hMCAcc1M",binLow,binHi);
	hMCAcc1M->SetTitle(Form("hMCAcc1M_Pt_%3.1f_%3.1f",PtMin,PtMax));
	hMCAcc1M->GetYaxis()->SetTitle("dN/dM");

	hUpMCAcc1MvsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hUpMCAcc1MvsPt->GetYaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	TH1D *hUpMCAcc1M = (TH1D *)hUpMCAcc1MvsPt->ProjectionY("hUpMCAcc1M",binLow,binHi);
	hUpMCAcc1M->SetTitle(Form("hUpMCAcc1M_Pt_%3.1f_%3.1f",PtMin,PtMax));
	hUpMCAcc1M->GetYaxis()->SetTitle("dN/dM");

	hLowMCAcc1MvsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hLowMCAcc1MvsPt->GetYaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	TH1D *hLowMCAcc1M = (TH1D *)hLowMCAcc1MvsPt->ProjectionY("hLowMCAcc1M",binLow,binHi);
	hLowMCAcc1M->SetTitle(Form("hLowMCAcc1M_Pt_%3.1f_%3.1f",PtMin,PtMax));
	hLowMCAcc1M->GetYaxis()->SetTitle("dN/dM");

	hRCAcc1MvsPt3D->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hRCAcc1MvsPt3D->GetYaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	TH1D *hRCAcc1M = (TH1D *)hRCAcc1MvsPt3D->ProjectionY("hRCAcc1M",binLow,binHi);
	hRCAcc1M->SetTitle(Form("hRCAcc1M_Pt_%3.1f_%3.1f",PtMin,PtMax));
	hRCAcc1M->GetYaxis()->SetTitle("dN/dM");

	c1->cd();
	gPad->SetLogy(1);
	hMCAcc0M->SetLineColor(1);
	hMCAcc0M->GetXaxis()->SetRangeUser(0,4.);
	hMCAcc0M->GetXaxis()->SetTitle("M_{ee} (Gev/c^{2})");
	hMCAcc0M->GetXaxis()->SetLabelFont(22);
	hMCAcc0M->GetXaxis()->SetLabelSize(0.04);
	hMCAcc0M->GetXaxis()->SetTitleFont(22);
	hMCAcc0M->GetXaxis()->SetTitleSize(0.05);
	hMCAcc0M->GetYaxis()->SetTitle("dN/dM (c^{2}/GeV)");
	hMCAcc0M->SetMinimum(1.e-7);
	hMCAcc1M->SetLineColor(2);
	hRCAcc1M->SetLineColor(4);
	hMCAcc0M->Draw("chist");
	hMCAcc1M->Draw("chistsame");
	hRCAcc1M->Draw("chistsame");
	leg->AddEntry(hMCAcc0M,"cocktail (|Y_{ee}|<1)","pl");
	leg->AddEntry(hMCAcc1M,"cocktail (|Y_{ee}|<1, p_{T}^{e}>0.2 GeV/c, |#eta_{e}|<1)","pl");
	leg->AddEntry(hRCAcc1M,"cocktail (|Y_{ee}|<1, p_{T}^{e}>0.2 GeV/c, |#eta_{e}|<1, w/ eff.)","pl");
	leg->Draw("same");
	drawLatex(0.18,0.2,"Au+Au #sqrt{s_{NN}}=62 GeV (MinBias)",22,0.05,1);
	pdfAction(c1,ps);
	if(withRho){
		c1->SaveAs(Form("%s/cen%d%d_cocktail_withRho.png", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
		c1->SaveAs(Form("%s/cen%d%d_cocktail_withRho.pdf", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
		c1->SaveAs(Form("%s/cen%d%d_cocktail_withRho.eps", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	}
	else {
		c1->SaveAs(Form("%s/cen%d%d_cocktail_withoutRho.png", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
		c1->SaveAs(Form("%s/cen%d%d_cocktail_withoutRho.pdf", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
		c1->SaveAs(Form("%s/cen%d%d_cocktail_withoutRho.eps", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]));
	}
	hMCAcc0M->GetXaxis()->UnZoom();

	Int_t nColumns = 2;
	Int_t nRaws = 2;
	Int_t nPads = nColumns*nRaws;

	c1->Clear();
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);

	const Int_t nPtBins = 7;  
	Double_t PtLow[nPtBins] = {0.,   0.15, 1, 0,   0.5, 1,   1.5};
	Double_t PtHi[nPtBins]  = {0.15, 1,   10, 0.5, 1,   1.5, 5.0};

	TH1D *hMCAcc0M_Pt[nPtBins];
	TH1D *hMCAcc1M_Pt[nPtBins];
	TH1D *hUpMCAcc1M_Pt[nPtBins];
	TH1D *hLowMCAcc1M_Pt[nPtBins];
	TH1D *hRCAcc1M_Pt[nPtBins];
	for(Int_t i=0;i<nPtBins;i++){
		binLow = hMCAcc0MvsPt->GetXaxis()->FindBin(PtLow[i]+1.e-6);
		binHi = hMCAcc0MvsPt->GetXaxis()->FindBin(PtHi[i]-1.e-6);
		hMCAcc0M_Pt[i] = (TH1D *)hMCAcc0MvsPt->ProjectionY(Form("hMCAcc0M_PtBin%d",i),binLow,binHi);
		hMCAcc0M_Pt[i]->SetTitle(Form("hMCAcc0M_Pt_%3.1f_%3.1f",PtLow[i],PtHi[i]));
		hMCAcc1M_Pt[i] = (TH1D *)hMCAcc1MvsPt->ProjectionY(Form("hMCAcc1M_PtBin%d",i),binLow,binHi);
		hMCAcc1M_Pt[i]->SetTitle(Form("hMCAcc1M_Pt_%3.1f_%3.1f",PtLow[i],PtHi[i]));
		hUpMCAcc1M_Pt[i] = (TH1D *)hUpMCAcc1MvsPt->ProjectionY(Form("hUpMCAcc1M_PtBin%d",i),binLow,binHi);
		hUpMCAcc1M_Pt[i]->SetTitle(Form("hUpMCAcc1M_Pt_%3.1f_%3.1f",PtLow[i],PtHi[i]));
		hLowMCAcc1M_Pt[i] = (TH1D *)hLowMCAcc1MvsPt->ProjectionY(Form("hLowMCAcc1M_PtBin%d",i),binLow,binHi);
		hLowMCAcc1M_Pt[i]->SetTitle(Form("hLowMCAcc1M_Pt_%3.1f_%3.1f",PtLow[i],PtHi[i]));
		hRCAcc1M_Pt[i] = (TH1D *)hRCAcc1MvsPt3D->ProjectionY(Form("hRCAcc1M_PtBin%d",i),binLow,binHi);
		hRCAcc1M_Pt[i]->SetTitle(Form("hRCAcc1M_Pt_%3.1f_%3.1f",PtLow[i],PtHi[i]));
		hMCAcc0M_Pt[i]->SetLineColor(1);
		hMCAcc0M_Pt[i]->GetXaxis()->SetRangeUser(0,4.);
		hMCAcc0M_Pt[i]->GetYaxis()->SetTitle("dN/dM");
		hMCAcc0M_Pt[i]->SetMinimum(1.e-7);
		hMCAcc1M_Pt[i]->SetLineColor(2);
		hRCAcc1M_Pt[i]->SetLineColor(4);
		c1->cd(i%nPads+1);
		gPad->SetLogy(1);
		hMCAcc0M_Pt[i]->Draw("chist");
		hMCAcc1M_Pt[i]->Draw("chistsame");
		hRCAcc1M_Pt[i]->Draw("chistsame");
		leg->Clear();
		leg->AddEntry(hMCAcc0M_Pt[i],"cocktail (|Y_{ee}|<1)","pl");
		leg->AddEntry(hMCAcc1M_Pt[i],"cocktail (|Y_{ee}|<1, p_{T}^{e}>0.2, |#eta_{e}|<1)","pl");
		leg->AddEntry(hRCAcc1M_Pt[i],"cocktail (|Y_{ee}|<1, p_{T}^{e}>0.2, |#eta_{e}|<1, w/ eff.)","pl");
		leg->Draw("same");
		drawLatex(0.2,0.2,Form("%3.1f<p_{T}^{ee}<%3.1f GeV/c",PtLow[i],PtHi[i]),22,0.05,1);
		if(i%nPads==nPads-1) pdfAction(c1,ps);
	}
	if(nPtBins%nPads != 0) pdfAction(c1,ps);
	for(Int_t i=0;i<nPtBins;i++){
		hMCAcc0M_Pt[i]->GetXaxis()->UnZoom();
	}

	ps->On();
	ps->Close();

	TFile *f;
	if(withRho) {
		f= new TFile(Form("%s/cen%d%d_cocktail_withRho.root", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]),"recreate");
	}
	else {
		f= new TFile(Form("%s/cen%d%d_cocktail_withoutRho.root", dir.Data(), CentralityLow[mCenIdx], CentralityHi[mCenIdx]),"recreate");
	}
	f->cd();
	hMCAcc0MvsPt->Write();
	hMCAcc1MvsPt->Write();
	hUpMCAcc1MvsPt->Write();
	hLowMCAcc1MvsPt->Write();
	hMCAcc1MvsPtpi0dalitz->Write();
	hMCAcc1MvsPtetadalitz->Write();
	hMCAcc1MvsPtomegadalitz->Write();
	hMCAcc1MvsPtetaprimdalitz->Write();
	hMCAcc1MvsPtphidalitz->Write();
	//if(withRho) hMCAcc1MvsPtrho2ee->Write();
	hMCAcc1MvsPtomega2ee->Write();
	hMCAcc1MvsPtphi2ee->Write();
	hMCAcc1MvsPtjpsi2ee->Write();
	//hMCAcc1MvsPtpsi2ee->Write();
	hMCAcc1MvsPtcc2ee->Write();
	//hMCAcc1MvsPtbb2ee->Write();
	//hMCAcc1MvsPtdy2ee->Write();
	hMCAcc1MvsPtomega->Write();
	hMCAcc1MvsPtphi->Write();
	hRCAcc1MvsPt3D->Write();
	hMCAcc0M->Write();
	hMCAcc1M->Write();
	hUpMCAcc1M->Write();
	hLowMCAcc1M->Write();
	hRCAcc1M->Write();
	hRCAcc1MvsPt3Dpi0dalitz->Write();
	hRCAcc1MvsPt3Detadalitz->Write();
	hRCAcc1MvsPt3Domegadalitz->Write();
	hRCAcc1MvsPt3Detaprimdalitz->Write();
	hRCAcc1MvsPt3Dphidalitz->Write();
	hRCAcc1MvsPt3Domega2ee->Write();
	hRCAcc1MvsPt3Dphi2ee->Write();
	hRCAcc1MvsPt3Djpsi2ee->Write();
	hRCAcc1MvsPt3Dcc2ee->Write();
	for(Int_t i=0;i<nPtBins;i++){
		hMCAcc0M_Pt[i]->Write();
		hMCAcc1M_Pt[i]->Write();
		hUpMCAcc1M_Pt[i]->Write();
		hLowMCAcc1M_Pt[i]->Write();
		hRCAcc1M_Pt[i]->Write();
	}
	f->Close();

	cout<<"The program has completed!!"<<endl;

}
