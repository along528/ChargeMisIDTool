#include <TFile.h>
#include <TMath.h>
#include <string>

#include <ChargeMisIDTool/ChargeMisIDTool.h>
#include <iostream>


using namespace std;


namespace WWW{
  ChargeMisIDTool::ChargeMisIDTool(string chargeMisIDFilename,bool useMCCentralValue,string energyCorrectionFilename) : m_bIsInitialized(false),m_eSample(OtherSample),m_eEvent(OtherEvent),m_bVerbosity(false),m_dMatchDeltaR(0.1),m_pRandom(0),m_bDoChargeFlip(false),m_iLeptonChargeFlipIndex(-1),m_hRates(0),m_hEnergyCorrection(0),m_hEnergySmearing(0),m_pRatesFile(0),m_pEnergyCorrectionFile(0),m_eStatus(OK),m_eSysFlag(Central),m_bDoEnergyCorrections(false),m_bUseMC(useMCCentralValue)  {
  	m_pRatesFile= new TFile(chargeMisIDFilename.c_str());
	if(!m_pRatesFile->IsOpen()) {
		cout << "WWW::ChargeMisIDTool ERROR Couldn't open rates file: " << chargeMisIDFilename << endl;
		return;
	}

	if(energyCorrectionFilename!=""){
		m_pEnergyCorrectionFile = new TFile(energyCorrectionFilename.c_str());
		if(!m_pEnergyCorrectionFile->IsOpen()) {
			cout << "WWW::ChargeMisIDTool ERROR Couldn't open energy correction file: " << energyCorrectionFilename << endl;
			return;
		}
	}


	//string histname = "MisIDRates";
	string histname = "MisIDRates_TotError";
	if(m_bUseMC) histname = "MisIDRates_MC_TotError";
	std::cout << histname << std::endl;
	m_hRates = (TH2F*)((TH2F*)m_pRatesFile->Get(histname.c_str()))->Clone("rates");


	if (!m_hRates){
		cout << "WWW::ChargeMisIDTool ERROR Couldn't configure using: " << chargeMisIDFilename << endl;
		return;
	}
	cout << "WWW::ChargeMisIDTool configured to use rates from: "<< chargeMisIDFilename << endl;

	if(m_pEnergyCorrectionFile && m_pEnergyCorrectionFile->GetListOfKeys()->Contains("hEnergyCorrection")) m_hEnergyCorrection = (TH1F*)((TH1F*)m_pEnergyCorrectionFile->Get("hEnergyCorrection"))->Clone("EnergyCorr");
	if(m_pEnergyCorrectionFile && m_pEnergyCorrectionFile->GetListOfKeys()->Contains("hEnergySmearing")) m_hEnergySmearing = (TH1F*)((TH1F*)m_pEnergyCorrectionFile->Get("hEnergySmearing"))->Clone("EnergySmear");

	if(!m_hEnergyCorrection) cout << "WWW::ChargeMisIDTool Couldn't find hEnergyCorrection" << endl;
	if(!m_hEnergySmearing)   cout << "WWW::ChargeMisIDTool Couldn't find hEnergySmearing" << endl;

	if(!m_hEnergyCorrection || !m_hEnergySmearing) {
		cout << "WWW::ChargeMisIDTool will NOT apply energy corrections" << endl;
		m_bDoEnergyCorrections = false;
	}
	else {
		cout << "WWW::ChargeMisIDTool energy corrections configured" << endl;
		m_bDoEnergyCorrections = true;
	}


	m_vTruthLeptonPDGID.clear();
	m_vTruthLeptonIndex.clear();
    	m_vTruthLeptonTLV.clear();

	m_pRandom = new TRandom();

	m_bIsInitialized = true;


  }


  ChargeMisIDTool::~ChargeMisIDTool() {
  	delete m_pRandom;
  	m_pRandom = NULL;
	m_pRatesFile->Close();
	delete m_pRatesFile;
	m_pRatesFile = NULL;
	if (m_pEnergyCorrectionFile){
		m_pEnergyCorrectionFile->Close();
		delete m_pEnergyCorrectionFile;
		m_pEnergyCorrectionFile = NULL;
	}
	delete m_pRandom;
	m_pRandom = NULL;
  }


  double ChargeMisIDTool::getWeight(int mc_channel_number, int EventNumber,D3PDReader::TruthParticleD3PDObject &mc, vector<TLorentzVector> leptonsTLV, vector<int> leptonsCharge, vector<bool> leptonsIsElectron  ){
    //units in MeV
    //initialize assuming nothing will happen
    m_eStatus = OK;
    setDoChargeFlip(false);
    double weight = 1.;
    //check to see that this is a sample which should be processed
    //if not then do nothing
    SampleType sampleType = classifySample(mc_channel_number);
    if(sampleType!=PowhegWlZll and sampleType!=PowhegWtauZll and sampleType!=PowhegZllZll) return 1.;


    //check to see if this event is an event which can
    //pass into the 0 SFOS region. If not then this 
    //event gets a weight of 0
    EventType eventType = classifyEvent(mc);
    if(eventType==OtherEvent){
    	printVerbose("Event is of type other");
	return 0.;
    }
    
    //If the event could migrate to the 0 SFOS region with a charge flip
    //we force the charge flip of the likely lepton
    
    //we search for those leptons which should be charge flipped
    vector<unsigned int> candidateLeptonsToChargeFlip = findLeptonsToChargeFlip(m_vTruthLeptonTLV,m_vTruthLeptonPDGID,leptonsTLV,leptonsCharge,leptonsIsElectron);

    if(!checkStatus()) return 0.;

    if(candidateLeptonsToChargeFlip.size()==0){
    	printVerbose("WARNING return 0 candidate leptons to charge flip");
	return 0.;
    }

    

	
    //the weight for the event is determined by summing the non-zero rates
    //to go to the 0 SFOS region for all of the leptons.  
    double ratesum = 0.;
    for(unsigned int i = 0;i<candidateLeptonsToChargeFlip.size();i++){
    	unsigned int lepIndex = candidateLeptonsToChargeFlip[i];
    	if(lepIndex >= leptonsTLV.size()) {
		cout << "WWW:ChargeMisIDTool::getWeight incorrect leptons indices" << endl;
		return 0;
	}
	ratesum += getChargeMisIDRate(leptonsTLV[lepIndex]);
    }


    //if more than 1 lepton could flip and bring to the 0 SFOS region, 
    //we pick one of the leptons randomly based on the rates of each lepton
    unsigned int chargeFlippedLeptonIndex = OK;
    if(candidateLeptonsToChargeFlip.size()==1) chargeFlippedLeptonIndex = candidateLeptonsToChargeFlip[0];
    else chargeFlippedLeptonIndex = pickRandomIndex(EventNumber,leptonsTLV,candidateLeptonsToChargeFlip);

    if(chargeFlippedLeptonIndex == Misc){
		cout << "WWW:ChargeMisIDTool::getWeight failure to pick lepton index" << endl;
		return 0;
    }

    

    weight = ratesum;
    if(weight < 0.) {
    	cout << "WWW::ChargeMisIDTool::getWeight WARNING returned bad weight" << endl;
	return 0.;
    }
    //Only after full processing do we determine if a charge flip will occur
    setDoChargeFlip(true,chargeFlippedLeptonIndex); 
    return weight;

  }

  double ChargeMisIDTool::getChargeMisIDRate(TLorentzVector momentum){
  	//expects MeV
	double Pt = momentum.Pt()/1000.;  //convert to GeV
	double AbsEta =  TMath::Abs(momentum.Eta());

	//cout << "Pt = " << Pt << " GeV, |Eta| = " << AbsEta << endl;

	//overflow values if Pt > 140 GeV, or Pt < 15 GeV
	//assumes last bin runs to at least 141 GeV
	if (Pt > 140) Pt = 141.;
	if (Pt < 15 ) Pt = 16.;

	//int binnum = m_hRates->FindBin(Pt,AbsEta); //flipped for ruiqi's rates
	int binnum = m_hRates->FindBin(AbsEta,Pt); 
	//cout << "bin# = " << binnum << endl;
	//
	
	double rate = 0.;
	if(m_eSysFlag == Central) rate = m_hRates->GetBinContent(binnum);
	else if(m_eSysFlag == SysUp) rate = m_hRates->GetBinContent(binnum) + m_hRates->GetBinError(binnum);
	else if(m_eSysFlag == SysDown) {
		rate = m_hRates->GetBinContent(binnum) - m_hRates->GetBinError(binnum);
		if (rate<=0.) rate = 0.;
	}
	else cout << "WWW::ChargeMisIDTool::getChargeMisIDRate unexpected variation" << endl;



	//cout << "rate = " << rate << endl;
	return rate;
  }

  vector<unsigned int> ChargeMisIDTool::findLeptonsToChargeFlip(vector<TLorentzVector> truthLeptonsTLV, vector<int> truthLeptonsPDG,  vector<TLorentzVector> recoLeptonsTLV, vector<int> recoLeptonsCharge, vector<bool> recoLeptonsIsElectron  ){
        if (truthLeptonsTLV.size()!=truthLeptonsPDG.size() ){
		cout << "WWW::ChargeMisIDTool::getLeptonToChargeFlip truth lepton indices don't match" << endl;
		m_eStatus = Misc;
		return vector<unsigned int> ();
	}
        if (recoLeptonsTLV.size()!=recoLeptonsCharge.size() || recoLeptonsTLV.size()!=recoLeptonsIsElectron.size() ){
		cout << "WWW::ChargeMisIDTool::getLeptonToChargeFlip reco lepton indices don't match" << endl;
		m_eStatus = Misc;
		return vector<unsigned int> ();
	}

        
	//search for reco leptons that match to the truth leptons we found before
	//when classifying the event
	vector<int> matchIndices ;
  	for (unsigned int recoIndex = 0 ;recoIndex < recoLeptonsTLV.size();recoIndex++){
	  int matchedTruthIndex = -1;
  	  for (unsigned int truthIndex = 0 ;truthIndex < truthLeptonsTLV.size();truthIndex++){
	  	double deltaR = recoLeptonsTLV[recoIndex].DeltaR(truthLeptonsTLV[truthIndex]);
		if (TMath::Abs(deltaR) < m_dMatchDeltaR)  {
			matchedTruthIndex = truthIndex;
			break;
		}
	  }
	  //check for unmatched reco leptons
	  if(matchedTruthIndex < 0) {
	  	m_eStatus = UnmatchedLepton;
		return vector<unsigned int> ();
	  }
	  //check for mismatching flavors
	  bool truthIsElectron = (TMath::Abs(truthLeptonsPDG[matchedTruthIndex]) == 11);
	  if(recoLeptonsIsElectron[recoIndex]!=truthIsElectron){
	  	m_eStatus = LeptonFlavorMismatch;
		return vector<unsigned int> ();
	  }

	  //check for overlaps
	  for(unsigned int i=0;i<matchIndices.size();i++)
	  	if (matchIndices[i] == matchedTruthIndex) {
			m_eStatus = MatchOverlap;
			return vector<unsigned int> ();
		}

	  //check to see if event already has charge flip
	  int truthCharge = (truthLeptonsPDG[matchedTruthIndex] > 0 ? -1 : 1);
	  if(recoLeptonsCharge[recoIndex]!=truthCharge) {
	  	m_eStatus = ChargeMisID;
		return vector<unsigned int> ();
	  }

	  matchIndices.push_back(matchedTruthIndex);
	}

	//check that matched truth leptons still agree with 
	//event type check
	vector<int> pdgs;
	for(unsigned int i =0;i<matchIndices.size();i++){
		int index = matchIndices[i];
		pdgs.push_back(truthLeptonsPDG[index]);
	}
	if( mapEventType(pdgs) != m_eEvent) {
		printVerbose("Re-classifying event");
		m_eEvent = mapEventType(pdgs);
		//cout << "WWW::ChargeMisIDTool::getLeptonToChargeFlip couldn't match to same event type" << endl;
		//m_eStatus = Misc;
		//return vector<unsigned int> ();
	}
	if(m_eEvent == OtherEvent) {
		printVerbose("WWW::ChargeMisIDTool::getLeptonToChargeFlip event type not recognized");
		return vector<unsigned int> ();
	}

    	vector<unsigned int> leptons;
	leptons.clear();
	if(m_eSample == PowhegWlZll || m_eSample == PowhegWtauZll || m_eSample == PowhegZllZll ){
		if(pdgs.size()!=3) {
			m_eStatus = Misc;
			cout << "WWW::ChargeMisIDTool::getLeptonToChargeFlip couldn't match 3 leptons" << endl;
			return vector<unsigned int> ();
		}
		if(pdgs.size()!=recoLeptonsTLV.size()) {
			m_eStatus = Misc;
			cout << "WWW::ChargeMisIDTool::getLeptonToChargeFlip couldn't match leptons" << endl;
			return vector<unsigned int> ();
		}

		if(m_eEvent == EmEmEp) {
			leptons.push_back( getMatchingIndex(-11,pdgs)); //return e+ index
			return leptons;
		}
		if(m_eEvent == EpEmEp) {
			leptons.push_back( getMatchingIndex(11,pdgs)); //return e- index
			return leptons;
		}
		//if(m_eEvent == MmEmEp) return getMatchingIndex(11,pdgs); //return e- index (bias to charge sum = +/- 1)
		//if(m_eEvent == MpEmEp) return getMatchingIndex(-11,pdgs); //return e+ index (bias to charge sum = +/- 1)
		//if you want to fill in the --- and +++ bins
		if(m_eEvent == MmEmEp || m_eEvent == MpEmEp) {
			leptons.push_back(getMatchingIndex(11,pdgs));
			leptons.push_back(getMatchingIndex(-11,pdgs));
			return leptons;
		}
	}


        

	m_eStatus = Misc;
	return vector<unsigned int> ();



  }
  unsigned int ChargeMisIDTool::pickRandomIndex(int RandomSeed,vector<TLorentzVector> recoLeptons, vector<unsigned int> leptonIndicesToConsider){
  	double rateSum = 0.;
	//randomly pick lepton based on how likely it is to flip 
	//For example: if one electron has a charge flip rate of 1 %
	//and the other has a charge flip rate of 0.5%
	//then given that there is a charge flip, 
	//the first will be flipped 2/3 of the time
	//and the other 1/3 of the time
	
	//first determine the rates for each lepton
  	vector<double> rates;
	for(unsigned int i =0;i<leptonIndicesToConsider.size();i++){
		unsigned int index = leptonIndicesToConsider[i];
		double rate = getChargeMisIDRate(recoLeptons[index]);
		//cout << "rate " << i<< ": " << rate << endl;
		rateSum+=rate;
		rates.push_back(rate);
	}

        //avoid divide by zero.
	//if total rate is zero then should get zero weight.
	//could occur with downward systematic
	if(rateSum<=0.) {
		cout << "WWW:ChargeMisIDTool::pickRandomIndex ratesum <=0 " << endl;
		m_eStatus = Misc;
		return Misc;
	}

	//normalize the rates with respect to each other
	//to determine the probability of one of them 
	//flipping given there is a charge flip
	double normalizedRateSum = 0.;
	vector<double> thresholds;
	for(unsigned int i = 0;i<rates.size();i++){
		normalizedRateSum += rates[i]/rateSum;
		//cout << "normalized rate sum = " << normalizedRateSum << endl;
		thresholds.push_back(normalizedRateSum);
	}
	//normalizedRateSum should sum to 1
	
	//pick a random number from 0 to 1
	m_pRandom->SetSeed(RandomSeed);
	double random = m_pRandom->Uniform();  
	//cout << "random num = " << random << endl;

	//select the number lepton based on the random number
	for (unsigned int i = 0; i< thresholds.size();i++)
		if(random < thresholds[i]) return leptonIndicesToConsider[i];
	

	cout << "WWW::ChargeMisIDTool::pickRandomIndex ERROR didn't pick an index" << endl;
	m_eStatus = Misc;
	return Misc;
  }
  EventType ChargeMisIDTool::mapEventType(vector<int> pdgs){
  	int sum = 0;
	//cout << "size: " << pdgs.size()<<endl;
  	for (unsigned  int i =0;i< pdgs.size();i++){
		if(pdgs[i]==11)             sum+=1;
		else if(pdgs[i]==-11)      sum+=10;
		else if(pdgs[i]==13)      sum+=100;
		else if(pdgs[i]==-13)    sum+=1000;
		else if(pdgs[i]==15)    sum+=10000;
		else if(pdgs[i]==-15)  sum+=100000;
		else{
			cout << "WWW::ChargeMisIDTool ERROR!  Unexpected saved truth particle with pdg = " << pdgs[i]<< endl;
			sum = 0;
			break;
		}
	}

	if(sum==EmEmEp){
		printVerbose("Event is e-e-e+");
		return EmEmEp;
	}
	else if(sum==EpEmEp) {
		printVerbose("Event is e+e-e+");
		return EpEmEp;
	}
	else if(sum==MmEmEp) {
		printVerbose("Event is mu-e-e+");
		return MmEmEp;
	}
	else if(sum==MpEmEp) {
		printVerbose("Event is mu+e-e+");
		return MpEmEp;
	}
	else if(sum==TmEmEp) {
		printVerbose("Event is tau-e-e+");
		return TmEmEp;
	}
	else if(sum==TpEmEp) {
		printVerbose("Event is tau+e-e+");
		return TpEmEp;
	}
	else if(sum==EmEpEmEp) {
		printVerbose("Event is e-e+e-e+");
		return EmEpEmEp;
	}
	else if(sum==MmMpEmEp) {
		printVerbose("Event is mu-mu+e-e+");
		return MmMpEmEp;
	}
	printVerbose("Event does not match expected");
	//cout << "sum = "<<sum << endl;
	return OtherEvent;


  }
  EventType ChargeMisIDTool::classifyEvent(D3PDReader::TruthParticleD3PDObject &mc){
	m_vTruthLeptonPDGID.clear();
	m_vTruthLeptonIndex.clear();
    	m_vTruthLeptonTLV.clear();
  	if(m_eSample == PowhegWlZll || m_eSample == PowhegWtauZll || m_eSample == PowhegZllZll ){
		for(int i = 0 ; i < mc.n(); i++){
			int absPDG = TMath::Abs(mc[i].pdgId());
			//if(! (absPDG == 11 || absPDG == 13 || absPDG == 15) ) continue;
			if(! (absPDG == 11 || absPDG == 13 ) ) continue;
			/*
			cout << i<< " "<<mc[i].pdgId() << " " << mc[i].status() << " " << endl;
			if(mc[i].parent_index().size() > 0) {
				int p = mc[i].parent_index()[0];
				cout << "\t"<<mc[p].pdgId() << " " << mc[p].status() << endl;
			}
			*/
			//e,mu from W,Z decay are status 1
			//taus from W decay are status 2
			//if(mc[i].status()!=1 &&  !(absPDG==15 && mc[i].status()==2)) continue;
			if(mc[i].status()!=1 ) continue;
			if(mc[i].parent_index().size() <= 0) continue;
			int parentIndex = mc[i].parent_index()[0];
			int parentAbsPDG = TMath::Abs(mc[parentIndex].pdgId());
			if(! (parentAbsPDG == 23 || parentAbsPDG == 24) && 
			   ! (m_eSample == PowhegWtauZll && parentAbsPDG == 15)) continue;
			m_vTruthLeptonPDGID.push_back(mc[i].pdgId());
			m_vTruthLeptonIndex.push_back(i);
			TLorentzVector tlv;
			tlv.SetPtEtaPhiM(mc[i].pt(),mc[i].eta(),mc[i].phi(),mc[i].m());
			m_vTruthLeptonTLV.push_back(tlv);
		}

		if((m_eSample == PowhegWlZll || m_eSample == PowhegWtauZll)  && m_vTruthLeptonPDGID.size()!=3) printVerbose("WARNGING: Event doesn't have 3 leptons");
		if((m_eSample == PowhegZllZll )  && m_vTruthLeptonPDGID.size()!=4) printVerbose("WARNGING: Event doesn't have 4 leptons");
  		m_eEvent = mapEventType(m_vTruthLeptonPDGID);
	}
	else m_eEvent = OtherEvent;
		
	return m_eEvent;

  	

  }
  double  ChargeMisIDTool::getChargeFlippedLeptonEnergy (int randomSeed,TLorentzVector momentum) const {
  	//please pass the uncorrected energy in the TLorentzVector!
	//according to the following: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/VectorBosonScattering2012#Charge_Flip_Background
	m_pRandom->SetSeed(randomSeed);

	if(!m_bDoEnergyCorrections) return momentum.E();

	double oldEnergy = momentum.E();

	//parameters are a function of |eta|
	//Method when using TGraph
	//double energyCorr = m_gEnergyCorrection->Eval(TMath::Abs(momentum.Eta()));
	//double energySmear = m_gEnergySmearing->Eval(TMath::Abs(momentum.Eta()));
	int energyBin = m_hEnergyCorrection->FindBin(TMath::Abs(momentum.Eta()));
	//put in units of MeV
	double energyCorr = m_hEnergyCorrection->GetBinContent(energyBin)*1000.; 
	double energySmear = m_hEnergySmearing->GetBinContent(energyBin)*1000.;
	double newEnergy = oldEnergy - 1.35 * energyCorr + m_pRandom->Gaus(0,1.25*energySmear);
	//std::cout << "energy corr: " << energyCorr << " smear: " << energySmear << " oldEnergy: " << oldEnergy << " newEnergy: "<< newEnergy << std::endl;

	return newEnergy;
	
  };

  SampleType ChargeMisIDTool::classifySample(int mc_channel_number){
    if(
      mc_channel_number == 185795 || //mc12_8TeV.185795.PowhegPythia8_AU2CT10_WmZ_3e_mll0p25_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.9655,1,0.051928,749996
      mc_channel_number == 185796 || //mc12_8TeV.185796.PowhegPythia8_AU2CT10_WmZ_e2mu_mll0p4614_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.6326,1,0.073874,749998
      mc_channel_number == 185798 || //mc12_8TeV.185798.PowhegPythia8_AU2CT10_WmZ_mu2e_mll0p25_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.9687,1,0.054302,749997
      mc_channel_number == 185799 || //mc12_8TeV.185799.PowhegPythia8_AU2CT10_WmZ_3mu_mll0p4614_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.6479,1,0.071268,749999
      mc_channel_number == 185804 || //mc12_8TeV.185804.PowhegPythia8_AU2CT10_WpZ_3e_mll0p25_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,1.416,1,0.053051,744897
      mc_channel_number == 185805 || //mc12_8TeV.185805.PowhegPythia8_AU2CT10_WpZ_e2mu_mll0p4614_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.9421,1,0.075904,749999
      mc_channel_number == 185807 || //mc12_8TeV.185807.PowhegPythia8_AU2CT10_WpZ_mu2e_mll0p25_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,1.412,1,0.055296,749998
      mc_channel_number == 185808 || //mc12_8TeV.185808.PowhegPythia8_AU2CT10_WpZ_3mu_mll0p4614_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.9572,1,0.073362,749499
      mc_channel_number == 129477 || //mc12_8TeV.129477.PowhegPythia8_AU2CT10_WZ_Wm11Z11_mll0p250d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129478 || //mc12_8TeV.129478.PowhegPythia8_AU2CT10_WZ_Wm11Z13_mll0p4614d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129480 || //mc12_8TeV.129480.PowhegPythia8_AU2CT10_WZ_Wm13Z11_mll0p250d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129481 || //mc12_8TeV.129481.PowhegPythia8_AU2CT10_WZ_Wm13Z13_mll0p4614d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129486 || //mc12_8TeV.129486.PowhegPythia8_AU2CT10_WZ_W11Z11_mll0p250d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129487 || //mc12_8TeV.129487.PowhegPythia8_AU2CT10_WZ_W11Z13_mll0p4614d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129489 || //mc12_8TeV.129489.PowhegPythia8_AU2CT10_WZ_W13Z11_mll0p250d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129490) //mc12_8TeV.129490.PowhegPythia8_AU2CT10_WZ_W13Z13_mll0p4614d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
    			m_eSample = PowhegWlZll;
    else if(
      mc_channel_number == 185797 || //mc12_8TeV.185797.PowhegPythia8_AU2CT10_WmZ_e2tau_mll3p804_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.1125,1,0.012544,1479993
      mc_channel_number == 185800 || //mc12_8TeV.185800.PowhegPythia8_AU2CT10_WmZ_mu2tau_mll3p804_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.1125,1,0.01258,1494893
      mc_channel_number == 185806 || //mc12_8TeV.185806.PowhegPythia8_AU2CT10_WpZ_e2tau_mll3p804_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.1755,1,0.013867,1499989
      mc_channel_number == 185809 || //mc12_8TeV.185809.PowhegPythia8_AU2CT10_WpZ_mu2tau_mll3p804_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.1755,1,0.013891,1499898
      mc_channel_number == 129479 || //mc12_8TeV.129479.PowhegPythia8_AU2CT10_WZ_Wm11Z15_mll3p804d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129482 || //mc12_8TeV.129482.PowhegPythia8_AU2CT10_WZ_Wm13Z15_mll3p804d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129488 || //mc12_8TeV.129488.PowhegPythia8_AU2CT10_WZ_W11Z15_mll3p804d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129491 ) //mc12_8TeV.129491.PowhegPythia8_AU2CT10_WZ_W13Z15_mll3p804d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
    			m_eSample = PowhegWlZtautau;
    else if(
      mc_channel_number == 185801 || //mc12_8TeV.185801.PowhegPythia8_AU2CT10_WmZ_tau2e_mll0p25_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.9687,1,0.012075,749998
      mc_channel_number == 185802 || //mc12_8TeV.185802.PowhegPythia8_AU2CT10_WmZ_tau2mu_mll0p4614_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.6326,1,0.01664,749998
      mc_channel_number == 185810 || //mc12_8TeV.185810.PowhegPythia8_AU2CT10_WpZ_tau2e_mll0p25_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,1.412,1,0.012105,749997
      mc_channel_number == 185811 || //mc12_8TeV.185811.PowhegPythia8_AU2CT10_WpZ_tau2mu_mll0p4614_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.9421,1,0.016718,749997
      mc_channel_number == 129483 || //mc12_8TeV.129483.PowhegPythia8_AU2CT10_WZ_Wm15Z11_mll0p250d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129484 || //mc12_8TeV.129484.PowhegPythia8_AU2CT10_WZ_Wm15Z13_mll0p4614d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129492 || //mc12_8TeV.129492.PowhegPythia8_AU2CT10_WZ_W15Z11_mll0p250d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129493 ) //mc12_8TeV.129493.PowhegPythia8_AU2CT10_WZ_W15Z13_mll0p4614d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
    			m_eSample = PowhegWtauZll;
    else if(
      mc_channel_number == 185803 || //mc12_8TeV.185803.PowhegPythia8_AU2CT10_WmZ_3tau_mll3p804_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.1108,1,0.0034037,1499994
      mc_channel_number == 185812 || //mc12_8TeV.185812.PowhegPythia8_AU2CT10_WpZ_3tau_mll3p804_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.172,1,0.0036427,1499893
      mc_channel_number == 129485 || //mc12_8TeV.129485.PowhegPythia8_AU2CT10_WZ_Wm15Z15_mll3p804d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 129494) //mc12_8TeV.129494.PowhegPythia8_AU2CT10_WZ_W15Z15_mll3p804d0_2LeptonFilter5.merge.NTUP_SMWZ.e1300_s1469_s1470_r3542_r3549_p1328/	
    			m_eSample = PowhegWtauZtautau;
    else if(
      mc_channel_number == 147197) //mc12_8TeV.147197.Sherpa_CT10_lllnu_WZ_l10.merge.NTUP_SMWZ.e1614_s1499_s1504_r3658_r3549_p1328/
			m_eSample = SherpaWZ;
    else if(
      mc_channel_number == 185813 || //mc12_8TeV.185813.PowhegPythia8_AU2CT10_ZZ_4e_mll4_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.07677,1,0.57204,3499890
      mc_channel_number == 185814 || //mc12_8TeV.185814.PowhegPythia8_AU2CT10_ZZ_2e2mu_mll4_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.1757,1,0.49893,3489992
      mc_channel_number == 185816 || //mc12_8TeV.185816.PowhegPythia8_AU2CT10_ZZ_4mu_mll4_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.07677,1,0.58293,3499995
      mc_channel_number == 126937 || //mc12_8TeV.126937.PowhegPythia8_AU2CT10_ZZ_4e_mll4_2pt5.merge.NTUP_SMWZ.e1280_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 126938 || //mc12_8TeV.126938.PowhegPythia8_AU2CT10_ZZ_2e2mu_mll4_2pt5.merge.NTUP_SMWZ.e1280_s1469_s1470_r3752_r3549_p1328/	
      mc_channel_number == 126940 ) //mc12_8TeV.126940.PowhegPythia8_AU2CT10_ZZ_4mu_mll4_2pt5.merge.NTUP_SMWZ.e1280_s1469_s1470_r3752_r3549_p1328/	
      			m_eSample = PowhegZllZll;
    else if(
      mc_channel_number == 185815 || //mc12_8TeV.185815.PowhegPythia8_AU2CT10_ZZ_2e2tau_mll4_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.1757,1,0.086032,3494790
      mc_channel_number == 185817 || //mc12_8TeV.185817.PowhegPythia8_AU2CT10_ZZ_2mu2tau_mll4_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.1757,1,0.087166,3488996
      mc_channel_number == 126939 || //mc12_8TeV.126939.PowhegPythia8_AU2CT10_ZZ_2e2tau_mll4_2pt5.merge.NTUP_SMWZ.e1280_s1469_s1470_r3542_r3549_p1328/	
      mc_channel_number == 126941 ) //mc12_8TeV.126941.PowhegPythia8_AU2CT10_ZZ_2mu2tau_mll4_2pt5.merge.NTUP_SMWZ.e1280_s1469_s1470_r3752_r3549_p1328/	
      			m_eSample = PowhegZllZtautau;
    else if(
      mc_channel_number == 185818 || //mc12_8TeV.185818.PowhegPythia8_AU2CT10_ZZ_4tau_mll4_TriLeptonFilter.merge.NTUP_SMWZ.e3069_s1773_s1776_r4485_r4540_p1328/,0.07677,1,0.0076557,3499992
      mc_channel_number == 126942 ) //mc12_8TeV.126942.PowhegPythia8_AU2CT10_ZZ_4tau_mll4_2pt5.merge.NTUP_SMWZ.e1280_s1469_s1470_r3542_r3549_p1328/	
      			m_eSample = PowhegZtautauZtautau;
    else m_eSample = OtherSample;

    //cout<< "Sample Classified as "<<m_eSample << endl;





    //


    return m_eSample;
    //not considering Sherpa ZZ

  }

  bool ChargeMisIDTool::checkStatus(){
  	if(m_eStatus == OK) return true;
        if(m_eStatus == UnmatchedLepton)
    	   printVerbose("WARNING couldn't match all leptons");
    	if(m_eStatus == LeptonFlavorMismatch)
    		printVerbose("WARNING matched leptons are of mis-matched flavor");
    	if(m_eStatus == MatchOverlap)
    		printVerbose("WARNING couldn't matched leptons overlap");
    		//ignornig  events with actual charge flip for now
		//is this what we want to do??
    	if(m_eStatus == ChargeMisID)
    		printVerbose("WARNING event already has charge mis-id");
    	if(m_eStatus == Misc)
    		printVerbose("WARNING couldn't get weight");
	
	return false;


  }

};
