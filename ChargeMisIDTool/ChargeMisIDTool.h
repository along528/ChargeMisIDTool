#ifndef WWW_CHARGEMISIDTOOL_H
#define WWW_CHARGEMISIDTOOL_H
#include "Definitions.h"
#include "D3PDReader/TruthParticleD3PDObject.h"
#include <string>
#include <vector>
#include "TH2F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TRandom.h"
#include <iostream>
#include <TLorentzVector.h>

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>;
#endif


using namespace std;
namespace WWW{


  class ChargeMisIDTool
  {

  public:
    ChargeMisIDTool(string chargeMisIDFilename,bool useMCCentralValue = false, string energyCorrectionFilename = "");//,const D3PDReader::MCEventD3PDObjectElement &mc);  
    ~ChargeMisIDTool();
    bool isInitialized() const { return m_bIsInitialized;};
    void setVerbosity(bool verbose){ 
    	cout << "WWW::ChargeMisIDTool setting verbosity to "<< (verbose ? "True" : "False") << endl;
	m_bVerbosity = verbose;
    }
    void printVerbose(string out){
    	if(!m_bVerbosity) return ;
	cout << "WWW::ChargeMisIDTool "<< out<< endl;
    }

    double getWeight(int mc_channel_number, int EventNumber,D3PDReader::TruthParticleD3PDObject &mc, vector<TLorentzVector> leptonsTLV, vector<int> leptonsCharge, vector<bool> leptonsIsElectron  );
    void setMatchDeltaR(double deltar){m_dMatchDeltaR = deltar;};
    bool doChargeFlip () const {return m_bDoChargeFlip;};
    double  getChargeFlippedLeptonEnergy (int,TLorentzVector) const ;
    int getChargeFlippedLeptonIndex () const {
    	if(!doChargeFlip()) cout << "WWW::ChargeMisIDTool::getChargeFlippedLeptonIndex WARNING This event does not have a charge flip"<<endl;
    	return m_iLeptonChargeFlipIndex;
    };
    bool checkStatus() ;
    void setSystematics(SysFlag sys){
    	m_eSysFlag = sys;
	if(m_bVerbosity){
		cout << "WWW::ChargeMisTool ";
		if (m_eSysFlag==Central) cout << "using Central Value";
		else if (m_eSysFlag==SysUp) cout << "using Systematic Up Variation;";
		else if (m_eSysFlag==SysDown) cout << "using Systematic Down Variation";
		else{
			cout << "ERROR uncexpected variation.  Using Central value.";
			m_eSysFlag = Central;
		}
		cout << endl;
	}
    }
    

  protected:
    double getChargeMisIDRate(TLorentzVector momentum);
    void setDoChargeFlip(bool dochargeflip, unsigned int index=-1){ 
    	m_bDoChargeFlip = dochargeflip;
	if(m_bDoChargeFlip) m_iLeptonChargeFlipIndex = index;
	else m_iLeptonChargeFlipIndex = -1;
    };
    unsigned int pickRandomIndex(int RandomSeed,vector<TLorentzVector>, vector<unsigned int>);
    vector<unsigned int> findLeptonsToChargeFlip(vector<TLorentzVector> truthLeptonsTLV, vector<int> truthLeptonsPDG,  vector<TLorentzVector> recoLeptonsTLV, vector<int> recoLeptonsCharge, vector<bool> recoLeptonsIsElectron  );
    SampleType classifySample(int );
    EventType classifyEvent(D3PDReader::TruthParticleD3PDObject &);
    EventType mapEventType(vector<int> );
    unsigned int getMatchingIndex(int tomatch, vector<int> list){
	unsigned int index = 999;
    	for (unsigned int i = 0;i< list.size();i++){
		if (list[i]==tomatch) {
			if(index!=999) {
				cout << "WWW::ChargeMisIDTool::getMatchingIndex ERROR found multiple matches" << endl;
				return 999;
			}
			index = i;
		}
	}

	return index;

    }
    





  private:
    bool m_bIsInitialized;
    SampleType m_eSample;
    EventType m_eEvent;
    bool m_bVerbosity;
    double m_dMatchDeltaR;
    TRandom *m_pRandom;
    bool m_bDoChargeFlip;
    int m_iLeptonChargeFlipIndex;
    TH2F *m_hRates;
    TH1F *m_hEnergyCorrection;
    TH1F *m_hEnergySmearing;
    TFile *m_pRatesFile;
    TFile *m_pEnergyCorrectionFile;
    ErrorFlag m_eStatus;
    SysFlag m_eSysFlag;
    bool m_bDoEnergyCorrections;
    bool m_bUseMC;

    vector<int> m_vTruthLeptonPDGID;
    vector<TLorentzVector> m_vTruthLeptonTLV;
    vector<int> m_vTruthLeptonIndex;





    

  };

};
#endif // WWW_CHARGEMISIDTOOL_H
