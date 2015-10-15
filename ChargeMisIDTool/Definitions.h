#ifndef WWW_DEFINITIONS_H
#define WWW_DEFINITIONS_H

using namespace std;
namespace WWW{
  typedef enum SysFlag{
  	Central =0,
	SysUp,
	SysDown
  } SysFlag;

  typedef enum ErrorFlag{
  	OK = 500,
  	UnmatchedLepton = 901,
  	LeptonFlavorMismatch,
	MatchOverlap,
	ChargeMisID,
	Misc =999
  } ErrorFlag;
  typedef enum SampleType{
  	PowhegWlZll = 0,
  	PowhegWlZtautau,
	PowhegWtauZll,
	PowhegWtauZtautau,
	PowhegZllZll,
	PowhegZllZtautau,
	PowhegZtautauZtautau,
	SherpaWZ,
	SherpaZZ,
	OtherSample
  } SampleType;

  /* Based on Truth ID of leptons in event
   * currently ignoring ZZ, any form of tau decay
   * and Z->mumu decay.
   * These are the processes which should give the largest rates to 
   * the 0 SFOS region.
   * The mapping is from summing
   * nElectronMinus*1. + nElectronPlus*100. + nMuonMinus*10000. + nMuonPlus*1000000.
   */

  typedef enum EventType{
  	OtherEvent  = 0,
	EmEmEp = 12,
	EpEmEp = 21,
	MmEmEp = 111,
	MpEmEp = 1011,
	TmEmEp = 10011,
	TpEmEp = 100011,
	EmEpEmEp = 22,
	MmMpEmEp = 1111
  } EventType;

  typedef enum ParticleType{
  	OtherParticle = 0,
  	Electron = 11,
	Muon = 13,
	Tau = 15
  } ParticleType;


};
#endif // WWW_DEFINTIONS_H
