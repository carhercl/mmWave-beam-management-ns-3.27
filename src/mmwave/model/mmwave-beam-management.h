/*
 * mmwave-beam-management.h
 *
 *  Created on: 12 oct. 2017
 *      Author: carlos
 */

#ifndef MMWAVE_BEAM_MANAGEMENT_H_
#define MMWAVE_BEAM_MANAGEMENT_H_


#include "ns3/object.h"
#include <ns3/spectrum-value.h>
#include <string.h>
#include "ns3/uinteger.h"
#include <complex>
#include <ns3/nstime.h>
#include <ns3/simple-ref-count.h>
#include <ns3/ptr.h>
#include <ns3/net-device-container.h>
#include <map>
#include <ns3/spectrum-signal-parameters.h>
#include <ns3/mobility-model.h>
#include <ns3/spectrum-propagation-loss-model.h>
#include <ns3/mmwave-phy-mac-common.h>
#include <ns3/random-variable-stream.h>


namespace ns3
{

typedef std::vector< std::complex<double> > complexVector_t;
typedef std::vector<complexVector_t> complex2DVector_t;
typedef std::pair<uint16_t, uint16_t > sinrKey;



struct BeamSweepingParams
{
	uint16_t m_currentBeamId;		// The current (or last) beam id used for beam sweeping
	Time m_steerBeamInterval;		// The time period the beam is changed to the next one
	complex2DVector_t m_codebook;	// The codebook used for beam sweeping
};

struct BeamPairInfoStruct
{
	Ptr<NetDevice> m_targetNetDevice;
	uint16_t m_txBeamId;
	uint16_t m_rxBeamId;
	SpectrumValue m_sinrPsd;
	double m_avgSinr;

	BeamPairInfoStruct();
};

class MmWaveBeamManagement : public Object
{

public:

	MmWaveBeamManagement ();
	virtual ~MmWaveBeamManagement ();

	static TypeId GetTypeId (void);

	void InitializeBeamSweepingTx(Time beamChangeTime);
	void InitializeBeamSweepingRx(Time beamChangeTime);

	void SetBeamChangeInterval (Time beamChangePeriod);
	void SetBeamSweepCodebook (complex2DVector_t codebook);

	complexVector_t GetBeamSweepVector ();
	complexVector_t GetBeamSweepVector (uint16_t index);

	void BeamSweepStepTx ();
	void BeamSweepStepRx ();
	void BeamSweepStep ();

	void DisplayCurrentBeamId ();

	uint16_t GetCurrentBeamId ();

	uint16_t GetNumBlocksSinceLastBeamSweepUpdate ();

	uint16_t IncreaseNumBlocksSinceLastBeamSweepUpdate ();

	void ResetNumBlocksSinceLastBeamSweepUpdate ();

	void AddEnbSinr (Ptr<NetDevice> enbNetDevice, uint16_t enbBeamId, uint16_t ueBeamId, SpectrumValue sinr);

	BeamPairInfoStruct FindBestScannedBeamPair ();
	BeamPairInfoStruct GetBestScannedBeamPair ();

	void UpdateBestScannedEnb();

	void ScheduleSsSlotSetStart(MmWavePhyMacCommon::SsBurstPeriods period);

	Time GetNextSsBlockTransmissionTime (Ptr<MmWavePhyMacCommon> mmWaveCommon, uint16_t currentSsBlock);

private:

	std::complex<double> ParseComplex (std::string strCmplx);

	/**
	* \brief Loads a generic beamforming matrix
	* \param Path to the input file to read
	* \return Array of complex values
	*/
	complex2DVector_t LoadCodebookFile (std::string inputFilename);

	void SetBestScannedEnb(BeamPairInfoStruct bestEnbBeamInfo);

	BeamSweepingParams m_beamSweepParams;

	Time m_lastBeamSweepUpdate;			// Timestamp of the last beam id modification

	uint16_t m_ssBlocksLastBeamSweepUpdate;	// Counter of ss blocks since the last beam update (when counter was reset)

	BeamPairInfoStruct m_bestScannedEnb;

	std::map <Ptr<NetDevice>,std::map <sinrKey,SpectrumValue>> m_enbSinrMap;	//Map to all the eNBs
//	std::map <Ptr<NetDevice>,std::map <sinrKey,float>> m_ueSinrMap;	//Map to all the UEs


};

}

#endif /* MMWAVE_BEAM_MANAGEMENT_H_ */
