/*
 * mmwave-beam-management.cc
 *
 *  Created on: 12 oct. 2017
 *      Author: Carlos Herranz
 */


#include "mmwave-beam-management.h"

#include <ns3/log.h>
#include <fstream>
#include <ns3/simulator.h>
#include <ns3/abort.h>
#include <ns3/mmwave-enb-net-device.h>
#include <ns3/mmwave-ue-net-device.h>
#include <ns3/mmwave-ue-phy.h>
#include <ns3/antenna-array-model.h>
#include <ns3/node.h>
#include <algorithm>
#include <ns3/double.h>
#include <ns3/boolean.h>
#include "mmwave-spectrum-value-helper.h"

namespace ns3{

NS_LOG_COMPONENT_DEFINE ("MmWaveBeamManagement");

NS_OBJECT_ENSURE_REGISTERED (MmWaveBeamManagement);




BeamPairInfoStruct::BeamPairInfoStruct()
{
	m_avgSinr = -100;
	m_txBeamId = 65535;
	m_rxBeamId = 65535;
	m_sinrPsd = -1.0;
}


MmWaveBeamManagement::MmWaveBeamManagement()
{
	m_beamSweepParams.m_currentBeamId = 0;
	m_ssBlocksLastBeamSweepUpdate = 0;
	m_maxNumBeamPairCandidates = 5;
}


MmWaveBeamManagement::~MmWaveBeamManagement()
{
	std::map< Ptr<NetDevice>, std::map <sinrKey,SpectrumValue>>::iterator it1;
	for (it1 = m_enbSinrMap.begin(); it1 != m_enbSinrMap.end(); ++it1)
	{
		it1->second.clear();
	}
	m_enbSinrMap.clear();
}


TypeId
MmWaveBeamManagement::GetTypeId (void)
{
	static TypeId tid = TypeId ("ns3::MmWaveBeamManagement")
		.SetParent<Object> ()
	;
  	return tid;
}


void
MmWaveBeamManagement::InitializeBeamSweepingTx(Time beamChangeTime)
{
	m_beamSweepParams.m_currentBeamId = 0;

//	std::string txFilePath = "src/mmwave/model/BeamFormingMatrix/TxCodebook.txt";
	std::string txFilePath = "src/mmwave/model/BeamFormingMatrix/KronCodebook16h4v.txt";
//	std::string txFilePath = "../InputFiles/KronCodebook16h4v.txt";

	m_beamSweepParams.m_codebook = LoadCodebookFile(txFilePath);
	this->SetBeamChangeInterval(beamChangeTime);
	m_lastBeamSweepUpdate = Simulator::Now();
	NS_LOG_INFO ("InitializeBeamSweepingTx");
//	DisplayCurrentBeamId();

}

void
MmWaveBeamManagement::InitializeBeamSweepingRx(Time beamChangeTime)
{
	m_beamSweepParams.m_currentBeamId = 0;

//	std::string rxFilePath = "src/mmwave/model/BeamFormingMatrix/RxCodebook.txt";
	std::string rxFilePath = "src/mmwave/model/BeamFormingMatrix/KronCodebook8h2v.txt";
//	std::string rxFilePath = "../InputFiles/KronCodebook8h2v.txt";

	m_beamSweepParams.m_codebook = LoadCodebookFile(rxFilePath);
	this->SetBeamChangeInterval(beamChangeTime);
	NS_LOG_INFO ("InitializeBeamSweepingRx");
	m_lastBeamSweepUpdate = Simulator::Now();
//	DisplayCurrentBeamId();

}

void MmWaveBeamManagement::SetBeamChangeInterval (Time period)
{
	m_beamSweepParams.m_steerBeamInterval = period;
}

void MmWaveBeamManagement::SetBeamSweepCodebook (complex2DVector_t codebook)
{
	m_beamSweepParams.m_codebook = codebook;
}

complexVector_t MmWaveBeamManagement::GetBeamSweepVector ()
{
	return m_beamSweepParams.m_codebook.at(m_beamSweepParams.m_currentBeamId);

}

complexVector_t MmWaveBeamManagement::GetBeamSweepVector (uint16_t index)
{
	return m_beamSweepParams.m_codebook.at(index);

}

void MmWaveBeamManagement::BeamSweepStep()
{

	Time currentTime = Simulator::Now();
	uint16_t numCodes = m_beamSweepParams.m_codebook.size();
//	NS_LOG_INFO ("[" << currentTime << "] Beam id " << m_beamSweepParams.m_currentBeamId << " of " << numCodes - 1);
	m_beamSweepParams.m_currentBeamId = (m_beamSweepParams.m_currentBeamId + 1) % numCodes;
//	m_lastBeamSweepUpdate = currentTime;
//	NS_LOG_INFO ("[" << currentTime << "] Beam id " << m_beamSweepParams.m_currentBeamId << " of " << numCodes - 1);

}

void MmWaveBeamManagement::BeamSweepStepTx()
{
//	Time currentTime = Simulator::Now();
//	if (currentTime == 0 || currentTime >= m_lastBeamSweepUpdate + m_beamSweepParams.m_steerBeamInterval)
	{
		BeamSweepStep();
	}
//	else
//	{
//		NS_LOG_INFO ("MmWaveBeamManagement::BeamSweepStepTx() This should not happen");
//	}
}

void MmWaveBeamManagement::BeamSweepStepRx()
{
//	Time currentTime = Simulator::Now();
//	if (currentTime >= m_lastBeamSweepUpdate + m_beamSweepParams.m_steerBeamInterval)
	{
		BeamSweepStep();
	}
}


void MmWaveBeamManagement::DisplayCurrentBeamId ()
{
	Time currentTime = Simulator::Now();
	uint16_t numCodes = m_beamSweepParams.m_codebook.size();
	NS_LOG_INFO ("[" << currentTime << "] Beam id " << m_beamSweepParams.m_currentBeamId << " of " << numCodes - 1);
}

uint16_t
MmWaveBeamManagement::GetCurrentBeamId ()
{
	return m_beamSweepParams.m_currentBeamId;
}


uint16_t
MmWaveBeamManagement::GetNumBlocksSinceLastBeamSweepUpdate ()
{
	return m_ssBlocksLastBeamSweepUpdate;
}

uint16_t
MmWaveBeamManagement::IncreaseNumBlocksSinceLastBeamSweepUpdate ()
{
	m_ssBlocksLastBeamSweepUpdate++;
	return m_ssBlocksLastBeamSweepUpdate;

}

void MmWaveBeamManagement::ResetNumBlocksSinceLastBeamSweepUpdate ()
{
	m_ssBlocksLastBeamSweepUpdate = 0;
}

void
MmWaveBeamManagement::AddEnbSinr (Ptr<NetDevice> enbNetDevice, uint16_t enbBeamId, uint16_t ueBeamId, SpectrumValue sinr)
{
	sinrKey key = std::make_pair(enbBeamId,ueBeamId);
	std::map< Ptr<NetDevice>, std::map <sinrKey,SpectrumValue>>::iterator it1 = m_enbSinrMap.find(enbNetDevice);
	if (it1 != m_enbSinrMap.end ())
	{
		std::map <sinrKey,SpectrumValue>::iterator it2 = it1->second.find(key);
		if (it2 != it1->second.end())
		{
			it2->second = sinr;
		}
		else
		{
			it1->second.insert(make_pair(key,sinr));
		}
//		m_enbSinrMap.erase (it1);
	}
	else
	{
		std::map <sinrKey,SpectrumValue> m2;
		m2.insert(make_pair(key,sinr));
		m_enbSinrMap.insert(make_pair(enbNetDevice,m2));
	}
//	std::cout << Simulator::Now().GetNanoSeconds() << " " << enbBeamId << " " << ueBeamId << " "
//			<< Sum(sinr)/sinr.GetSpectrumModel()->GetNumBands() << std::endl;

}

void
MmWaveBeamManagement::FindBeamPairCandidatesSinr ()
{
	for (std::map< Ptr<NetDevice>, std::map <sinrKey,SpectrumValue>>::iterator it1 = m_enbSinrMap.begin();
				it1 != m_enbSinrMap.end();
				++it1)
	{
		Ptr<NetDevice> pDevice = it1->first;
		std::vector<BeamPairInfoStruct> candidateBeamPairs;

		// Internal SINR values to control the size of the candidate beam vector.
		double minSinr=0;
		uint8_t numBeamPairs = 0;
		//double maxSinr=0;

		for (std::map <sinrKey,SpectrumValue>::iterator it2 = it1->second.begin();
						it2 != it1->second.end();
						++it2)
		{
			int nbands = it2->second.GetSpectrumModel ()->GetNumBands ();
			double avgSinr = Sum (it2->second)/nbands;

			// Skip the current pair of beams if they provide lower SINR than the current candidate pairs
			if(avgSinr < minSinr)
				continue;

			BeamPairInfoStruct beamPair;
			beamPair.m_avgSinr = avgSinr;
			beamPair.m_sinrPsd = it2->second;
			beamPair.m_targetNetDevice = it1->first;
			beamPair.m_txBeamId = it2->first.first;
			beamPair.m_rxBeamId = it2->first.second;

			if(candidateBeamPairs.empty())
			{
				candidateBeamPairs.push_back(beamPair);
				numBeamPairs = 1;
			}
//			if(candidateBeams.size() < m_numMaxCandidateBeams)
//			{
//
//			}
			else
			{
				bool iterate = true;
				std::vector<BeamPairInfoStruct>::iterator itV = candidateBeamPairs.begin();
				while (iterate && itV != candidateBeamPairs.end())
				// The candidate beam pairs are stored in SINR descending order.
//				for (std::vector<BeamPairInfoStruct>::iterator itV = candidateBeamPairs.begin();
//						itV != candidateBeamPairs.end();
//						itV++)
				{
					if (itV->m_avgSinr < avgSinr)
					{
						candidateBeamPairs.insert(itV,beamPair);
						numBeamPairs++;
						iterate = false;
//						if(candidateBeams.size() > m_numMaxCandidateBeams)
//						{
//							candidateBeams.erase(candidateBeams.end());
//						}
					}
					++itV;
				}
				// Only get best m_numMaxCandidateBeams pair of beams
				if(candidateBeamPairs.size() > m_maxNumBeamPairCandidates)
				{
					candidateBeamPairs.resize(m_maxNumBeamPairCandidates);
					numBeamPairs = m_maxNumBeamPairCandidates;
				}
			}

			//Update internal sinr values
			//maxSinr = candidateBeamPairs.at(0).m_avgSinr;
			minSinr = candidateBeamPairs.at(numBeamPairs-1).m_avgSinr;
		}

		//to map:
		std::map <Ptr<NetDevice>,std::vector<BeamPairInfoStruct>>::iterator itMap = m_candidateBeamsMap.find(pDevice);
		if(itMap == m_candidateBeamsMap.end())
		{
			m_candidateBeamsMap.insert(std::pair<Ptr<NetDevice>,std::vector<BeamPairInfoStruct>>(pDevice,candidateBeamPairs));
		}
		else
		{
			itMap->second = candidateBeamPairs;
		}
	}
}


void
MmWaveBeamManagement::FindBeamPairCandidatesVicinity ()
{
	int nbands;
	for (std::map< Ptr<NetDevice>, std::map <sinrKey,SpectrumValue>>::iterator it1 = m_enbSinrMap.begin();
				it1 != m_enbSinrMap.end();
				++it1)
	{
		Ptr<NetDevice> pDevice = it1->first;
		std::vector<BeamPairInfoStruct> candidateBeamPairs;

		// Internal SINR value to control the size of the candidate beam vector.
		double minSinr=0;
		BeamPairInfoStruct beamPair, beamPairExtra1, beamPairExtra2, beamPairExtra3, beamPairExtra4, beamPairExtra5;

		for (std::map <sinrKey,SpectrumValue>::iterator it2 = it1->second.begin();
						it2 != it1->second.end();
						++it2)
		{
			nbands = it2->second.GetSpectrumModel ()->GetNumBands ();
			double avgSinr = Sum (it2->second)/nbands;

			// Skip the current pair of beams if they provide lower SINR than the current candidate pairs
			if(avgSinr > minSinr)
			{
				beamPair.m_avgSinr = avgSinr;
				beamPair.m_sinrPsd = it2->second;
				beamPair.m_targetNetDevice = pDevice;
				beamPair.m_txBeamId = it2->first.first;
				beamPair.m_rxBeamId = it2->first.second;
				minSinr = avgSinr;
			}

		}

		candidateBeamPairs.push_back(beamPair);
		uint16_t txBeamId = beamPair.m_txBeamId;
		uint16_t rxBeamId = beamPair.m_rxBeamId;
		// FIXME: Hard-coded for the 8x2 array at the UE side. Make it UE-independent
		uint16_t beamIdRight = (rxBeamId + 1)%16;
		uint16_t beamIdLeft = (rxBeamId - 1)%16;
		uint16_t beamIdVertical = (rxBeamId + 8)%16;
		uint16_t beamIdVerticalRight = (rxBeamId + 8 + 1)%16;
		uint16_t beamIdVerticalLeft = (rxBeamId + 8 - 1)%16;

		std::map <sinrKey,SpectrumValue>::iterator itExtra1 =
				it1->second.find(std::make_pair(txBeamId,beamIdRight));
		beamPairExtra1.m_avgSinr = Sum(itExtra1->second)/nbands;
		beamPairExtra1.m_sinrPsd = itExtra1->second;
		beamPairExtra1.m_targetNetDevice = pDevice;
		beamPairExtra1.m_txBeamId = itExtra1->first.first;
		beamPairExtra1.m_rxBeamId = itExtra1->first.second;
		candidateBeamPairs.push_back(beamPairExtra1);

		std::map <sinrKey,SpectrumValue>::iterator itExtra2 =
				it1->second.find(std::make_pair(txBeamId,beamIdLeft));
		beamPairExtra2.m_avgSinr = Sum(itExtra2->second)/nbands;
		beamPairExtra2.m_sinrPsd = itExtra2->second;
		beamPairExtra2.m_targetNetDevice = pDevice;
		beamPairExtra2.m_txBeamId = itExtra2->first.first;
		beamPairExtra2.m_rxBeamId = itExtra2->first.second;
		candidateBeamPairs.push_back(beamPairExtra2);

		std::map <sinrKey,SpectrumValue>::iterator itExtra3 =
				it1->second.find(std::make_pair(txBeamId,beamIdVertical));
		beamPairExtra3.m_avgSinr = Sum(itExtra3->second)/nbands;
		beamPairExtra3.m_sinrPsd = itExtra3->second;
		beamPairExtra3.m_targetNetDevice = pDevice;
		beamPairExtra3.m_txBeamId = itExtra3->first.first;
		beamPairExtra3.m_rxBeamId = itExtra3->first.second;
		candidateBeamPairs.push_back(beamPairExtra3);

		std::map <sinrKey,SpectrumValue>::iterator itExtra4 =
				it1->second.find(std::make_pair(txBeamId,beamIdVerticalRight));
		beamPairExtra4.m_avgSinr = Sum(itExtra4->second)/nbands;
		beamPairExtra4.m_sinrPsd = itExtra4->second;
		beamPairExtra4.m_targetNetDevice = pDevice;
		beamPairExtra4.m_txBeamId = itExtra4->first.first;
		beamPairExtra4.m_rxBeamId = itExtra4->first.second;
		candidateBeamPairs.push_back(beamPairExtra4);

		std::map <sinrKey,SpectrumValue>::iterator itExtra5 =
				it1->second.find(std::make_pair(txBeamId,beamIdVerticalLeft));
		beamPairExtra5.m_avgSinr = Sum(itExtra5->second)/nbands;
		beamPairExtra5.m_sinrPsd = itExtra5->second;
		beamPairExtra5.m_targetNetDevice = pDevice;
		beamPairExtra5.m_txBeamId = itExtra5->first.first;
		beamPairExtra5.m_rxBeamId = itExtra5->first.second;
		candidateBeamPairs.push_back(beamPairExtra5);


		//to map:
		std::map <Ptr<NetDevice>,std::vector<BeamPairInfoStruct>>::iterator itMap = m_candidateBeamsMap.find(pDevice);
		if(itMap == m_candidateBeamsMap.end())
		{
			m_candidateBeamsMap.insert(std::pair<Ptr<NetDevice>,std::vector<BeamPairInfoStruct>>(pDevice,candidateBeamPairs));
		}
		else
		{
			itMap->second = candidateBeamPairs;
		}
	}
}

uint16_t
MmWaveBeamManagement::GetMaxNumBeamPairCandidates()
{
	return m_maxNumBeamPairCandidates;
}

void MmWaveBeamManagement::SetMaxNumBeamPairCandidates(uint16_t nBeamPairs)
{
	m_maxNumBeamPairCandidates = nBeamPairs;
}


BeamPairInfoStruct
MmWaveBeamManagement::FindBestScannedBeamPair ()
{
	BeamPairInfoStruct bestPairInfo;
	double bestAvgSinr = -1.0;

	// First create all beam pair candidates
	//FindBeamPairCandidatesSinr();
	FindBeamPairCandidatesVicinity();

	// Now iterate along the map and find the best gNB providing the best beam pairs in terms of SINR
	for (std::map <Ptr<NetDevice>,std::vector<BeamPairInfoStruct>>::iterator it = m_candidateBeamsMap.begin();
			it != m_candidateBeamsMap.end();
			++it)
	{
		// The candidate beam pairs are ordered in SINR descending order.
		double currentBeamPairAvgSinr = it->second.at(0).m_avgSinr;
		if (currentBeamPairAvgSinr > bestAvgSinr)
		{
			bestPairInfo = it->second.at(0);
			bestAvgSinr = currentBeamPairAvgSinr;
		}
	}

	return bestPairInfo;
}


BeamPairInfoStruct
MmWaveBeamManagement::GetBestScannedBeamPair()
{
	return m_bestScannedEnb;
}


void MmWaveBeamManagement::UpdateBestScannedEnb()
{
	// TODO: Modify this function to support multiple beams monitoring.
	// Update the candidate beam set from all available gNBs.


	// Get the strongest pair of beams from the best gNB.
	BeamPairInfoStruct bestScannedBeamPair = FindBestScannedBeamPair();
	SetBestScannedEnb(bestScannedBeamPair);
	Time currentTime = Simulator::Now();
//	NS_LOG_INFO("[" << currentTime.GetSeconds() <<"]Best beam pair update: tx=" << bestScannedBeamPair.m_txBeamId <<
//			" rx=" << bestScannedBeamPair.m_rxBeamId <<
//			" avgSinr=" << bestScannedBeamPair.m_avgSinr);
	std::cout << "[" << currentTime.GetSeconds() <<"]Best beam pair update: tx=" << bestScannedBeamPair.m_txBeamId <<
				" rx=" << bestScannedBeamPair.m_rxBeamId <<
				" avgSinr=" << bestScannedBeamPair.m_avgSinr << std::endl;
}


void MmWaveBeamManagement::ScheduleSsSlotSetStart(MmWavePhyMacCommon::SsBurstPeriods period)
{
	m_ssBlocksLastBeamSweepUpdate = 0;
	Simulator::Schedule(MicroSeconds(1000*period)-NanoSeconds(1),&MmWaveBeamManagement::ScheduleSsSlotSetStart,this,period);
}


std::complex<double>
MmWaveBeamManagement::ParseComplex (std::string strCmplx)
{
    double re = 0.00;
    double im = 0.00;
    size_t findj = 0;
    std::complex<double> out_complex;

    findj = strCmplx.find("i");
    if( findj == std::string::npos )
    {
        im = -1.00;
    }
    else
    {
        strCmplx[findj] = '\0';
    }
    if( ( strCmplx.find("+",1) == std::string::npos && strCmplx.find("-",1) == std::string::npos ) && im != -1 )
    {
        /* No real value */
        re = -1.00;
    }
    std::stringstream stream( strCmplx );
    if( re != -1.00 )
    {
        stream>>re;
    }
    else
    {
        re = 0;
    }
    if( im != -1 )
    {
        stream>>im;
    }
    else
    {
        im = 0.00;
    }
    //  std::cout<<" ---  "<<re<<"  "<<im<<std::endl;
    out_complex = std::complex<double>(re,im);
    return out_complex;
}


complex2DVector_t
MmWaveBeamManagement::LoadCodebookFile (std::string inputFilename)
{
	//std::string filename = "src/mmwave/model/BeamFormingMatrix/TxAntenna.txt";
	complex2DVector_t output;
	NS_LOG_FUNCTION (this << "Loading TxAntenna file " << inputFilename);
	std::ifstream singlefile;
	std::complex<double> complexVar;
	singlefile.open (inputFilename.c_str (), std::ifstream::in);

	NS_LOG_INFO (this << " File: " << inputFilename);
	NS_ASSERT_MSG(singlefile.good (), inputFilename << " file not found");
    std::string line;
    std::string token;
    while( std::getline(singlefile, line) ) //Parse each line of the file
    {
    	complexVector_t txAntenna;
        std::istringstream stream(line);
        while( getline(stream,token,',') ) //Parse each comma separated string in a line
        {
        	complexVar = ParseComplex(token);
		    txAntenna.push_back(complexVar);
		}
        output.push_back(txAntenna);
	}
    return output;
//    NS_LOG_INFO ("TxAntenna[instance:"<<g_enbAntennaInstance.size()<<"][antennaSize:"<<g_enbAntennaInstance[0].size()<<"]");
}

void MmWaveBeamManagement::SetBestScannedEnb(BeamPairInfoStruct bestEnbBeamInfo)
{
	if (bestEnbBeamInfo.m_targetNetDevice)
	{
		m_bestScannedEnb = bestEnbBeamInfo;
	}
}

Time MmWaveBeamManagement::GetNextSsBlockTransmissionTime (Ptr<MmWavePhyMacCommon> mmWaveCommon, uint16_t currentSsBlock)
{
	Time nextScheduledSsBlock;
	//m_currentSsBlockSlotId
	uint16_t beamPatternLength = mmWaveCommon->GetSsBlockPatternLength();
	uint16_t nextSsBlockId = currentSsBlock + 1;
	uint16_t currentOfdmSymbol = mmWaveCommon->GetSsBurstOfdmIndex(currentSsBlock);
	uint16_t nextOfdmSymbol;
	Time now = Simulator::Now();
	if(nextSsBlockId < beamPatternLength)
	{
		nextOfdmSymbol = mmWaveCommon->GetSsBurstOfdmIndex(nextSsBlockId);
		uint16_t numSymbols = nextOfdmSymbol-currentOfdmSymbol;
		nextScheduledSsBlock = MicroSeconds(mmWaveCommon->GetSymbolPeriod()*numSymbols);
	}
	else
	{
		nextOfdmSymbol = mmWaveCommon->GetSsBurstOfdmIndex(0);
		m_lastBeamSweepUpdate += MilliSeconds(mmWaveCommon->GetSsBurstSetPeriod());
		Time t2 = NanoSeconds(1000.0 * mmWaveCommon->GetSymbolPeriod() * nextOfdmSymbol);
		nextScheduledSsBlock = m_lastBeamSweepUpdate - now + t2;
	}


	return nextScheduledSsBlock;
}

}

