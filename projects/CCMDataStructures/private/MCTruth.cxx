/*!**********************************************
 * \file MCTruth.cxx
 * \brief Source code for the #MCTruth class
 * \author R. T. Thornton (LANL)
 * \date November 11, 2020
 ***********************************************/

#include <memory>
#include <numeric>
#include <algorithm>

#include "CCMAnalysis/CCMDataStructures/MCTruth.h"

#include "CCMAnalysis/CCMUtils/MsgLog.h"

ClassImp(MCTruth)

/*!**********************************************
 * \fn MCTruth::MCTruth()
 * \brief Default Constructor
 ***********************************************/
MCTruth::MCTruth() : TObject()
{
	Reset(0);
}

/*!**********************************************
 * \fn MCTruth::MCTruth(int particleID)
 * \brief Constructor (calls #Reset)
 * \param[in] particleID The PDG number for the particle that was simulated
 ***********************************************/
MCTruth::MCTruth(int particleID)
: TObject()
{
  Reset(particleID);
}

/*!**********************************************
 * \fn void MCTruth::Reset(int particleID)
 * \brief Resets the values to what is passed
 * \param[in] particleID The PDG number for the particle that was simulated
 ***********************************************/
void MCTruth::Reset(int particleID)
{
	fParticleID = particleID;
	fMomX = 0;
	fPosX = 0;
	fMomY = 0;
	fPosY = 0;
	fMomZ = 0;
	fPosZ = 0;

	fHitEnergy.clear();
	fHitAngle.clear();
	fHitUncoated.clear();
	fHitRow.clear();
	fHitCol.clear();
	fHitTime.clear();
  fPassedQF.clear();
}

/*!**********************************************
 * \fn MCTruth::MCTruth(const MCTruth & rhs)
 * \brief Copy constructor
 * \param[in] rhs The object to copy
 ***********************************************/
MCTruth::MCTruth(const MCTruth & rhs) : TObject(rhs)
{
  this->operator=(rhs);
}

/*!**********************************************
 * \fn MCTruth::~MCTruth()
 * \brief Destructor
 ***********************************************/
MCTruth::~MCTruth()
{
  // empty
	fHitEnergy.clear();
	fHitAngle.clear();
	fHitUncoated.clear();
	fHitRow.clear();
	fHitCol.clear();
	fHitTime.clear();
  fPassedQF.clear();
}

/*!**********************************************
 * \fn void MCTruth::SetParticleMomentum(float x, float y, float z)
 * \brief Sets the particles initial momentum 
 * \param[in] x The x coordinate 
 * \param[in] y The Y coordinate
 * \param[in] z The Z coordinate
 ***********************************************/
void MCTruth::SetParticleMomentum(float x, float y, float z)
{
	fMomX = x;
	fMomY = y;
	fMomZ = z;
	return;
}

/*!**********************************************
 * \fn void MCTruth::SetParticlePosition(float x, float y, float z)
 * \brief Sets the particles initial position
 * \param[in] x The x coordinate 
 * \param[in] y The Y coordinate
 * \param[in] z The Z coordinate
 ***********************************************/
void MCTruth::SetParticlePosition(float x, float y , float z)
{
	fPosX = x;
	fPosY = y;
	fPosZ = z;
	return;
}

void MCTruth::SetQuenchingFactor(double qfvalue)
{
  fQFValue = qfvalue;
  return;
}

/*!**********************************************
 * \fn void MCTruth::AddHitInformation(int row, int col, bool uncoated, float energy, float time, float angle)
 * \brief Adds new hit information to the current vectors
 * \param[in] row The row number
 * \param[in] col The column number
 * \param[in] uncoated A flag saying if the PMT is coated or uncoated
 * \param[in] energy The energy of the hit in eV
 * \param[in] time The time of the hit in ns
 * \param[in] angle The angle the PMT was hit
 ***********************************************/
void MCTruth::AddHitInformation(int row, int col, bool uncoated, float energy, float time, float angle, bool quench)
{
	fHitRow.emplace_back(row);
	fHitCol.emplace_back(col);
	fHitUncoated.push_back(uncoated);
	fHitEnergy.emplace_back(energy);
	fHitTime.emplace_back(time);
	fHitAngle.emplace_back(angle);
	fPassedQF.push_back(quench);

	return;
}

/*!**********************************************
 * \fn void MCTruth::SetHitInformation(size_t hitIndex, int row, int col, bool uncoated, float energy, float time, float angle)
 * \brief Sets the hit information at the passed hit index
 * \param[in] hitIndex The index to save the information to
 * \param[in] row The row number
 * \param[in] col The column number
 * \param[in] uncoated A flag saying if the PMT is coated or uncoated
 * \param[in] energy The energy of the hit in eV
 * \param[in] time The time of the hit in ns
 * \param[in] angle The angle the PMT was hit
 * \param[in] quench True if the event is not quenched, false if it is
 *
 * If hitIndex is greater than the current size of the vectors
 * the information will be appended to the vectors and an
 * error message will be printed
 ***********************************************/
void MCTruth::SetHitInformation(size_t hitIndex, int row, int col, bool uncoated, float energy, float time, float angle, bool quench)
{
	if (hitIndex >= fHitRow.size()) {
		MsgError(MsgLog::Form("The index passed %zu is greater than current size %zu. Appending the hit information",
					hitIndex,fHitRow.size()));
		AddHitInformation(row,col,uncoated,energy,time,angle,quench);
		return;
	}

	fHitRow.at(hitIndex) = row;
	fHitCol.at(hitIndex) = col;
	fHitUncoated.at(hitIndex) = uncoated;
	fHitEnergy.at(hitIndex) = energy;
	fHitTime.at(hitIndex) = time;
	fHitAngle.at(hitIndex) = angle;
	fPassedQF.at(hitIndex) = quench;

	return;
}


/*!**********************************************
 * \fn void MCTruth::SetHitQuench(size_t hitIndex, bool quench)
 * \brief changes the quenched value for a single hit.
 * \param[in] hitIndex The index to save the information to
 * \param[in] quench True if the event is not quenched, false if it is
 ***********************************************/
void MCTruth::SetHitQuench(size_t hitIndex, bool quench)
{
  if (fPassedQF.size() != fHitRow.size()) {
    fPassedQF.resize(fHitRow.size(),true);
  }

	if (hitIndex >= fPassedQF.size()) {
		MsgError(MsgLog::Form("The index passed %zu is greater than current size %zu. No hit will be changed.",
					hitIndex,fPassedQF.size()));
		return;
	}
  
	fPassedQF.at(hitIndex) = quench;

	return;
}

/*!**********************************************
 * \fn void MCTruth::SetNumberOfHits(size_t numHits)
 * \brief Sets the hit information vectors to the size passed
 * \param[in] numHits The size to set the vectors to
 ***********************************************/
void MCTruth::SetNumberHits(size_t numHits)
{
	fHitRow.resize(numHits);
	fHitCol.resize(numHits);
	fHitUncoated.resize(numHits);
	fHitEnergy.resize(numHits);
	fHitTime.resize(numHits);
	fHitAngle.resize(numHits);
  fPassedQF.resize(numHits,true);
	return;
}

/*!**********************************************
 * \fn void MCTruth::GetParticleMomentum(float & x, float & y, flato & z)
 * \brief Gets the particles initial momentum 
 * \param[out] x The x coordinate 
 * \param[out] y The Y coordinate
 * \param[out] z The Z coordinate
 ***********************************************/
void MCTruth::GetParticleMomentum(float & x, float & y , float & z)
{
	x = fMomX;
	y = fMomY;
	z = fMomZ;
	return;
}

/*!**********************************************
 * \fn void MCTruth::GetParticlePosition(float & x, float & y, flato & z)
 * \brief Gets the particles initial position 
 * \param[out] x The x coordinate 
 * \param[out] y The Y coordinate
 * \param[out] z The Z coordinate
 ***********************************************/
void MCTruth::GetParticlePosition(float & x, float & y , float & z)
{
	x = fPosX;
	y = fPosY;
	z = fPosZ;
	return;
}

/*!**********************************************
 * \fn void MCTruth::GetParticleMomentum(float & x, float & y, flato & z) const
 * \brief Gets the particles initial momentum 
 * \param[out] x The x coordinate 
 * \param[out] y The Y coordinate
 * \param[out] z The Z coordinate
 ***********************************************/
void MCTruth::GetParticleMomentum(float & x, float & y , float & z) const
{
	x = fMomX;
	y = fMomY;
	z = fMomZ;
	return;
}

/*!**********************************************
 * \fn void MCTruth::GetParticlePosition(float & x, float & y, flato & z) const
 * \brief Gets the particles initial position 
 * \param[out] x The x coordinate 
 * \param[out] y The Y coordinate
 * \param[out] z The Z coordinate
 ***********************************************/
void MCTruth::GetParticlePosition(float & x, float & y , float & z) const
{
	x = fPosX;
	y = fPosY;
	z = fPosZ;
	return;
}

/*!**********************************************
 * \fn MCTruth & MCTruth::operator=(const MCTruth & rhs)
 * \brief Copy assignment for the #MCTruth class
 * \return Return reference to the current class
 ***********************************************/
MCTruth & MCTruth::operator=(const MCTruth & rhs)
{

	this->fParticleID = rhs.fParticleID;
	this->fEnergyDeposited = rhs.fEnergyDeposited;
	this->fMomX = rhs.fMomX;
	this->fMomY = rhs.fMomY;
	this->fMomZ = rhs.fMomZ;
	this->fPosX = rhs.fPosX;
	this->fPosY = rhs.fPosY;
	this->fPosZ = rhs.fPosZ;

	this->fHitEnergy = rhs.fHitEnergy;
	this->fHitTime = rhs.fHitTime;
	this->fHitRow = rhs.fHitRow;
	this->fHitCol = rhs.fHitCol;
	this->fHitUncoated = rhs.fHitUncoated;
	this->fHitAngle = rhs.fHitAngle;
	this->fPassedQF = rhs.fPassedQF;

  this->fQFValue = rhs.fQFValue;

  return *this;
}



