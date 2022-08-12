/*!**********************************************
 * \file MCTruth.h
 * \brief Header file for the #MCTruth class
 * \author R. T. Thornton (LANL)
 * \date November 11, 2020
 ***********************************************/
#ifndef MCTruth_h
#define MCTruth_h

#include <vector>
#include <iostream>

#include "TObject.h"

/*!**********************************************
 * \class MCTruth
 * \brief Container for the raw data coming out for the detector
 ***********************************************/
class MCTruth : public TObject
{
  public:
    MCTruth();
    MCTruth(int particleID);
    MCTruth(const MCTruth & rhs);
    ~MCTruth();

    void Reset(int particleID);

		void SetParticleID(int particleID) { fParticleID = particleID; }
		void SetTotalEnergyDeposited(float energy) { fEnergyDeposited = energy; }
		void SetParticleMomentum(float x, float y, float z);
		void SetParticlePosition(float x, float y, float z);
		void SetQuenchingFactor(double qfvalue);
		void AddHitInformation(int row, int col, bool uncoated, float energy, float time, float angle, bool quench);
		void SetHitInformation(size_t hitIndex, int row, int col, bool uncoated, float energy, float time, float angle, bool quench);
		void SetHitQuench(size_t hitIndex, bool quench);
		void SetNumberHits(size_t hitIndex);

		int GetParticleID() { return fParticleID; }
		float GetTotalEnergyDeposited() { return fEnergyDeposited; }
		float GetParticleMomentumX() { return fMomX; }
		float GetParticleMomentumY() { return fMomY; }
		float GetParticleMomentumZ() { return fMomZ; }
		float GetParticlePositionX() { return fPosX; }
		float GetParticlePositionY() { return fPosY; }
		float GetParticlePositionZ() { return fPosZ; }
		double GetQuenchingFactor() { return fQFValue; }
		void GetParticleMomentum(float & x, float & y , float & z);
		void GetParticlePosition(float & x, float & y , float & z);

		size_t GetHitNumber() { return fHitRow.size(); }
		int GetHitRow(size_t index) { return fHitRow.at(index); }
		int GetHitCol(size_t index) { return fHitCol.at(index); }
		bool GetHitUncoated(size_t index) { return fHitUncoated.at(index); }
		float GetHitEnergy(size_t index) { return fHitEnergy.at(index); }
		float GetHitAngle(size_t index) { return fHitAngle.at(index); }
		float GetHitTime(size_t index) { return fHitTime.at(index); }
    bool GetPassedQF(size_t index) { return fPassedQF.at(index); }

		int GetParticleID() const { return fParticleID; }
		float GetTotalEnergyDeposited() const { return fEnergyDeposited; }
		float GetParticleMomentumX() const { return fMomX; }
		float GetParticleMomentumY() const { return fMomY; }
		float GetParticleMomentumZ() const { return fMomZ; }
		float GetParticlePositionX() const { return fPosX; }
		float GetParticlePositionY() const { return fPosY; }
		float GetParticlePositionZ() const { return fPosZ; }
		double GetQuenchingFactor() const { return fQFValue; }
		void GetParticleMomentum(float & x, float & y , float & z) const;
		void GetParticlePosition(float & x, float & y , float & z) const;

		size_t GetHitNumber() const { return fHitRow.size(); }
		int GetHitRow(size_t index) const { return fHitRow.at(index); }
		int GetHitCol(size_t index) const { return fHitCol.at(index); }
		bool GetHitUncoated(size_t index) const { return fHitUncoated.at(index); }
		float GetHitEnergy(size_t index) const { return fHitEnergy.at(index); }
		float GetHitAngle(size_t index) const { return fHitAngle.at(index); }
		float GetHitTime(size_t index) const { return fHitTime.at(index); }
    bool GetPassedQF(size_t index) const { return fPassedQF.at(index); }

    MCTruth & operator=(const MCTruth & rhs);

  private:

		/// initial particle id using pdg encoding
		int fParticleID;

    /// Total Energy Deposited in keV
    float fEnergyDeposited;
    /// X Momentum
    float fMomX;
    /// Y Momentum
    float fMomY;
    /// Z Momentum
    float fMomZ;
    /// X Position
    float fPosX;
    /// Y Position
    float fPosY;
    /// Z Position
    float fPosZ;

		/// The energy of every hit in eV
		std::vector<float> fHitEnergy;
		/// The time of every hit in ns
		std::vector<float> fHitTime;
		/// The row number of the PMT that was hit
		std::vector<int> fHitRow;
		/// The column number of the PMT that was hit
		std::vector<int> fHitCol;
		/// A flag if the PMT is coated or uncoated
		std::vector<bool> fHitUncoated;
		/// The angle at which the PMT was hit
		std::vector<float> fHitAngle;
		
		std::vector<bool> fPassedQF;
		double fQFValue;


  ClassDef(MCTruth,2)
};

#endif // MCTruth_h

