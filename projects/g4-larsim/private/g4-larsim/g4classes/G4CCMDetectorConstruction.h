/*
  header file for the detector construction

initiates the methods used in G4CCMDetectorConstruction.cc and defines a number of global variables
Methods include: geometry modifications, pmt placement, material definition, and boolean swtichers
Global variables include most of the important volumes and materials.
*/
#ifndef G4CCMDetectorConstruction_H
#define G4CCMDetectorConstruction_H 1

#include "icetray/CCMPMTKey.h"
#include "dataclasses/I3Map.h"
#include "g4-larsim/g4classes/G4CCMDetectorMessenger.h"
#include "g4-larsim/g4classes/G4CCMReadout.h"
#include "g4-larsim/g4classes/G4CCMTreeTracker.h"
#include "g4-larsim/g4classes/G4CCMEDepSD.h"

#include <G4Cache.hh>
#include <G4Material.hh>
#include <G4RunManager.hh>
#include <G4VisAttributes.hh>
#include <G4RotationMatrix.hh>
#include <CLHEP/Units/SystemOfUnits.h>
#include <G4VUserDetectorConstruction.hh>
#include <G4SystemOfUnits.hh>

class G4CCMMainVolume;
class G4CCMPMTSD;

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4Material;
class G4MaterialPropertiesTable;
class G4Sphere;
class G4Tubs;
class G4VPhysicalVolume;

class G4CCMDetectorConstruction : public G4VUserDetectorConstruction {
  public:

    // constructor and destructor
    G4CCMDetectorConstruction(G4double UVAbsA, G4double UVAbsB, G4double UVAbsD, G4double UVAbsScaling,
                              G4double WLSNPhotonsEndCapFoil, G4double WLSNPhotonsSideFoil, G4double WLSNPhotonsPMT,
                              G4double EndCapFoilTPBThickness, G4double SideFoilTPBThickness, G4double PMTTPBThickness,
                              G4double Rayleigh128, G4double TPBAbsTau, G4double TPBAbsNorm, G4double TPBAbsScale,
                              G4double Mie_GG, G4double Mie_Ratio, G4double Normalization, G4double PhotonSamplingFactor);
    ~G4CCMDetectorConstruction() override;

    void SetReadout(G4CCMReadout * readout);

    // build the detector
    G4VPhysicalVolume* Construct() override;

    // Add SD
    void ConstructSDandField() override;

    //Functions to modify the geometry
    void SetDefaults();//Method to set the default values of all geometry modifications

    //Construct main volume
    void SetMainVolumeOn(G4bool b);
    G4bool GetMainVolumeOn() const { return fMainVolumeOn; }

    void SetSaveAllEnergyLossesVector(G4bool b) { SaveAllEnergyLossesVector_ = b; }
    void SetSaveAllEnergyLossesTree(G4bool b) { SaveAllEnergyLossesTree_ = b; }

    void SetVetoSDSaveEnergyLossesVector(G4bool b) { VetoSDSaveEnergyLossesVector_ = b; }
    void SetVetoSDSaveEnergyLossesTree(G4bool b) { VetoSDSaveEnergyLossesTree_ = b; }
    void SetVetoSDPruneTree(G4bool b) { VetoSDPruneTree_ = b; }

    void SetInteriorSDSaveEnergyLossesVector(G4bool b) { InteriorSDSaveEnergyLossesVector_ = b; }
    void SetInteriorSDSaveEnergyLossesTree(G4bool b) { InteriorSDSaveEnergyLossesTree_ = b; }
    void SetInteriorSDPruneTree(G4bool b) { InteriorSDPruneTree_ = b; }

    void SetKillNeutrinos(G4bool b) { KillNeutrinos_ = b; }
    void SetKillPhotons(G4bool b) { KillPhotons_ = b; }
    void SetKillScintillation(G4bool b) { KillScintillation_ = b; }
    void SetKillCherenkov(G4bool b) { KillCherenkov_ = b; }

    void SetTimeCut(G4bool b) { TimeCut_ = b; }
    void SetDetailedPhotonTracking(G4bool b) { DetailedPhotonTracking_ = b; }
    void SetTrackParticles(G4bool b) { TrackParticles_ = b; }
    void SetTrackEnergyLosses(G4bool b) { TrackEnergyLosses_ = b; }

    void SetG4RangeCut(G4double G4RangeCut) { G4RangeCut_ = G4RangeCut; }
    void SetG4EDepMin(G4double G4EDepMin) { G4EDepMin_ = G4EDepMin; }
    void SetG4ETrackingMin(G4double G4ETrackingMin) { G4ETrackingMin_ = G4ETrackingMin; }

    // set SD configuration
    void SetRecordHits(bool RecordHits) { RecordHits_ = RecordHits; }

    bool GetSaveAllEnergyLossesVector() { return SaveAllEnergyLossesVector_; }
    bool GetSaveAllEnergyLossesTree() { return SaveAllEnergyLossesTree_; }

    bool GetVetoSDSaveEnergyLossesVector() { return VetoSDSaveEnergyLossesVector_; }
    bool GetVetoSDSaveEnergyLossesTree() { return VetoSDSaveEnergyLossesTree_; }
    bool GetVetoSDPruneTree() { return VetoSDPruneTree_; }

    bool GetInteriorSDSaveEnergyLossesVector() { return InteriorSDSaveEnergyLossesVector_; }
    bool GetInteriorSDSaveEnergyLossesTree() { return InteriorSDSaveEnergyLossesTree_; }
    bool GetInteriorSDPruneTree() { return InteriorSDPruneTree_; }

    bool GetKillNeutrinos() { return KillNeutrinos_; }
    bool GetKillPhotons() { return KillPhotons_; }
    bool GetKillScintillation() { return KillScintillation_; }
    bool GetKillCherenkov() { return KillCherenkov_; }

    bool GetTimeCut() { return TimeCut_; }
    bool GetDetailedPhotonTracking() { return DetailedPhotonTracking_; }
    bool GetTrackParticles() { return TrackParticles_; }
    bool GetTrackEnergyLosses() { return TrackEnergyLosses_; }

    G4double GetG4RangeCut() { return G4RangeCut_; }
    G4double GetG4EDepMin() { return G4EDepMin_; }
    G4double GetG4ETrackingMin() { return G4ETrackingMin_; }

    // get SD configuration
    bool GetRecordHits() { return RecordHits_; }

    // set sodium source calibration run status
    void InitializeSodiumSourceRun(bool SourceRodIn, G4double SourceRodLocation, bool CobaltSourceRun,
                                  bool SodiumSourceRun, bool TrainingSource, G4double DecayX, G4double DecayY, G4double DecayZ) {
        SourceRodIn_ = SourceRodIn;
        SourceRodLocation_ = SourceRodLocation;
        CobaltSourceRun_ = CobaltSourceRun;
        SodiumSourceRun_ = SodiumSourceRun;
        TrainingSource_ = TrainingSource;
        DecayX_ = DecayX;
        DecayY_ = DecayY;
        DecayZ_ = DecayZ;
    }

    double HarmonicOscillatorRefractiveIndex(double a0, double aUV, double gammaUV, double wlUV, double wl);
    double HarmonicOscillatorRefractiveIndexDerivative(double aUV, double lambda, double lambdaUV, double gammaUV);

  private:

    G4CCMReadout * readout_ = nullptr;

    void DefineMaterials();
    G4CCMDetectorMessenger* fDetectorMessenger = nullptr;

    G4Box* fExperimentalHall_box = nullptr;
    G4LogicalVolume* fExperimentalHall_log = nullptr;
    G4VPhysicalVolume* fExperimentalHall_phys = nullptr;

    // Materials & Elements
    G4Element* fH = nullptr;
    G4Element* fC = nullptr;
    G4Element* fN = nullptr;
    G4Element* fO = nullptr;
    G4Element* elFe = nullptr;
    G4Element* elCr = nullptr;
    G4Element* elNi = nullptr;
    G4Element* elC = nullptr;
    G4Element* elMn = nullptr;
    G4Element* elF = nullptr;
    G4Material* fGlass = nullptr;
    G4Material* fPlastic = nullptr;
    G4Material* fBlackPlastic = nullptr;
    G4Material* fLAr = nullptr;
    G4Material* fAlum = nullptr;
    G4Material* fSteel = nullptr;
    G4Material* fVacuum = nullptr;
    G4Material* fPTFE = nullptr;
    G4Material* fTPBFoilSides = nullptr;
    G4Material* fTPBFoilTopBottom = nullptr;
    G4Material* fTPBPMT = nullptr;

    // Geometry
    G4bool fMainVolumeOn = true;

    G4CCMMainVolume* fMainVolume = nullptr;

    G4MaterialPropertiesTable* fGlass_mt = nullptr;
    G4MaterialPropertiesTable* fPlastic_mt = nullptr;
    G4MaterialPropertiesTable* fBlackPlastic_mt = nullptr;
    G4MaterialPropertiesTable* fLAr_mt = nullptr;
    G4MaterialPropertiesTable* fAlum_mt = nullptr;
    G4MaterialPropertiesTable* fSteel_mt = nullptr;
    G4MaterialPropertiesTable* fVacuum_mt = nullptr;
    G4MaterialPropertiesTable* fPTFE_mt = nullptr;
    G4MaterialPropertiesTable* fTPBFoilSides_mt = nullptr;
    G4MaterialPropertiesTable* fTPBFoilTopBottom_mt = nullptr;
    G4MaterialPropertiesTable* fTPBPMT_mt = nullptr;

    // Sensitive Detector
    G4Cache<G4CCMPMTSD*> fPMT_SD;
    G4Cache<G4CCMTreeTracker*> fTreeTracker_SD;
    G4Cache<G4CCMEDepSD*> fVeto_SD;
    G4Cache<G4CCMEDepSD*> fInterior_SD;

    bool SaveAllEnergyLossesVector_ = false;
    bool SaveAllEnergyLossesTree_ = false;
    bool VetoSDSaveEnergyLossesVector_ = false;
    bool VetoSDSaveEnergyLossesTree_ = false;
    bool VetoSDPruneTree_ = false;
    bool InteriorSDSaveEnergyLossesVector_ = false;
    bool InteriorSDSaveEnergyLossesTree_ = false;
    bool InteriorSDPruneTree_ = false;

    bool KillNeutrinos_ = false;
    bool KillPhotons_ = false;
    bool KillScintillation_ = false;
    bool KillCherenkov_ = false;

    bool TimeCut_ = false;
    bool DetailedPhotonTracking_ = false;
    bool TrackParticles_ = false;
    bool TrackEnergyLosses_ = false;



    // controls to turn SD on/off (set via G4Interface)
    bool RecordHits_ = true; // turn hit recording on/off

    G4double G4RangeCut_ = 1e-6 * m; // range cut for all particles
    G4double G4EDepMin_ = 0.01 * keV; // minimum energy deposit for hits
    G4double G4ETrackingMin_ = 0.1 * keV; // minimum energy for tracking

    // controls to turn sodium source on/off
    bool SourceRodIn_ = false;
    G4double SourceRodLocation_ = 0.0 * cm;
    bool CobaltSourceRun_ = false;
    bool SodiumSourceRun_ = false;
    bool TrainingSource_ = false;
    G4double DecayX_ = 0.0 * cm;
    G4double DecayY_ = 0.0 * cm;
    G4double DecayZ_ = 0.0 * cm;

    G4double Rayleigh128_ = 95.0 * cm;
    G4double UVAbsA_ = 0.234 * (1.0/nm);
    G4double UVAbsB_ = 113.02 * nm;
    G4double UVAbsD_ = 5.8 * cm;
    G4double UVAbsScaling_ = 1.0;
    G4double WLSNPhotonsEndCapFoil_ = 0.605;
    G4double WLSNPhotonsSideFoil_ = 0.605;
    G4double WLSNPhotonsPMT_ = 0.605;
    G4double EndCapFoilTPBThickness_ = 0.00278035 * mm;
    G4double SideFoilTPBThickness_ = 0.00278035 * mm;
    G4double PMTTPBThickness_ = 0.00203892 * mm;
    G4double TPBAbsTau_ = 0.13457;
    G4double TPBAbsNorm_ = 8.13914e-21;
    G4double TPBAbsScale_ = 1.0;
    G4double Mie_GG_ = 0.99;
    G4double Mie_Ratio_ = 1.0;
    G4double Normalization_ = 1.0;
    G4double PhotonSamplingFactor_ = 0.5;
};

#endif
