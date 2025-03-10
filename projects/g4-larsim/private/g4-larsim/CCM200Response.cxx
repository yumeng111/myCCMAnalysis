// standard library stuff

#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include "g4-larsim/CCM200Response.h"
#include "g4-larsim/CCMDetectorResponse.h"
#include "g4-larsim/g4classes/G4Interface.h"

#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3FrameObject.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3ConditionalModule.h"
#include "icetray/I3SingleServiceFactory.h"

#include "phys-services/I3RandomService.h"
#include "simclasses/CCMMCPE.h"

#include <vector>
#include <string>
#include <algorithm>

    bool KillNeutrinos_ = false;
    bool KillPhotons_ = false;
    bool KillScintillation_ = false;
    bool KillCherenkov_ = false;

    bool TimeCut_ = false;
    bool DetailedPhotonTracking_ = false;
    bool TrackParticles_ = false;
    bool TrackEnergyLosses_ = false;

CCM200Response::CCM200Response(const I3Context& context) :
    CCMDetectorResponse(context),
    SaveAllEnergyLossesTree_(false),
    VetoSDSaveEnergyLossesVector_(false), VetoSDSaveEnergyLossesTree_(false), VetoSDPruneTree_(false),
    InteriorSDSaveEnergyLossesVector_(false), InteriorSDSaveEnergyLossesTree_(false), InteriorSDPruneTree_(false),
    KillNeutrinos_(false), KillPhotons_(false), KillScintillation_(false), KillCherenkov_(false),
    TimeCut_(false), DetailedPhotonTracking_(false), TrackParticles_(false), TrackEnergyLosses_(false),
    RecordHits_(true), SourceRodIn_(false), SourceRodLocation_(0.0 * I3Units::cm),
    CobaltSourceRun_(false), SodiumSourceRun_(false), TrainingSource_(false), DecayX_(0.0 * I3Units::cm), DecayY_(0.0 * I3Units::cm), DecayZ_(0.0 * I3Units::cm),
    SingletTau_(8.2 * I3Units::nanosecond), TripletTau_(743.0 * I3Units::nanosecond),
    Rayleigh128_(95.0 * I3Units::cm), UVAbsA_(0.234 * (1.0/I3Units::nanometer)), UVAbsB_(113.02 * I3Units::nanometer), UVAbsD_(5.8 * I3Units::cm), UVAbsScaling_(1.0),
    WLSNPhotonsEndCapFoil_(0.605), WLSNPhotonsSideFoil_(0.605), WLSNPhotonsPMT_(0.605),
    EndCapFoilTPBThickness_(0.00278035 * I3Units::mm), SideFoilTPBThickness_(0.00278035 * I3Units::mm), PMTTPBThickness_(0.00203892 * I3Units::mm),
    TPBAbsTau_(0.13457), TPBAbsNorm_(8.13914e-21), TPBAbsScale_(1.0), MieGG_(0.99), MieRatio_(0.8), Normalization_(1.0), PhotonSampling_(0.5), RandomSeed_(0) {

    AddParameter("SaveAllEnergyLossesTree", "save all energy losses to tree", SaveAllEnergyLossesTree_);
    AddParameter("VetoSDSaveEnergyLossesVector", "save energy losses in veto sensitive detector to vector", VetoSDSaveEnergyLossesVector_);
    AddParameter("VetoSDSaveEnergyLossesTree", "save energy losses in veto sensitive detector to tree", VetoSDSaveEnergyLossesTree_);
    AddParameter("VetoSDPruneTree", "prune tree in veto sensitive detector", VetoSDPruneTree_);
    AddParameter("InteriorSDSaveEnergyLossesVector", "save energy losses in interior sensitive detector to vector", InteriorSDSaveEnergyLossesVector_);
    AddParameter("InteriorSDSaveEnergyLossesTree", "save energy losses in interior sensitive detector to tree", InteriorSDSaveEnergyLossesTree_);
    AddParameter("InteriorSDPruneTree", "prune tree in interior sensitive detector", InteriorSDPruneTree_);

    AddParameter("KillNeutrinos", "kill neutrinos", KillNeutrinos_);
    AddParameter("KillPhotons", "kill photons", KillPhotons_);
    AddParameter("KillScintillation", "kill scintillation", KillScintillation_);
    AddParameter("KillCherenkov", "kill cherenkov", KillCherenkov_);
    AddParameter("TimeCut", "only track events up to 200nsec", TimeCut_);
    AddParameter("DetailedPhotonTracking", "track all optical photons to get distance travelled", DetailedPhotonTracking_);
    AddParameter("TrackParticles", "track particles", TrackParticles_);
    AddParameter("TrackEnergyLosses", "track energy losses", TrackEnergyLosses_);

    AddParameter("RecordHits", "record hits", RecordHits_);
    AddParameter("SourceRodIn", "true if we want to simulate the sodium source rod", SourceRodIn_);
    AddParameter("SourceRodLocation", "z location of the end of the sodium source rod", SourceRodLocation_);
    AddParameter("CobaltSourceRun", "true if we want to simulate cobalt source pellet", CobaltSourceRun_);
    AddParameter("SodiumSourceRun", "true if we want to simulate sodium source pellet", SodiumSourceRun_);
    AddParameter("TrainingSource", "true if we want to simulate training source events", TrainingSource_);
    AddParameter("DecayX", "if generating training source data, provid X position", DecayX_);
    AddParameter("DecayY", "if generating training source data, provid Y position", DecayY_);
    AddParameter("DecayZ", "if generating training source data, provid Z position", DecayZ_);
    AddParameter("SingletTimeConstant", "LAr singlet tau", SingletTau_);
    AddParameter("TripletTimeConstant", "LAr triplet tau", TripletTau_);
    AddParameter("Rayleigh128Length", "Rayleigh scattering length for 128nm light", Rayleigh128_);
    AddParameter("UVAbsA", "Set UV absorption slope [1/nm]", UVAbsA_ / (1.0/I3Units::nanometer));
    AddParameter("UVAbsB", "Set UV absorption offset [nm]", UVAbsB_ / I3Units::nanometer);
    AddParameter("UVAbsD", "Set UV absorption reference distance [m]", UVAbsD_ / I3Units::meter);
    AddParameter("UVAbsScaling", "Set UV absorption scale [dimensionless]", UVAbsScaling_);
    AddParameter("WLSNPhotonsEndCapFoil", "mean number of photons produced per WLS for TPB foils on the end caps of the detector", WLSNPhotonsEndCapFoil_);
    AddParameter("WLSNPhotonsSideFoil", "mean number of photons produced per WLS for TPB foils on the sides of the detector", WLSNPhotonsSideFoil_);
    AddParameter("WLSNPhotonsPMT", "mean number of photons produced per WLS for TPB on PMTs", WLSNPhotonsPMT_);
    AddParameter("EndCapFoilTPBThickness", "thickness of TPB on the endcap foil", EndCapFoilTPBThickness_);
    AddParameter("SideFoilTPBThickness", "thickness of TPB on the side foil", SideFoilTPBThickness_);
    AddParameter("PMTTPBThickness", "thickness of TPB on the PMTs", PMTTPBThickness_);
    AddParameter("TPBAbsorptionTau", "factor in exponential for visible part of tpb absorption spectrum", TPBAbsTau_);
    AddParameter("TPBAbsorptionNorm", "normalization for exponential for visible part of tpb absorption spectrum", TPBAbsNorm_);
    AddParameter("TPBAbsorptionScale", "overall scaling for entire tpb absorption curve", TPBAbsScale_);
    AddParameter("MieGG", "used to calculate angular spread of outgoing photons from mie scattering in tpb", MieGG_);
    AddParameter("MieRatio", "ratio of forward : backwards scattering from mie scattering in tpb", MieRatio_);
    AddParameter("Normalization", "normalization of number of photons produced in liquid argon", Normalization_);
    AddParameter("PhotonSampling", "scaling of number of photons produced in liquid argon", PhotonSampling_);
    AddParameter("RandomSeed", "seed for geant4 random generator", RandomSeed_);
}

void CCM200Response::Configure() {
    GetParameter("SaveAllEnergyLossesTree", SaveAllEnergyLossesTree_);
    GetParameter("VetoSDSaveEnergyLossesVector", VetoSDSaveEnergyLossesVector_);
    GetParameter("VetoSDSaveEnergyLossesTree", VetoSDSaveEnergyLossesTree_);
    GetParameter("VetoSDPruneTree", VetoSDPruneTree_);
    GetParameter("InteriorSDSaveEnergyLossesVector", InteriorSDSaveEnergyLossesVector_);
    GetParameter("InteriorSDSaveEnergyLossesTree", InteriorSDSaveEnergyLossesTree_);
    GetParameter("InteriorSDPruneTree", InteriorSDPruneTree_);

    GetParameter("KillNeutrinos", KillNeutrinos_);
    GetParameter("KillPhotons", KillPhotons_);
    GetParameter("KillScintillation", KillScintillation_);
    GetParameter("KillCherenkov", KillCherenkov_);
    GetParameter("TimeCut", TimeCut_);
    GetParameter("DetailedPhotonTracking", DetailedPhotonTracking_);
    GetParameter("TrackParticles", TrackParticles_);
    GetParameter("TrackEnergyLosses", TrackEnergyLosses_);

    GetParameter("RecordHits", RecordHits_);
    GetParameter("SourceRodIn", SourceRodIn_);
    GetParameter("SourceRodLocation", SourceRodLocation_);
    GetParameter("CobaltSourceRun", CobaltSourceRun_);
    GetParameter("SodiumSourceRun", SodiumSourceRun_);
    GetParameter("TrainingSource", TrainingSource_);
    GetParameter("DecayX", DecayX_);
    GetParameter("DecayY", DecayY_);
    GetParameter("DecayZ", DecayZ_);
    GetParameter("SingletTimeConstant", SingletTau_);
    GetParameter("TripletTimeConstant", TripletTau_);
    GetParameter("Rayleigh128Length", Rayleigh128_);
    GetParameter("UVAbsA", UVAbsA_);
    GetParameter("UVAbsB", UVAbsB_);
    GetParameter("UVAbsD", UVAbsD_);
    GetParameter("UVAbsScaling", UVAbsScaling_);
    GetParameter("WLSNPhotonsEndCapFoil", WLSNPhotonsEndCapFoil_);
    GetParameter("WLSNPhotonsSideFoil", WLSNPhotonsSideFoil_);
    GetParameter("WLSNPhotonsPMT", WLSNPhotonsPMT_);
    GetParameter("EndCapFoilTPBThickness", EndCapFoilTPBThickness_);
    GetParameter("SideFoilTPBThickness", SideFoilTPBThickness_);
    GetParameter("PMTTPBThickness", PMTTPBThickness_);
    GetParameter("TPBAbsorptionTau", TPBAbsTau_);
    GetParameter("TPBAbsorptionNorm", TPBAbsNorm_);
    GetParameter("TPBAbsorptionScale", TPBAbsScale_);
    GetParameter("MieGG", MieGG_);
    GetParameter("MieRatio", MieRatio_);
    GetParameter("Normalization", Normalization_);
    GetParameter("PhotonSampling", PhotonSampling_);
    GetParameter("RandomSeed", RandomSeed_);

    UVAbsA_ *= (1.0/I3Units::nanometer);
    UVAbsB_ *= I3Units::nanometer;
    UVAbsD_ *= I3Units::meter;
}

CCM200Response::~CCM200Response() {}

void CCM200Response::Initialize() {

    if(g4Interface_ == nullptr) {
        g4Interface_ = G4Interface::GetInstance();
    }

    // let's let's construct the detector
    g4Interface_->InstallDetector(
                                  SaveAllEnergyLossesTree_,
                                  VetoSDSaveEnergyLossesVector_, VetoSDSaveEnergyLossesTree_, VetoSDPruneTree_,
                                  InteriorSDSaveEnergyLossesVector_, InteriorSDSaveEnergyLossesTree_, InteriorSDPruneTree_,
                                  KillNeutrinos_, KillPhotons_, KillScintillation_, KillCherenkov_,
                                  TimeCut_, DetailedPhotonTracking_, TrackParticles_, TrackEnergyLosses_,
                                  RecordHits_, SourceRodIn_, SourceRodLocation_, CobaltSourceRun_, SodiumSourceRun_, TrainingSource_, 
                                  DecayX_, DecayY_, DecayZ_, SingletTau_, TripletTau_, Rayleigh128_, UVAbsA_, UVAbsB_, UVAbsD_, UVAbsScaling_,
                                  WLSNPhotonsEndCapFoil_, WLSNPhotonsSideFoil_, WLSNPhotonsPMT_,
                                  EndCapFoilTPBThickness_, SideFoilTPBThickness_, PMTTPBThickness_, TPBAbsTau_, TPBAbsNorm_, TPBAbsScale_,
                                  MieGG_, MieRatio_, Normalization_, PhotonSampling_, RandomSeed_);
}

I3FrameObjectPtr CCM200Response::GetSimulationConfiguration() {
    DetectorResponseConfigPtr config = boost::make_shared<DetectorResponseConfig>();
    config->rayleigh_scattering_length_ = Rayleigh128_;
    config->uv_absorption_a_ = UVAbsA_;
    config->uv_absorption_b_ = UVAbsB_;
    config->uv_absorption_d_ = UVAbsD_;
    config->uv_absorption_scaling_ = UVAbsScaling_;
    config->pmt_tpb_qe_ = WLSNPhotonsPMT_;
    config->endcap_tpb_qe_ = WLSNPhotonsEndCapFoil_;
    config->side_tpb_qe_ = WLSNPhotonsSideFoil_;
    config->pmt_tpb_thickness_ = PMTTPBThickness_;
    config->endcap_tpb_thickness_ = EndCapFoilTPBThickness_;
    config->side_tpb_thickness_ = SideFoilTPBThickness_;
    config->tpb_abs_tau_ = TPBAbsTau_;
    config->tpb_abs_norm_ = TPBAbsNorm_;
    config->tpb_abs_scale_ = TPBAbsScale_;
    config->mie_gg_ = MieGG_;
    config->mie_ratio_ = MieRatio_;
    config->normalization_ = Normalization_;
    config->photon_sampling_factor_ = PhotonSampling_;
    return config;
}

void CCM200Response::SimulateEvent(const I3Particle& primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries, I3MCTreePtr veto_tree, I3MCTreePtr inner_tree, I3VectorI3ParticlePtr veto_vector, I3VectorI3ParticlePtr inner_vector) {
    g4Interface_->SimulateEvent(primary, tree, mcpeseries, veto_tree, inner_tree, veto_vector, inner_vector);
}

void CCM200Response::SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> veto_trees, std::vector<I3MCTreePtr> inner_trees, std::vector<I3VectorI3ParticlePtr> veto_vectors, std::vector<I3VectorI3ParticlePtr> inner_vectors) {
    g4Interface_->SimulateEvents(primaries, trees, mcpeseries, veto_trees, inner_trees, veto_vectors, inner_vectors);
}

void CCM200Response::DestroyInterface() {

    g4Interface_->DestroyInstance();
}

typedef I3SingleServiceFactory<CCM200Response,CCMDetectorResponse> CCM200ResponseFactory;

I3_SERVICE_FACTORY(CCM200ResponseFactory);


