
#include "g4-larsim/g4classes/G4CCMEDepSD.h"
#include "g4-larsim/g4classes/G4CCMScintHit.h"

#include <G4ios.hh>
#include <G4Step.hh>
#include <G4Event.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>
#include <G4SDManager.hh>
#include <G4VTouchable.hh>
#include <G4ParticleTypes.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4TouchableHistory.hh>
#include <G4ParticleDefinition.hh>
#include <G4EventManager.hh>
#include <G4TrackingManager.hh>


G4CCMEDepSD::G4CCMEDepSD(G4String name, G4CCMReadout::VolumeType volume) : G4VSensitiveDetector(name), volume_(volume) {
    collectionName.insert(name + "ScintCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMEDepSD::Initialize(G4HCofThisEvent* hitsCE) {
    fScintCollection = new G4CCMScintHitsCollection(SensitiveDetectorName, collectionName[0]);

    if(fHitsCID < 0) {
        fHitsCID = G4SDManager::GetSDMpointer()->GetCollectionID(fScintCollection);
    }
    hitsCE->AddHitsCollection(fHitsCID, fScintCollection);

    event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    primary_ = readout_->GetPrimary(event_id);

    if(SaveEnergyLossesVector_) {
        output_energy_losses_vector = readout_->GetVolumeEDepVector(event_id, volume_);
    }
    if(SaveEnergyLossesTree_) {
        output_energy_losses_tree = readout_->GetVolumeEDepMCTree(event_id, volume_);
        if(not I3MCTreeUtils::Has(*output_energy_losses_tree, primary_.GetID())) {
            I3MCTreeUtils::AddPrimary(*output_energy_losses_tree, primary_);
        }
    }
}

void G4CCMEDepSD::EndOfEvent(G4HCofThisEvent*) {
    if(SaveEnergyLossesTree_) {
        std::map<I3ParticleID, std::tuple<bool, bool, std::vector<I3ParticleID>>> energy_loss_map;
        std::vector<I3ParticleID> queue;
        std::vector<I3Particle*> primaries;
        if(tree_tracker->mcTree->size() > 0)
            primaries = I3MCTreeUtils::GetPrimariesPtr(tree_tracker->mcTree);

        // Add primaries to the queue
        for(const I3Particle* primary : primaries) {
            I3ParticleID primary_id = primary->GetID();
            bool is_loss = G4CCMTreeTracker::energyLossPDGCodes.find(primary->GetPdgEncoding()) != G4CCMTreeTracker::energyLossPDGCodes.cend();
            bool our_loss = energy_loss_ids.find(primary_id) != energy_loss_ids.cend();
            energy_loss_map[primary_id] = std::make_tuple(is_loss, our_loss, std::vector<I3ParticleID>());
            queue.push_back(primary_id);
        }

        // Iteratively process particles in the queue and add their daughters to the queue
        while(!queue.empty()) {
            I3ParticleID particle_id = queue.back();
            const I3Particle* particle = I3MCTreeUtils::GetParticlePtr(tree_tracker->mcTree, particle_id);
            queue.pop_back();
            std::vector<I3Particle*> const daughters = I3MCTreeUtils::GetDaughtersPtr(tree_tracker->mcTree, particle_id);
            for(const I3Particle* daughter : daughters) {
                I3ParticleID daughter_id = daughter->GetID();
                bool is_loss = G4CCMTreeTracker::energyLossPDGCodes.find(daughter->GetPdgEncoding()) != G4CCMTreeTracker::energyLossPDGCodes.cend();
                bool our_loss = energy_loss_ids.find(daughter_id) != energy_loss_ids.cend();
                energy_loss_map.insert(std::make_pair(daughter_id, std::make_tuple(is_loss, our_loss, std::vector<I3ParticleID>({particle_id}))));
                queue.push_back(daughter_id);

                // if we have an energy loss then we need to mark all parents as having an energy loss
                if(our_loss) {
                    I3ParticleID parent_id = std::get<2>(energy_loss_map[daughter_id]).back();
                    while(true) {
                        std::tuple<bool, bool, std::vector<I3ParticleID>> & parent = energy_loss_map[parent_id];
                        // if we've already marked this parent as having an energy loss then we can break
                        if(std::get<1>(parent)) {
                            break;
                        }
                        // mark parent as having an energy loss
                        std::get<1>(parent) = true;
                        // if parent has no parents then we can break
                        if(std::get<2>(parent).empty()) {
                            break;
                        }
                        parent_id = std::get<2>(parent).back();
                    }
                }
            }
        }

        // Now we can build the tree
        // Start by adding the primaries
        for(const I3Particle* primary : primaries) {
            I3ParticleID primary_id = primary->GetID();
            std::tuple<bool, bool, std::vector<I3ParticleID>> & primary_energy_loss = energy_loss_map[primary_id];
            // Keep the particle if it has a relevant energy loss as a descendant
            // Or keep the particle if it is not an energy loss but we are not pruning the tree
            // When pruning the tree we discard all particles not related to the relevant energy losses
            if(std::get<1>(primary_energy_loss) or ((not PruneTree_) and (not std::get<0>(primary_energy_loss)))) {
                if(not I3MCTreeUtils::Has(*output_energy_losses_tree, primary_id))
                    I3MCTreeUtils::AddPrimary(*output_energy_losses_tree, *primary);
                queue.push_back(primary_id);
            }
        }

        // Now add the daughters
        while(!queue.empty()) {
            I3ParticleID particle_id = queue.back();
            const I3Particle* particle = I3MCTreeUtils::GetParticlePtr(tree_tracker->mcTree, particle_id);
            queue.pop_back();
            std::vector<I3Particle*> daughters = I3MCTreeUtils::GetDaughtersPtr(tree_tracker->mcTree, particle_id);
            for(const I3Particle* daughter : daughters) {
                I3ParticleID daughter_id = daughter->GetID();
                std::tuple<bool, bool, std::vector<I3ParticleID>> & daughter_energy_loss = energy_loss_map[daughter_id];
                // Keep the particle if it has a relevant energy loss as a descendant
                // Or keep the particle if it is not an energy loss but we are not pruning the tree
                // When pruning the tree we discard all particles not related to the relevant energy losses
                if(std::get<1>(daughter_energy_loss) or ((not PruneTree_) and (not std::get<0>(daughter_energy_loss)))) {
                    if(not I3MCTreeUtils::Has(*output_energy_losses_tree, daughter_id))
                        I3MCTreeUtils::AppendChild(*output_energy_losses_tree, particle_id, *daughter);
                    queue.push_back(daughter_id);
                }
            }
        }
    }

    if(SaveEnergyLossesVector_) {
        bool source_from_tree = SaveEnergyLossesTree_ or (SaveEnergyLossesVector_ and tree_tracker->GetTrackEnergyLosses());
        // Only need to fill the vector if we are referencing ids from the tree
        if(source_from_tree) {
            for(auto const & energy_loss_id : energy_loss_ids) {
                output_energy_losses_vector->push_back(I3MCTreeUtils::GetParticle(*tree_tracker->mcTree, energy_loss_id));
            }
        }
        // Otherwise vector will already be filled
    }

    Reset();
}

G4bool G4CCMEDepSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    // note -- this chunk of code resets global time to 0 in the case of radioactive decays
    // very important for retaining time structure of scintillation photons!!!
    G4Track* track = aStep->GetTrack();

    if(track->GetTrackStatus() == fStopAndKill) {
        return false;
    }

    // Check the type of particle
    G4ParticleDefinition * particle_definition = aStep->GetTrack()->GetDefinition();

    if(particle_definition == G4OpticalPhoton::OpticalPhotonDefinition()) {
        return false;
    }

    // get pdg code
    G4int pdg = particle_definition->GetPDGEncoding();

    G4int parent_id = aStep->GetTrack()->GetParentID();
    G4int track_id = aStep->GetTrack()->GetTrackID();

    // now we want to grab energy deposited, location, direction, time, and process type to save to MCTree

    bool source_from_tree = SaveEnergyLossesTree_ or (SaveEnergyLossesVector_ and tree_tracker->GetTrackEnergyLosses());

    if(source_from_tree) {
        auto it = tree_tracker->DaughterParticleMap.find(track_id);
        if(it != tree_tracker->DaughterParticleMap.cend()) {
            std::vector<I3Particle> daughters = I3MCTreeUtils::GetDaughters(*(tree_tracker->mcTree), it->second);
            if(daughters.empty()) {
                std::cout << "Warning: no daughters found for track_id " << track_id << std::endl;
                return false;
            }
            I3Particle energy_loss = *std::max_element(daughters.begin(), daughters.end(), [](I3Particle const & p0, I3Particle const & p1) -> bool {return p0.GetID() < p1.GetID();});
            if(SaveEnergyLossesTree_) {
                energy_loss_ids.insert(energy_loss.GetID());
            }
            if(SaveEnergyLossesVector_) {
                output_energy_losses_vector->push_back(energy_loss);
            }
        }
        return false;
    }

    // now let's check energy deposited
    G4double edep = aStep->GetTotalEnergyDeposit() / electronvolt * I3Units::eV;
    G4double ekin = aStep->GetTrack()->GetKineticEnergy() / electronvolt * I3Units::eV;

    // position
    G4ThreeVector prePosition = aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector postPosition = aStep->GetPostStepPoint()->GetPosition();
    I3Position position(prePosition.x() / mm * I3Units::mm, prePosition.y() / mm * I3Units::mm, prePosition.z() / mm * I3Units::mm);
    position += I3Position(postPosition.x() / mm * I3Units::mm, postPosition.y() / mm * I3Units::mm, postPosition.z() / mm * I3Units::mm);
    position /= 2.;

    // direction
    G4ThreeVector photonDirection = aStep->GetPostStepPoint()->GetMomentumDirection();
    I3Direction direction(photonDirection.x(), photonDirection.y(), photonDirection.z());

    // time
    G4double time = aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;

    // process name -- use for parent id == 0!
    G4VProcess const * process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    std::string processName = (process) ? process->GetProcessName() : "Unknown";

    // get name and pdg code
    G4ParticleDefinition* fParticleDefinition = aStep->GetTrack()->GetDefinition();
    G4String particleName = fParticleDefinition->GetParticleName();

    if(G4CCMTreeTracker::energyLossToI3ParticlePDGCode.find(processName) != G4CCMTreeTracker::energyLossToI3ParticlePDGCode.end()) {
        // now add energy loss
        I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(G4CCMTreeTracker::energyLossToI3ParticlePDGCode.at(processName));
        I3Particle daughter(daughter_type);
        daughter.SetEnergy(edep);
        daughter.SetPos(position);
        daughter.SetDir(direction);
        daughter.SetTime(time);
        output_energy_losses_vector->push_back(daughter);
    } else {
        std::cout << "oops! no conversion for " << processName << std::endl;
    }

    // now back to scint hit things
    G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
    auto theTouchable = (G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();

    G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

    // Get the average position of the hit
    G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
    pos /= 2.;

    auto scintHit = new G4CCMScintHit(thePrePV);

    scintHit->SetEdep(edep);
    scintHit->SetPos(pos);

    fScintCollection->insert(scintHit);

    return true;
}

void G4CCMEDepSD::Reset() {
    energy_loss_ids.clear();
    // Don't need the output objects anymore
    output_energy_losses_tree = nullptr;
    output_energy_losses_vector = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
