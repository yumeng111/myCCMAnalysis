
#include "g4-larsim/g4classes/G4CCMMainVolume.h"

#include <map>
#include <set>
#include <math.h>
#include <tuple>

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <globals.hh>
#include <G4Colour.hh>
#include <G4Sphere.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>
#include <G4SystemOfUnits.hh>
#include <G4LogicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4SubtractionSolid.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMMainVolume::G4CCMMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                                 G4LogicalVolume* pMotherLogical, G4bool pMany,
                                 G4int pCopyNo, G4CCMDetectorConstruction* c,
                                 G4bool SourceRodIn, G4double SourceRodLocation, G4bool CobaltSourceRun, G4bool SodiumSourceRun,
                                 G4bool TrainingSource, G4double DecayX, G4double DecayY, G4double DecayZ,
                                 G4double EndCapFoilTPBThickness, G4double SideFoilTPBThickness, G4double PMTTPBThickness)
  // Pass info to the G4PVPlacement constructor
  : G4PVPlacement(pRot, tlate,
            new G4LogicalVolume(
                new G4Tubs("Cryogen", 0*cm, 138*cm, 131*cm, 0*deg,360*deg),
                G4Material::GetMaterial("Steel"), "Cryogen"),
            "Cryogen", pMotherLogical, false, 0, true)
  , fConstructor(c) {

    // now let's build our detector
    // Outer cryogen
    fCryoVessel_log = GetLogicalVolume();

    // 3cm of steel until we get to vacuum

    // Vacuum jacket
    fVacuum = new G4Tubs("Vacuum", 0*cm, 135*cm, 126*cm, 0*deg, 360*deg);
    fVacuum_log = new G4LogicalVolume(fVacuum, G4Material::GetMaterial("Vacuum"), "Vacuum");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fVacuum_log, "Vacuum", fCryoVessel_log, false, 0, true);

    // 10cm of vacuum until we get to inner steel jacket

    // Inner cryogen
    fInnerJacket = new G4Tubs("InnerJacket", 0*cm, 125*cm, 120*cm, 0*deg, 360*deg);
    fInnerJacket_log = new G4LogicalVolume(fInnerJacket, G4Material::GetMaterial("Steel"), "InnerJacket");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fInnerJacket_log, "InnerJacket", fVacuum_log, false, 0, true);

    // 5cm of steel until we get to liquid argon

    // Argon outside the fiducial volume
    fOuterLAr = new G4Tubs("OuterLiquidArgon", 0*cm, 120*cm, 115*cm, 0*deg, 360*deg);
    fOuterLAr_log = new G4LogicalVolume(fOuterLAr, G4Material::GetMaterial("LAr"), "OuterLiquidArgon");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fOuterLAr_log, "OuterLiquidArgon", fInnerJacket_log, false, 0, true);

    // Radius of inner frame is 1037mm, inner radius is 1036mm
    // Total height of inner frame is 1239.6mm (half height is 619.8mm)
    G4double frame_thickness = 1.0 * mm;
    G4double ptfe_thickness = 0.5 * mm;
    //G4double tpb_thickness = 0.00278035 * mm; // from Vincent Basque's thesis https://pure.manchester.ac.uk/ws/portalfiles/portal/205622566/FULL_TEXT.PDF

    std::cout << "endcap tpb thickness = " << EndCapFoilTPBThickness << ", side foil tpb thickness = " << SideFoilTPBThickness << ", and pmt tpb thickness = " << PMTTPBThickness << std::endl;

    G4double frame_height = 1239.6*mm + (frame_thickness + ptfe_thickness + EndCapFoilTPBThickness) * 2.0;
    G4double frame_half_height = frame_height / 2.0;
    G4double frame_radius = 1037.0*mm - 3.15*mm + (frame_thickness + ptfe_thickness + SideFoilTPBThickness); // reducing radius by ~1/2cm to account for flexing of PTFE sheets

    // Aluminum frame holding PMTs and instrumentation
    fInnerFrame = new G4Tubs("InnerFrame", 0*cm, frame_radius, frame_half_height, 0*deg, 360*deg);
    fInnerFrame_log= new G4LogicalVolume(fInnerFrame, G4Material::GetMaterial("Alum"), "InnerFrame");
    new G4PVPlacement(0, G4ThreeVector(0,0,0), fInnerFrame_log, "InnerFrame", fOuterLAr_log, false, 0, true);

    G4double ptfe_half_height = frame_half_height - frame_thickness;
    G4double ptfe_radius = frame_radius - frame_thickness;

    // Reflector foils
    fReflectorFoil = new G4Tubs("PTFEFoil", 0*cm, ptfe_radius, ptfe_half_height, 0*deg, 360*deg);
    fReflectorFoil_log = new G4LogicalVolume(fReflectorFoil, G4Material::GetMaterial("PTFE"), "PTFEFoil");
    fReflectorFoil_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0), fReflectorFoil_log, "PTFEFoil", fInnerFrame_log, false, 0, true);

    G4double tpb_half_height = ptfe_half_height - ptfe_thickness;
    G4double tpb_radius = ptfe_radius - ptfe_thickness;

    // TPB on the foils
    // note -- trying to make TPB side foils different material than TPB top/bottom foils
    // Create the side cylinder
    fTPBFoilSides = new G4Tubs("TPBFoilSides", 0*cm, tpb_radius, tpb_half_height, 0.*deg, 360.*deg);
    fTPBFoilSides_log = new G4LogicalVolume(fTPBFoilSides, G4Material::GetMaterial("TPBFoilSides"), "TPBFoilSidesLogical");

    // Create the top and bottom disks
    fTPBFoilTop = new G4Tubs("TPBFoilTop", 0*cm, tpb_radius, EndCapFoilTPBThickness/2.0, 0.*deg, 360.*deg);
    fTPBFoilBottom = new G4Tubs("TPBFoilBottom", 0*cm, tpb_radius, EndCapFoilTPBThickness/2.0, 0.*deg, 360.*deg);
    fTPBFoilTop_log = new G4LogicalVolume(fTPBFoilTop, G4Material::GetMaterial("TPBFoilTopBottom"), "TPBFoilTopLogical");
    fTPBFoilBottom_log = new G4LogicalVolume(fTPBFoilBottom, G4Material::GetMaterial("TPBFoilTopBottom"), "TPBFoilBottomLogical");

    // Place the side cylinder in the world volume
    fTPBFoilSides_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0), fTPBFoilSides_log, "TPBFoilSides", fReflectorFoil_log, false, 0, true);

    // Place the top and bottom disks
    G4double zTopPosition = tpb_half_height - EndCapFoilTPBThickness/2;
    G4double zBottomPosition = -tpb_half_height + EndCapFoilTPBThickness/2;

    fTPBFoilTop_phys = new G4PVPlacement(0, G4ThreeVector(0,0,zTopPosition), fTPBFoilTop_log, "TPBFoilTop", fTPBFoilSides_log, false, 0, true);
    fTPBFoilBottom_phys = new G4PVPlacement(0, G4ThreeVector(0,0,zBottomPosition), fTPBFoilBottom_log, "TPBFoilBottom", fTPBFoilSides_log, false, 0, true);

    //fTPBFoil = new G4Tubs("TPBFoil", 0*cm, tpb_radius, tpb_half_height, 0*deg, 360*deg);
    //fTPBFoil_log = new G4LogicalVolume(fTPBFoil, G4Material::GetMaterial("TPBFoil"), "TPBFoil");
    //fTPBFoil_phys = new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fTPBFoil_log, "TPBFoil", fReflectorFoil_log, false, 0, true);

    G4double fiducial_lar_half_height = tpb_half_height - EndCapFoilTPBThickness;
    G4double fiducial_lar_radius = tpb_radius - SideFoilTPBThickness;

    // now fiducial LAr!
    fFiducialLAr = new G4Tubs("FiducialArgon", 0*cm, fiducial_lar_radius, fiducial_lar_half_height, 0*deg, 360*deg);
    fFiducialLAr_log = new G4LogicalVolume(fFiducialLAr, G4Material::GetMaterial("LAr"), "FiducialArgon");
    //fFiducialLAr_phys = new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fFiducialLAr_log, "FiducialArgon", fTPBFoil_log, false, 0, true);
    fFiducialLAr_phys = new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fFiducialLAr_log, "FiducialArgon", fTPBFoilSides_log, false, 0, true);

    G4double pmt_protrusion_distance = 61.89 * mm;
    G4double bridle_width = 12.0 * mm;
    G4double bridle_radius = 209.55 / 2.0 * mm;
    G4double frill_width = 0.813 * mm;
    G4double frill_radius = 228.6 / 2.0 * mm;

    pmt_solid_maker.bridle_width = bridle_width;
    pmt_solid_maker.bridle_radius = bridle_radius;
    pmt_solid_maker.frill_width = frill_width;
    pmt_solid_maker.frill_radius = frill_radius;

    pmt_solid_maker.pmt_tpb_thickness = PMTTPBThickness;
    pmt_solid_maker.cylinder_half_height = fiducial_lar_half_height;
    pmt_solid_maker.cylinder_radius = fiducial_lar_radius;
    pmt_solid_maker.pmt_protrusion_distance = pmt_protrusion_distance;
    pmt_solid_maker.CreateBaseSolids();

    fPMTPositions.clear();
    for(std::pair<std::string const, G4ThreeVector> const &pmt_position : pmt_solid_maker.pmt_positions) {
        fPMTPositions.push_back(pmt_position.second);
    }

    // now let's build PMTs using J4SolidMaker

    fPMTVacuumWall = pmt_solid_maker.pmt_wall_vacuum_solid.get();
    fPMTVacuumCaps = pmt_solid_maker.pmt_cap_vacuum_solid.get();
    fPMTVacuumWall_log = new G4LogicalVolume(fPMTVacuumWall, G4Material::GetMaterial("Vacuum"), "PMTVacuumWallLog");
    fPMTVacuumCaps_log = new G4LogicalVolume(fPMTVacuumCaps, G4Material::GetMaterial("Vacuum"), "PMTVacuumCapsLog");

    fPMTCoatedWall = pmt_solid_maker.pmt_wall_glass_solid.get();
    fPMTCoatedCaps = pmt_solid_maker.pmt_cap_glass_solid.get();
    fPMTCoatedWall_log = new G4LogicalVolume(fPMTCoatedWall, G4Material::GetMaterial("Glass"), "PMTCoatedWallLog");
    fPMTCoatedCaps_log = new G4LogicalVolume(fPMTCoatedCaps, G4Material::GetMaterial("Glass"), "PMTCoatedCapsLog");

    fPMTUncoatedWall = pmt_solid_maker.pmt_wall_glass_solid.get();
    fPMTUncoatedCaps = pmt_solid_maker.pmt_cap_glass_solid.get();
    fPMTUncoatedWall_log = new G4LogicalVolume(fPMTUncoatedWall, G4Material::GetMaterial("Glass"), "PMTUncoatedWallLog");
    fPMTUncoatedCaps_log = new G4LogicalVolume(fPMTUncoatedCaps, G4Material::GetMaterial("Glass"), "PMTUncoatedCapsLog");

    // now get TPB coating
    fTPBCoatingWall = pmt_solid_maker.pmt_wall_tpb_solid.get();
    fTPBCoatingCaps = pmt_solid_maker.pmt_cap_tpb_solid.get();
    fTPBCoatingWall_log = new G4LogicalVolume(fTPBCoatingWall, G4Material::GetMaterial("TPBPMT"), "TPBCoatingWallLog");
    fTPBCoatingCaps_log = new G4LogicalVolume(fTPBCoatingCaps, G4Material::GetMaterial("TPBPMT"), "TPBCoatingCapsLog");

    // let's build the shiny guy who is at C406R0
    // we are placing it in our loop over PMT position IDs!
    G4double shiny_radius = 205.0 / 2.0 * mm; // making the assumption that it is the same radius as the bridle
    G4double shiny_half_height = 0.125 * mm;
    fShinyC406R0 = new G4Tubs("ShinyC406R0", 0*cm, shiny_radius, shiny_half_height, 0*deg, 360*deg);
    fShinyC406R0_log= new G4LogicalVolume(fShinyC406R0, G4Material::GetMaterial("PTFE"), "ShinyC406R0");

    // speaking of shiny, let's build our shiny reflective circles that go on top and bottom of detector
    // documentation (https://docdb.lns.mit.edu/cgi-bin/captainmills/ShowDocument?docid=567) says it's 5.25in radius
    G4double top_shiny_radius = 13.335 * cm;
    fShinyTop = new G4Tubs("ShinyTop", 0*cm, top_shiny_radius, shiny_half_height, 0*deg, 360*deg);
    fShinyTop_log= new G4LogicalVolume(fShinyTop, G4Material::GetMaterial("PTFE"), "ShinyTop");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fiducial_lar_half_height - shiny_half_height), fShinyTop_log, "ShinyTop", fFiducialLAr_log, false, 0, true);
    fShinyBottom = new G4Tubs("ShinyBottom", 0*cm, top_shiny_radius, shiny_half_height, 0*deg, 360*deg);
    fShinyBottom_log= new G4LogicalVolume(fShinyBottom, G4Material::GetMaterial("PTFE"), "ShinyBottom");
    new G4PVPlacement(0, G4ThreeVector(0, 0, - fiducial_lar_half_height + shiny_half_height), fShinyBottom_log, "ShinyBottom", fFiducialLAr_log, false, 0, true);

    double pmt_radius_cm = (fiducial_lar_radius + (pmt_solid_maker.pmt_radius - pmt_protrusion_distance)) / cm;
    double pmt_height_cm = (fiducial_lar_half_height + (pmt_solid_maker.pmt_radius - pmt_protrusion_distance)) / cm;
    double bridle_height_cm = (fiducial_lar_half_height - (bridle_width / 2.0)) / cm;
    double frill_height_cm = (fiducial_lar_half_height - (frill_width / 2.0)) / cm;
    double bridle_radius_cm = fiducial_lar_radius / cm;
    double frill_radius_cm =  fiducial_lar_radius / cm;

    std::cout << "pmt_radius_cm = " << pmt_radius_cm << std::endl;
    std::cout << "pmt_height_cm = " << pmt_height_cm << std::endl;

    // let's start by giving the pmt glass some reflection properties
    G4OpticalSurface *UncoatedPMTGlassOpticalSurface = new G4OpticalSurface("UncoatedPMTGlassOpticalSurface");
    G4OpticalSurface *CoatedPMTGlassOpticalSurface = new G4OpticalSurface("CoatedPMTGlassOpticalSurface");

    // define uncoated pmts --> ground
    UncoatedPMTGlassOpticalSurface->SetModel(unified);
    UncoatedPMTGlassOpticalSurface->SetType(dielectric_dielectric);
    //UncoatedPMTGlassOpticalSurface->SetFinish(groundfrontpainted); // 100% Lambertian (diffuse) reflections
    UncoatedPMTGlassOpticalSurface->SetFinish(ground);
    UncoatedPMTGlassOpticalSurface->SetSigmaAlpha(0.5);

    // define coated pmts --> groundfrontpainted
    CoatedPMTGlassOpticalSurface->SetModel(unified);
    CoatedPMTGlassOpticalSurface->SetType(dielectric_dielectric);
    CoatedPMTGlassOpticalSurface->SetFinish(groundfrontpainted); // 100% Lambertian (diffuse) reflections

    // define reflectivity for PMT glass (to be used on both coated and uncoated pmts)
    std::vector<G4double> PMTGlassEnergy = {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV, 2.138*eV, 2.25*eV,  2.38*eV,
                                            2.48*eV,  2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV,  3.124*eV, 3.457*eV,
                                            3.643*eV, 3.812*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
                                            7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };

    G4double uvReflection = 0.05;
    G4double visUncoatedReflection = 0.10;

    std::vector<G4double> PMTUncoatedGlassReflection = { visUncoatedReflection, visUncoatedReflection, visUncoatedReflection, visUncoatedReflection, visUncoatedReflection,
                                                    visUncoatedReflection, visUncoatedReflection, visUncoatedReflection, visUncoatedReflection, visUncoatedReflection,
                                                    visUncoatedReflection, visUncoatedReflection, visUncoatedReflection, visUncoatedReflection,
                                                    visUncoatedReflection, visUncoatedReflection, visUncoatedReflection, visUncoatedReflection, uvReflection, uvReflection, uvReflection,
                                                    uvReflection, uvReflection, uvReflection, uvReflection};

    G4double uvTrans = 1.0 - uvReflection;
    G4double visUncoatedTrans = 1.0 - visUncoatedReflection;

    std::vector<G4double> PMTUncoatedGlassTrans = { visUncoatedTrans, visUncoatedTrans, visUncoatedTrans, visUncoatedTrans, visUncoatedTrans,
                                                    visUncoatedTrans, visUncoatedTrans, visUncoatedTrans, visUncoatedTrans, visUncoatedTrans,
                                                    visUncoatedTrans, visUncoatedTrans, visUncoatedTrans, visUncoatedTrans,
                                                    visUncoatedTrans, visUncoatedTrans, visUncoatedTrans, visUncoatedTrans, uvTrans, uvTrans, uvTrans,
                                                    uvTrans, uvTrans, uvTrans, uvTrans};

    std::vector<G4double> pp = {0.602*eV,  1.03*eV, 2.138*eV, 2.845*eV, 3.124*eV,  4.086*eV, 5.474*eV, 6.262*eV, 7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV};
    std::vector<G4double> specularlobe = {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4};
    std::vector<G4double> specularspike = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    std::vector<G4double> backscatter = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

    G4MaterialPropertiesTable* PMTUncoatedGlassMPT = new G4MaterialPropertiesTable();
    PMTUncoatedGlassMPT->AddProperty("SPECULARLOBECONSTANT", pp, specularlobe);
    PMTUncoatedGlassMPT->AddProperty("SPECULARSPIKECONSTANT", pp, specularspike);
    PMTUncoatedGlassMPT->AddProperty("BACKSCATTERCONSTANT", pp, backscatter);
    PMTUncoatedGlassMPT->AddProperty("REFLECTIVITY", PMTGlassEnergy, PMTUncoatedGlassReflection);
    PMTUncoatedGlassMPT->AddProperty("TRANSMITTANCE", PMTGlassEnergy, PMTUncoatedGlassTrans);
    UncoatedPMTGlassOpticalSurface->SetMaterialPropertiesTable(PMTUncoatedGlassMPT);

    G4MaterialPropertiesTable* PMTCoatedGlassMPT = new G4MaterialPropertiesTable();
    PMTCoatedGlassMPT->AddProperty("REFLECTIVITY", PMTGlassEnergy, PMTUncoatedGlassReflection);
    PMTCoatedGlassMPT->AddProperty("TRANSMITTANCE", PMTGlassEnergy, PMTUncoatedGlassTrans);
    CoatedPMTGlassOpticalSurface->SetMaterialPropertiesTable(PMTCoatedGlassMPT);

    // and surface properties for the frill + bridle
    // Definition of MPT for Plastic frills
    //std::vector<G4double> plastic_Energy = { 1.0*eV,1.2*eV,2.5*eV,3.0*eV,3.4*eV,6.5*eV,10.0*eV,12.6*eV };
    //std::vector<G4double> plastic_reflect = {0.10, 0.10, 0.25, 0.30, 0.10, 0.05, 0.01, 0.01};

    G4double uv_plastic_refl = 0.05;
    G4double vis_plastic_refl = 0.30;

    std::vector<G4double> plastic_reflection = { vis_plastic_refl, vis_plastic_refl, vis_plastic_refl, vis_plastic_refl, vis_plastic_refl,
                                                    vis_plastic_refl, vis_plastic_refl, vis_plastic_refl, vis_plastic_refl, vis_plastic_refl,
                                                    vis_plastic_refl, vis_plastic_refl, vis_plastic_refl, vis_plastic_refl,
                                                    vis_plastic_refl, vis_plastic_refl, vis_plastic_refl, vis_plastic_refl, uv_plastic_refl, uv_plastic_refl, uv_plastic_refl,
                                                    uv_plastic_refl, uv_plastic_refl, uv_plastic_refl, uv_plastic_refl};

    G4OpticalSurface *PlasticOpticalSurface = new G4OpticalSurface("PlasticOpticalSurface");

    PlasticOpticalSurface->SetModel(unified);
    PlasticOpticalSurface->SetType(dielectric_dielectric);
    PlasticOpticalSurface->SetFinish(groundfrontpainted); // 100% Lambertian (diffuse) reflections

    G4MaterialPropertiesTable *Plastic_MT = new G4MaterialPropertiesTable();
    //Plastic_MT->AddProperty("REFLECTIVITY", plastic_Energy, plastic_reflect);
    Plastic_MT->AddProperty("REFLECTIVITY", PMTGlassEnergy, plastic_reflection);
    PlasticOpticalSurface->SetMaterialPropertiesTable(Plastic_MT);

    // now let's make a surface for TPB
    G4double uv_tpb_refl = 0.01;
    G4double vis_tpb_refl = 0.05;

    std::vector<G4double> tpb_reflection = { vis_tpb_refl, vis_tpb_refl, vis_tpb_refl, vis_tpb_refl, vis_tpb_refl,
                                                    vis_tpb_refl, vis_tpb_refl, vis_tpb_refl, vis_tpb_refl, vis_tpb_refl,
                                                    vis_tpb_refl, vis_tpb_refl, vis_tpb_refl, vis_tpb_refl,
                                                    vis_tpb_refl, vis_tpb_refl, vis_tpb_refl, vis_tpb_refl, uv_tpb_refl, uv_tpb_refl, uv_tpb_refl,
                                                    uv_tpb_refl, uv_tpb_refl, uv_tpb_refl, uv_tpb_refl};
    G4double uv_tpb_trans = 1.0 - uv_tpb_refl;
    G4double vis_tpb_trans = 1.0 - vis_tpb_refl;

    std::vector<G4double> tpb_transmission = { vis_tpb_trans, vis_tpb_trans, vis_tpb_trans, vis_tpb_trans, vis_tpb_trans,
                                                    vis_tpb_trans, vis_tpb_trans, vis_tpb_trans, vis_tpb_trans, vis_tpb_trans,
                                                    vis_tpb_trans, vis_tpb_trans, vis_tpb_trans, vis_tpb_trans,
                                                    vis_tpb_trans, vis_tpb_trans, vis_tpb_trans, vis_tpb_trans, uv_tpb_trans, uv_tpb_trans, uv_tpb_trans,
                                                    uv_tpb_trans, uv_tpb_trans, uv_tpb_trans, uv_tpb_trans};

    G4OpticalSurface *TPBOpticalSurface = new G4OpticalSurface("TPBOpticalSurface");

    TPBOpticalSurface->SetModel(unified);
    TPBOpticalSurface->SetType(dielectric_dielectric);
    TPBOpticalSurface->SetFinish(groundfrontpainted); // 100% Lambertian (diffuse) reflections

    G4MaterialPropertiesTable *TPB_MT = new G4MaterialPropertiesTable();
    TPB_MT->AddProperty("REFLECTIVITY", PMTGlassEnergy, tpb_reflection);
    TPB_MT->AddProperty("TRANSMITTANCE", PMTGlassEnergy, tpb_transmission);
    TPBOpticalSurface->SetMaterialPropertiesTable(TPB_MT);

    // put tpb surface on our foils
    new G4LogicalBorderSurface("TPBFoilSides_Surface", fFiducialLAr_phys, fTPBFoilSides_phys, TPBOpticalSurface);
    new G4LogicalBorderSurface("TPBFoilTop_Surface", fFiducialLAr_phys, fTPBFoilTop_phys, TPBOpticalSurface);
    new G4LogicalBorderSurface("TPBFoilBottom_Surface", fFiducialLAr_phys, fTPBFoilBottom_phys, TPBOpticalSurface);

    // Place the frills
    for(std::pair<std::string const, std::shared_ptr<G4VSolid>> & p : pmt_solid_maker.frill_solids) {
        std::string name = p.first;
        std::shared_ptr<G4VSolid> frill_solid = p.second;
        G4ThreeVector & frill_position = pmt_solid_maker.frill_positions[name];
        std::shared_ptr<G4RotationMatrix> & frill_rotation = pmt_solid_maker.frill_rotations[name];
        G4LogicalVolume * frill_log;
        std::map<std::shared_ptr<G4VSolid>, std::tuple<unsigned int, std::shared_ptr<G4LogicalVolume>>>::iterator it = fFrillLogicalVolumes.find(frill_solid);
        unsigned int copy_number_k;
        if(it == fFrillLogicalVolumes.end()) {
            copy_number_k = 0;
            frill_log = new G4LogicalVolume(frill_solid.get(), G4Material::GetMaterial("Plastic"), name);
            fFrillLogicalVolumes[frill_solid] = {1, std::shared_ptr<G4LogicalVolume>(frill_log)};

            G4LogicalSkinSurface * frill_skin = new G4LogicalSkinSurface(name + "_Surface", frill_log, PlasticOpticalSurface);
            fLogicalSkinSurfaces.insert({{false, frill_solid}, std::shared_ptr<G4LogicalSkinSurface>(frill_skin)});
        } else {
            frill_log = std::get<1>(it->second).get();
            std::get<0>(it->second) += 1;
            copy_number_k = std::get<0>(it->second);
        }
        fPlacements.emplace_back(new G4PVPlacement(frill_rotation.get(), frill_position, frill_log, name, fFiducialLAr_log, false, copy_number_k, true));
    }

    // Place the bridles
    for(std::pair<std::string const, std::shared_ptr<G4VSolid>> & p : pmt_solid_maker.bridle_solids) {
        std::string name = p.first;
        std::shared_ptr<G4VSolid> bridle_solid = p.second;
        G4ThreeVector & bridle_position = pmt_solid_maker.bridle_positions[name];
        std::shared_ptr<G4RotationMatrix> & bridle_rotation = pmt_solid_maker.bridle_rotations[name];
        G4LogicalVolume * bridle_log;
        std::map<std::shared_ptr<G4VSolid>, std::tuple<unsigned int, std::shared_ptr<G4LogicalVolume>>>::iterator it = fBridleLogicalVolumes.find(bridle_solid);
        unsigned int copy_number_k;
        if(it == fBridleLogicalVolumes.end()) {
            copy_number_k = 0;
            bridle_log = new G4LogicalVolume(bridle_solid.get(), G4Material::GetMaterial("BlackPlastic"), "Bridle_" + name);
            fBridleLogicalVolumes[bridle_solid] = {1, std::shared_ptr<G4LogicalVolume>(bridle_log)};

            G4LogicalSkinSurface * bridle_skin = new G4LogicalSkinSurface("Bridle_" + name + "_Surface", bridle_log, PlasticOpticalSurface);
            fLogicalSkinSurfaces.insert({{false, bridle_solid}, std::shared_ptr<G4LogicalSkinSurface>(bridle_skin)});
        } else {
            bridle_log = std::get<1>(it->second).get();
            std::get<0>(it->second) += 1;
            copy_number_k = std::get<0>(it->second);
        }
        fPlacements.emplace_back(new G4PVPlacement(bridle_rotation.get(), bridle_position, bridle_log, name, fFiducialLAr_log, false, copy_number_k, true));
    }

    // Place the pmts
    for(std::pair<std::string const, std::shared_ptr<G4VSolid>> & p : pmt_solid_maker.pmt_solids) {
        std::string name = p.first;
        std::shared_ptr<G4VSolid> pmt_solid = p.second;
        G4ThreeVector pmt_position = pmt_solid_maker.pmt_positions[name];
        std::shared_ptr<G4RotationMatrix> pmt_rotation = pmt_solid_maker.pmt_rotations[name];
        G4LogicalVolume * pmt_log;

        std::shared_ptr<G4VSolid> tpb_solid = (pmt_solid_maker.tpb_solids.count(name) > 0) ? pmt_solid_maker.tpb_solids[name] : nullptr;
        G4ThreeVector tpb_position = (tpb_solid != nullptr) ? pmt_solid_maker.tpb_positions[name] : G4ThreeVector(0, 0, 0);
        std::shared_ptr<G4RotationMatrix> tpb_rotation = (tpb_solid != nullptr) ? pmt_solid_maker.tpb_rotations[name] : nullptr;
        G4LogicalVolume * tpb_log;

        bool pmt_on_cap, coated;
        int pmt_row, pmt_col, pmt_ring, pmt_number;
        pmt_solid_maker.ParsePMTID(name, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);
        std::string pmt_name = (coated) ? "Coated" : "Uncoated";
        pmt_name += (pmt_on_cap) ? "Cap" : "Wall";

        G4OpticalSurface * pmt_optical_surface;
        unsigned int copy_number_k;
        std::map<std::shared_ptr<G4VSolid>, std::tuple<unsigned int, std::shared_ptr<G4LogicalVolume>>>::iterator tpb_it;

        if(tpb_solid != nullptr) {
            tpb_it = fTPBLogicalVolumes.find(tpb_solid);
            if(tpb_it == fTPBLogicalVolumes.end()) {
                copy_number_k = 0;
                if(pmt_on_cap)
                    tpb_log = fTPBCoatingCaps_log;
                else
                    tpb_log = fTPBCoatingWall_log;
                fTPBLogicalVolumes[tpb_solid] = {1, std::shared_ptr<G4LogicalVolume>(tpb_log)};

                std::map<std::tuple<bool, std::shared_ptr<G4VSolid>>, std::tuple<unsigned int, std::shared_ptr<G4LogicalVolume>>>::iterator pmt_it;
                pmt_it = fPMTLogicalVolumes.find({coated, pmt_solid});
                if(pmt_it == fPMTLogicalVolumes.end()) {
                    copy_number_k = 0;
                    pmt_optical_surface = CoatedPMTGlassOpticalSurface;
                    if(pmt_on_cap)
                        pmt_log = fPMTCoatedCaps_log;
                    else
                        pmt_log = fPMTCoatedWall_log;
                    fPMTLogicalVolumes[{coated, pmt_solid}] = {1, std::shared_ptr<G4LogicalVolume>(pmt_log)};

                    G4LogicalSkinSurface * pmt_skin = new G4LogicalSkinSurface("PMTGlass_" + pmt_name + "_Surface", pmt_log, pmt_optical_surface);
                    fLogicalSkinSurfaces.insert({{coated, pmt_solid}, std::shared_ptr<G4LogicalSkinSurface>(pmt_skin)});

                    std::map<std::shared_ptr<G4VSolid>, std::tuple<unsigned int, std::shared_ptr<G4LogicalVolume>>>::iterator vac_it;
                    std::shared_ptr<G4VSolid> vacuum_solid = pmt_solid_maker.vacuum_solids[name];
                    G4LogicalVolume * vacuum_log;

                    vac_it = fVacuumLogicalVolumes.find(vacuum_solid);
                    if(vac_it == fVacuumLogicalVolumes.end()) {
                        copy_number_k = 0;
                        vacuum_log = new G4LogicalVolume(vacuum_solid.get(), G4Material::GetMaterial("Vacuum"), std::string("PMTVacuum_") + ((pmt_on_cap) ? "Cap" : "Wall") + "Log");
                        fVacuumLogicalVolumes[vacuum_solid] = {1, std::shared_ptr<G4LogicalVolume>(vacuum_log)};
                    } else {
                        vacuum_log = std::get<1>(vac_it->second).get();
                        std::get<0>(vac_it->second) += 1;
                        copy_number_k = std::get<0>(vac_it->second);
                    }
                    fPlacements.emplace_back(new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), vacuum_log, "PMTVacuum_" + name, pmt_log, false, copy_number_k, true));
                } else {
                    pmt_log = std::get<1>(pmt_it->second).get();
                    std::get<0>(pmt_it->second) += 1;
                    copy_number_k = std::get<0>(pmt_it->second);
                }
                fPlacements.emplace_back(new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), pmt_log, "PMTGlass_" + name, tpb_log, false, copy_number_k, true));
            } else {
                tpb_log = std::get<1>(tpb_it->second).get();
                std::get<0>(tpb_it->second) += 1;
                copy_number_k = std::get<0>(tpb_it->second);
            }
            fPlacements.emplace_back(new G4PVPlacement(tpb_rotation.get(), tpb_position, tpb_log, "PMTTPB_" + name, fFiducialLAr_log, false, copy_number_k, true));
            fLogicalBorderSurfaces.emplace_back(new G4LogicalBorderSurface("PMTTPB_" + name + "_Surface", fFiducialLAr_phys, fPlacements.back().get(), TPBOpticalSurface));
        } else {
            std::map<std::tuple<bool, std::shared_ptr<G4VSolid>>, std::tuple<unsigned int, std::shared_ptr<G4LogicalVolume>>>::iterator pmt_it;

            pmt_it = fPMTLogicalVolumes.find({coated, pmt_solid});
            if(pmt_it == fPMTLogicalVolumes.end()) {
                copy_number_k = 0;
                pmt_optical_surface = UncoatedPMTGlassOpticalSurface;
                if(pmt_on_cap)
                    pmt_log = fPMTUncoatedCaps_log;
                else
                    pmt_log = fPMTUncoatedWall_log;
                fPMTLogicalVolumes[{coated, pmt_solid}] = {1, std::shared_ptr<G4LogicalVolume>(pmt_log)};

                G4LogicalSkinSurface * pmt_skin = new G4LogicalSkinSurface("PMTGlass_" + pmt_name + "_Surface", pmt_log, pmt_optical_surface);
                fLogicalSkinSurfaces.insert({{coated, pmt_solid}, std::shared_ptr<G4LogicalSkinSurface>(pmt_skin)});

                std::map<std::shared_ptr<G4VSolid>, std::tuple<unsigned int, std::shared_ptr<G4LogicalVolume>>>::iterator vac_it;
                std::shared_ptr<G4VSolid> vacuum_solid = pmt_solid_maker.vacuum_solids[name];
                G4LogicalVolume * vacuum_log;

                vac_it = fVacuumLogicalVolumes.find(vacuum_solid);
                if(vac_it == fVacuumLogicalVolumes.end()) {
                    copy_number_k = 0;
                    vacuum_log = new G4LogicalVolume(vacuum_solid.get(), G4Material::GetMaterial("Vacuum"), std::string("PMTVacuum_") + ((pmt_on_cap) ? "Cap" : "Wall") + "Log");
                    fVacuumLogicalVolumes[vacuum_solid] = {1, std::shared_ptr<G4LogicalVolume>(vacuum_log)};
                } else {
                    vacuum_log = std::get<1>(vac_it->second).get();
                    std::get<0>(vac_it->second) += 1;
                    copy_number_k = std::get<0>(vac_it->second);
                }
                fPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0,0,0), vacuum_log, "PMTVacuum_" + name, pmt_log, false, copy_number_k, true));
            } else {
                pmt_log = std::get<1>(pmt_it->second).get();
                std::get<0>(pmt_it->second) += 1;
                copy_number_k = std::get<0>(pmt_it->second);
            }
            fPlacements.emplace_back(new G4PVPlacement(pmt_rotation.get(), pmt_position, pmt_log, "PMTGlass_" + name, fFiducialLAr_log, false, copy_number_k, true));
        }
    }

    std::string shiny_position_id = "C406R0";
    G4ThreeVector shiny_position = pmt_solid_maker.GetPMTPosition(shiny_position_id);
    shiny_position.setZ(fiducial_lar_half_height - shiny_half_height);
    // now place our shiny circle!
    fPlacements.emplace_back(new G4PVPlacement(0, shiny_position, fShinyC406R0_log, "ShinyC406R0", fFiducialLAr_log, false, 0, true));

    if (SourceRodIn){
        // now let's make our source rod
        // i'm just doing the bayonet for now since it is >1m tall, it should only matter for sodium runs above -50cm
        G4double rod_inner_radius = 0.0*cm;
        G4double rod_outer_radius = 6.31 * mm;
        G4double rod_height = fiducial_lar_half_height - SourceRodLocation - (shiny_half_height * 2.0);
        fSourceRod = new G4Tubs("SourceRod", rod_inner_radius, rod_outer_radius, rod_height/2, 0, 360*deg);
        G4ThreeVector rodPosition(0.0*cm, 0.0*cm, SourceRodLocation + rod_height/2);

        // let's make a source rod optical surface
        G4OpticalSurface *SourceRodOpticalSurface = new G4OpticalSurface("SourceRodOpticalSurface");

        SourceRodOpticalSurface->SetModel(unified);
        SourceRodOpticalSurface->SetType(dielectric_metal);
        SourceRodOpticalSurface->SetFinish(polished);

        // define reflectivity for stainless steel
        std::vector<G4double> StainlessSteelEnergy = {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV, 2.138*eV, 2.25*eV,  2.38*eV,
                                                2.48*eV,  2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV,  3.124*eV, 3.457*eV,
                                                3.643*eV, 3.812*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
                                                7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };

        G4double uvReflection = 0.30;
        G4double visReflection = 0.60;

        std::vector<G4double> StainlessSteelReflection = { visReflection, visReflection, visReflection, visReflection, visReflection,
                                                        visReflection, visReflection, visReflection, visReflection, visReflection,
                                                        visReflection, visReflection, visReflection, visReflection,
                                                        visReflection, visReflection, visReflection, visReflection, uvReflection, uvReflection, uvReflection,
                                                        uvReflection, uvReflection, uvReflection, uvReflection};

        G4MaterialPropertiesTable* StainlessSteelMPT = new G4MaterialPropertiesTable();
        StainlessSteelMPT->AddProperty("REFLECTIVITY", StainlessSteelEnergy, StainlessSteelReflection);
        SourceRodOpticalSurface->SetMaterialPropertiesTable(StainlessSteelMPT);

        if (CobaltSourceRun or SodiumSourceRun){
            // now make source pellet -- really guessing on these measurements
            G4double pellet_radius = 4.0*mm;
            G4double pellet_height = 3.0*mm;
            fSourcePellet = new G4Tubs("SourcePellet", 0, pellet_radius, pellet_height/2.0, 0, 360*deg);
            G4NistManager* nistManager = G4NistManager::Instance();
            G4Material* fSource = nistManager->FindOrBuildMaterial("G4_Na");
            if (CobaltSourceRun){
                fSource = nistManager->FindOrBuildMaterial("G4_Co");
            }
            fSourcePellet_log  = new G4LogicalVolume(fSourcePellet, fSource, "SourcePelletLog");

            if (TrainingSource) {
                // this is the special case where we are injecting sodium events randomly within the detector to get training data
                // so we have special logic for placing the source pellet and making steel housing
                G4double source_pellet_housing_height = 1.0 * cm;
                fSourcePelletHousing = new G4Tubs("SourcePelletHousing", rod_inner_radius, rod_outer_radius, source_pellet_housing_height/2, 0, 360*deg);

                // Subtract the pellet from the housing
                G4SubtractionSolid* housingWithPelletHole = new G4SubtractionSolid("HousingWithPelletHole", fSourcePelletHousing, fSourcePellet, nullptr, G4ThreeVector(0, 0, 0));
                fSourcePelletHousing_log = new G4LogicalVolume(housingWithPelletHole, G4Material::GetMaterial("Steel"), "fSourcePelletHousingLog");

                G4ThreeVector pelletPosition(DecayX, DecayY, DecayZ);
                std::cout << "placing sodium decay at " << DecayX << ", " << DecayY << ", " << DecayZ << std::endl;

                // First, subtract the pellet from the rod
                G4SubtractionSolid* rodWithPelletHole = new G4SubtractionSolid("RodWithPelletHole", fSourceRod, fSourcePellet, nullptr, pelletPosition);

                // Next, subtract the housing (with pellet hole) from the rod
                G4SubtractionSolid* rodWithPelletAndHousingHole = new G4SubtractionSolid("RodWithPelletAndHousingHole", rodWithPelletHole, housingWithPelletHole, nullptr, pelletPosition);

                // Create a logical volume for the rod with both subtractions applied
                fSourceRod_log = new G4LogicalVolume(rodWithPelletAndHousingHole, G4Material::GetMaterial("Steel"), "fSourceRodLogWithHoles");

                // Place the modified rod (with pellet and housing subtracted) in the fiducial argon volume
                new G4PVPlacement(nullptr, rodPosition, fSourceRod_log, "RodWithPelletAndHousingHole", fFiducialLAr_log, false, 0, true);

                // Place the source pellet and housing within the detector
                new G4PVPlacement(nullptr, pelletPosition, fSourcePellet_log, "SourcePellet", fFiducialLAr_log, false, 0, true);
                new G4PVPlacement(nullptr, pelletPosition, fSourcePelletHousing_log, "SourcePelletHousing", fFiducialLAr_log, false, 0, true);

            } else {
                // standard case where the source pellet is inserted 1/4 cm into the end of the source rod
                G4double inset = 1.0 * mm;
                G4ThreeVector pelletPosition(0.0*cm, 0.0*cm, SourceRodLocation + pellet_height/2.0 + inset);

                // again let's first subtract the pellet from the rod, make logical vol for rod, and place rod
                G4SubtractionSolid* rodWithPelletHole = new G4SubtractionSolid("RodWithPelletHole", fSourceRod, fSourcePellet, nullptr, pelletPosition);
                fSourceRod_log = new G4LogicalVolume(rodWithPelletHole, G4Material::GetMaterial("Steel"), "fSourceRodLogWithHole");
                new G4PVPlacement(nullptr, rodPosition, fSourceRod_log, "RodWithPelletHole", fFiducialLAr_log, false, 0, true);

                // now place sodium pellet
                new G4PVPlacement(nullptr, pelletPosition, fSourcePellet_log, "SourcePellet", fSourceRod_log, false, 0, true);
            }
        } else {
            // this is the case where we have the source rod in the detector but no pellet -- so let's place the rod like normal
            fSourceRod_log = new G4LogicalVolume(fSourceRod,  G4Material::GetMaterial("Steel"), "fSourceRodLog");
            new G4PVPlacement(nullptr, rodPosition, fSourceRod_log, "SourceRod", fFiducialLAr_log, false, 0, true);
        }

        // now add our optical surface to our source rod
        new G4LogicalSkinSurface("SourceRod_Surface", fSourceRod_log, SourceRodOpticalSurface);
    }

    VisAttributes(SourceRodIn);
    SurfaceProperties();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMMainVolume::VisAttributes(G4bool SourceRodIn)
{
    // let's make other layers of detector invisible
    auto lilac = new G4VisAttributes(G4Colour(153/255., 153/255., 153/255.));
    //lilac->SetForceSolid(true);
    auto salmon = new G4VisAttributes(G4Colour(255/255., 128/255., 128/255.));
    //salmon->SetForceSolid(true);
    auto dark_blue = new G4VisAttributes(G4Colour(0., 0., 128/255.));
    dark_blue->SetForceSolid(true);
    auto green = new G4VisAttributes(G4Colour(51/255., 153/255., 102/255.));
    //green->SetForceSolid(true);
    auto orange = new G4VisAttributes(G4Colour(255/255., 153/255., 0.));
    //orange->SetForceSolid(true);
    auto teal = new G4VisAttributes(G4Colour(0., 128/255., 128/255.));
    //teal->SetForceSolid(true);
    auto pink = new G4VisAttributes(G4Colour(255/255., 0., 255/255.));
    //pink->SetForceSolid(true);

    fCryoVessel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fVacuum_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fInnerJacket_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fOuterLAr_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fInnerFrame_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fTPBFoil_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fFiducialLAr_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    //fCryoVessel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fVacuum_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fInnerJacket_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fOuterLAr_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fInnerFrame_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fTPBFoil_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTCoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTUncoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fTPBCoating_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPhotocathCoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPhotocathUncoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto pmt_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8)); //grey idk
    pmt_va->SetForceSolid(true);
    //fPMTCoatedWall_log->SetVisAttributes(pmt_va);
    //fPMTCoatedCaps_log->SetVisAttributes(pmt_va);
    //fPMTUncoatedWall_log->SetVisAttributes(pmt_va);
    //fPMTUncoatedCaps_log->SetVisAttributes(pmt_va);
    pink->SetForceSolid(true);
    fPMTCoatedWall_log->SetVisAttributes(pink);
    fPMTCoatedCaps_log->SetVisAttributes(pink);
    fPMTUncoatedWall_log->SetVisAttributes(teal);
    fPMTUncoatedCaps_log->SetVisAttributes(teal);

    auto tpb_coating_va = new G4VisAttributes(G4Colour(0., 1., 0.)); //green
    tpb_coating_va->SetForceSolid(true);
    //fTPBCoatingWall_log->SetVisAttributes(tpb_coating_va);
    //fTPBCoatingCaps_log->SetVisAttributes(tpb_coating_va);
    fTPBCoatingWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fTPBCoatingCaps_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    salmon->SetForceSolid(true);
    teal->SetForceSolid(true);
    for(G4LogicalVolume * log : GetFrillLogicalVolumes()) {
        log->SetVisAttributes(dark_blue);
    }
    for(G4LogicalVolume * log : GetBridleLogicalVolumes()) {
        log->SetVisAttributes(salmon);
    }
    //fFrillWall_log->SetVisAttributes(teal);
    //fFrillCaps_log->SetVisAttributes(teal);
    //fBridleWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fBridleCaps_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fFrillWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fFrillCaps_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    // add shiny guy!
    auto red = new G4VisAttributes(G4Colour(1., 0., 0.));
    red->SetForceSolid(true);
    fShinyC406R0_log->SetVisAttributes(red);
    fShinyTop_log->SetVisAttributes(red);
    fShinyBottom_log->SetVisAttributes(red);

    // tpb foils!
    fTPBFoilSides_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fTPBFoilTop_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fTPBFoilBottom_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    // make our source rod green
    if (SourceRodIn)
        fSourceRod_log->SetVisAttributes(tpb_coating_va);
    //fSourceRod_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    //// and make sodium pellet blue
    //auto sodium_pellet_va = new G4VisAttributes(G4Colour(0., 0., 1.)); // blue
    //sodium_pellet_va->SetForceSolid(true);
    //fSodiumSourcePellet_log->SetVisAttributes(sodium_pellet_va);
    //fSourcePellet_log->SetVisAttributes(G4VisAttributes::GetInvisible());


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMMainVolume::SurfaceProperties()
{
    // now it's time to define optical surface properties
    G4OpticalSurface *ReflectorOpticalSurface = new G4OpticalSurface("ReflectorOpticalSurface");

    ReflectorOpticalSurface->SetModel(unified);
    ReflectorOpticalSurface->SetType(dielectric_metal);
    ReflectorOpticalSurface->SetFinish(polished);

    // let's add reflection properties as a function of energy
    std::vector<G4double> MylarReflectionEnergy = {1.239835823924519*eV, 2.1758440843642446*eV, 2.2155183016259334*eV, 2.2540941739092317*eV, 2.2965084198409396*eV,
            2.3398220945804553*eV, 2.3853127919447585*eV, 2.432339257160777*eV, 2.4816594267685304*eV, 2.5316156616832233*eV, 2.584054978035275*eV, 2.639307340578113*eV,
            2.6963833000212216*eV, 2.7566402380421113*eV, 2.820482241642379*eV, 2.8840814650497304*eV, 2.9538999741438468*eV, 3.0256914075766623*eV, 3.104173784836434*eV,
            3.1423952308218324*eV, 3.1815179074033235*eV, 3.2230571642902683*eV, 3.2634554851366535*eV, 3.306493319706959*eV, 3.352069231820206*eV, 3.445795378257322*eV,
            3.541505417163459*eV, 3.647364059603557*eV, 3.7597652297046125*eV, 3.876624049427723*eV, 4.000980170596548*eV, 4.135062961692673*eV, 4.276792392741077*eV,
            4.430351962448966*eV, 8.265572159496793*eV, 12.39835823924519*eV};

    std::vector<G4double> MylarReflection = {99.3226740868548, 99.3226740868548, 99.7628366612554, 99.7214862392674, 98.91393146613707, 98.91530448108364, 98.96965285731727, 98.76966534546983,
            98.1497114120479, 98.72409290812726, 98.32154569778484, 98.03407282162186, 97.82027684556721, 98.53168216346992, 97.48482826997746, 96.10633382322926, 96.51396948976577,
            96.65116783853793, 95.91151864322745, 95.30515359961677, 91.25574897695003, 81.90244998989932, 41.8931179140657, 24.101847317428096, 12.26957594631375, 12.06518399527252,
            12.58223308670786, 12.888602883272128, 14.106365509586524, 14.817742868938694, 15.529120228290893, 14.721537496610026, 13.103799846600765, 12.245584179933772, 0.0, 0.0};

    for (size_t ref_it = 0; ref_it < MylarReflection.size(); ref_it++){
        MylarReflection.at(ref_it) /= 100;
    }

    G4MaterialPropertiesTable *ReflectiveFoilMPT = new G4MaterialPropertiesTable();
    ReflectiveFoilMPT->AddProperty("REFLECTIVITY", MylarReflectionEnergy, MylarReflection);
    ReflectorOpticalSurface->SetMaterialPropertiesTable(ReflectiveFoilMPT);

    new G4LogicalSkinSurface("PTFE_Surface", fReflectorFoil_log, ReflectorOpticalSurface);
    new G4LogicalSkinSurface("ShinyC406R0_Surface", fShinyC406R0_log, ReflectorOpticalSurface);
    new G4LogicalSkinSurface("ShinyTop_Surface", fShinyTop_log, ReflectorOpticalSurface);
    new G4LogicalSkinSurface("ShinyBottom_Surface", fShinyBottom_log, ReflectorOpticalSurface);
}

