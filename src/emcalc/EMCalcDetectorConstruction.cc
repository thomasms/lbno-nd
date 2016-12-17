//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

//
// $Id: DetectorConstruction.cc,v 1.10 2009-11-14 18:04:20 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EMCalcDetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EMCalcDetectorConstruction::EMCalcDetectorConstruction()
:pBox(0), aMaterial(0) {
  BoxSize = 1*mm;
  DefineMaterials();
  SetMaterial("Iron");  
  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EMCalcDetectorConstruction::~EMCalcDetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* EMCalcDetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EMCalcDetectorConstruction::DefineMaterials()
{
  //
  // define Elements
  //

  G4NistManager* man = G4NistManager::Instance();

  concrete_ = man->FindOrBuildMaterial("G4_CONCRETE");
  iron_  = man->FindOrBuildMaterial("G4_Fe");
  lead_  = man->FindOrBuildMaterial("G4_Pb");
  air_   = man->FindOrBuildMaterial("G4_AIR");
  scint_ = man->FindOrBuildMaterial("G4_POLYSTYRENE");

  double density_iron, density_scint, density_MIND;
  density_iron = iron_->GetDensity();
  density_scint = scint_->GetDensity();
  density_MIND = (3.*density_iron + 2.*density_scint)/5.;

  MIND_ = new G4Material("MIND_mixture", density_MIND, 2);

  double fractionmass_iron = (3.*density_iron)/(5.*density_MIND);
  MIND_->AddMaterial(iron_, fractionmass_iron);

  double fractionmass_scint = (2.*density_scint)/(5.*density_MIND);
  MIND_->AddMaterial(scint_, fractionmass_scint);
/*
  G4double z,a;

  G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"   ,"C" , z= 6., a=  12.01*g/mole);  
  G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);
  G4Element* Ge = new G4Element("Germanium","Ge", z=32., a=  72.59*g/mole);
  G4Element* Pb = new G4Element("Lead"     ,"Pb", z=82., a= 207.19*g/mole);  
  G4Element* Bi = new G4Element("Bismuth"  ,"Bi", z=83., a= 208.98*g/mole);
  
  //
  // define materials
  //
  G4double density;
  G4int ncomponents, natoms;
  G4double fractionmass;  

  G4Material* Air = 
  new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=70.*perCent);
  Air->AddElement(O, fractionmass=30.*perCent);

  G4Material* H2l = 
  new G4Material("H2liquid", density= 70.8*mg/cm3, ncomponents=1);
  H2l->AddElement(H, fractionmass=1.);

  G4Material* H2O = 
  new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  ///H2O->SetChemicalFormula("H_2O");
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  G4Material* steam = 
  new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1);
  steam->AddMaterial(H2O, fractionmass=1.);
  steam->GetIonisation()->SetMeanExcitationEnergy(71.6*eV);  

  G4Material* BGO = 
  new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
  BGO->AddElement(O , natoms=12);
  BGO->AddElement(Ge, natoms= 3);
  BGO->AddElement(Bi, natoms= 4);  

  new G4Material("Aluminium"  , z=13., a= 26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"    , z=14., a= 28.09*g/mole, density= 2.330*g/cm3);
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  new G4Material("Iron"       , z=26., a= 55.85*g/mole, density= 7.870*g/cm3);
  new G4Material("Copper"     , z=29., a= 63.55*g/mole, density= 8.960*g/cm3);
  new G4Material("Germanium"  , z=32., a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Silver"     , z=47., a=107.87*g/mole, density= 10.50*g/cm3);
  new G4Material("Tungsten"   , z=74., a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Lead"       , z=82., a=207.19*g/mole, density= 11.35*g/cm3);
  new G4Material("Uranium"    , z=92., a=238.03*g/mole, density= 18.95*g/cm3);

  G4Material* ams = 
  new G4Material("ams", density= 7.409*g/cm3, ncomponents=3);
  ams->AddElement(Pb, fractionmass = 94.81*perCent);
  ams->AddElement(C , fractionmass =  4.79*perCent);
  ams->AddElement(H , fractionmass =  0.40*perCent);
*/
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


G4VPhysicalVolume* EMCalcDetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4Box* 
    sBox = new G4Box("Container",                       //its name
                     BoxSize/2,BoxSize/2,BoxSize/2);    //its dimensions

  G4LogicalVolume*
    lBox = new G4LogicalVolume(sBox,                    //its shape
                              aMaterial,                //its material
                              aMaterial->GetName());    //its name

  pBox = new G4PVPlacement(0,                           //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           lBox,                        //its logical volume
                           aMaterial->GetName(),        //its name
                           0,                           //its mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  //always return the root volume
  //
  return pBox;
}


void EMCalcDetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(BoxSize,"Length")
         << " of " << aMaterial->GetName() << G4endl;
}


#include "G4RunManager.hh"

void EMCalcDetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name, or build it from nist data base
  G4Material* pttoMaterial = 
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    aMaterial = pttoMaterial;
    if (pBox) G4RunManager::GetRunManager()
                             ->DefineWorldVolume(ConstructVolumes());
  } else {
    G4cout << "\n--> warning from EMCalcDetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;  
  }  
}
