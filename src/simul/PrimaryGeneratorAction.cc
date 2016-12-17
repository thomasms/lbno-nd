//____________________________________________________________________________
/*!

\class    PrimaryGeneratorAction

\author   Yordan Karadzhov <Yordan.Karadzhov \at cern.ch>
	  University of Geneva

\author	  Tom Stainer <tstainer \at liv.ac.uk>
          University of Liverpool

\created  Sep 2012
\last update June 2014

*/
//____________________________________________________________________________
#include "PrimaryGeneratorAction.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction(PhysicsList* p)
: verboseMode(true), verbose_(0), phys(p), beamCount(0),primVertex_(NULL){

   particleGun = new G4GeneralParticleSource();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {

  if(particleGun) delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  //if(primVertex_) delete primVertex_; // -- probably deleted by geant, gives seg fault if delete it here.
  primVertex_ = NULL;

  if(GunOn_ == 1) {

 	particleGun->GeneratePrimaryVertex(anEvent);
	energy = particleGun->GetParticleEnergy();
  }

  else{

  //set the initial energy
  energy = event_->getNuEnergy();

  TLorentzVector vertexPos = event_->getPosition();
  primVertex_ = new G4PrimaryVertex(vertexPos.X(),
                                    vertexPos.Y(),
                                    vertexPos.Z(), 0.);
  
  //fspl
  ParticleDescrShortRecord fspl = event_->getFspl();
  AddParticleToVertex(primVertex_, &fspl);

  //add all fss
  std::vector<ParticleDescrShortRecord> fss = event_->getFssVector();
  std::vector<ParticleDescrShortRecord>::iterator itr;

  if(verbose_>0){
    std::cout << "\n\n================= Event " << beamCount << " Summary =================";
    std::cout << "\n---------------------------"
	    << "\nFinal State Primary Lepton "
	    << "\n---------------------------\n";

    fspl.printToStream(std::cout);

    std::cout << "\n--------------------------------"
	    << "\nNumber of fss particles is: "
	    << fss.size()
	    << "\n--------------------------------\n";
  }
  
  for(itr = fss.begin(); itr < fss.end(); ++itr) {
	if((*itr).getPDG() < 100000){
		AddParticleToVertex(primVertex_, &(*itr));
		if(verbose_>0)(*itr).printToStream(std::cout);
	}
  }

  anEvent->AddPrimaryVertex(primVertex_);

  //event_->printToStream(std::cout);
  //vertex->Print();
  }
 
  beamCount++;
 
}

void PrimaryGeneratorAction::AddParticleToVertex(G4PrimaryVertex* v,
                                                 ParticleDescrShortRecord* p) {
  G4PrimaryParticle* primaryParticle;
  int pdg = p->getPDG();
  TLorentzVector p4 = p->getP4();
  primaryParticle = new G4PrimaryParticle(pdg);
  primaryParticle->Set4Momentum(p4.Px(), p4.Py(), p4.Pz(), p4.E());

  /* else {
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->GetIon(6,12,0,0);
    primaryParticle = new G4PrimaryParticle(particle);
    primaryParticle->SetMomentum(p->Px(), p->Py(), p->Pz());
    primaryParticle->SetMass(particle->GetPDGMass());
  }

    primaryParticle->SetCharge(p->Charge()*(eplus/3));
    TVector3 polz;
    p->GetPolarization(polz);
    primaryParticle->SetPolarization(polz.X(), polz.Y(), polz.Z());*/
    v->SetPrimary(primaryParticle);

    //delete primaryParticle;
}

void PrimaryGeneratorAction::printProgress(double percent, int i, int n, int m) {
  static int lastPrintedPercent = -1;
  int p = (int)(percent);
  if (p > lastPrintedPercent) {
    int width = 1 + (int)log10((double)n);
    std::cout << "Completion: "       << std::setw(3)     << p << "%;    ";
    std::cout << "Events read: "      << std::setw(width) << i << "/" << n << ";    ";
    std::cout << "Events simulated: " << std::setw(width) << m << ";" << std::endl;
    lastPrintedPercent = p;
  }
}
