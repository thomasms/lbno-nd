//____________________________________________________________________________
/*!

\class    lbnoReconstruct

\brief    Program to perform reconstruction techniques

\author	  Tom Stainer <tstainer \at liv.ac.uk>
          University of Liverpool

\created  June 2013
\last update November 2015

*/
//____________________________________________________________________________
#include <iostream>
#include <map>

#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TGeoVolume.h>

#include "LbnoPlotter.hh"
#include "ProgressBar.hh"
#include "TpcMomSmear.hh"
#include "Tracker.hh"

class myReconstructor : public LbnoPlotter
{
public:
    myReconstructor() {useMagFieldFromCard = true;}
    virtual ~myReconstructor() {};
    void reconstruct();
  
  void initHistos();
  void write();
  void initDataCards();
  void loadDataCards();

private:
  TFile file_;
  
  int hitCountCondition_;
  
  int badFileCount_;
  int pdg_,leptonPdg_;
  bool justTPC_,ccOnly_;
  double maxEnergy_,testMagField_;
  
  std::vector< std::vector<ReconData> > badFileVector_;
    
  TH1F muon_track_length_h1d;
  TH1F muon_track_length_above1GeV_h1d;
  TH1F muon_track_length_above2GeV_h1d;
  TH1F muon_track_length_above3GeV_h1d;
  TH1F muon_track_length_above4GeV_h1d;
  TH1F muon_track_length_above5GeV_h1d;
  TH2F muon_track_length_vs_mom_h2d;
  
  TH1F lower_bound_cc_muon_trans_mom_spectrum_h1d;
  TH1F upper_bound_cc_muon_trans_mom_spectrum_h1d;
  TH1F truth_cc_muon_trans_spectrum_h1d;

};

void myReconstructor::initDataCards() {
  cards_ = DataCards::getInstance();
  cards_->AddDataCardString("dataFileName", "nuData.root");
  cards_->AddDataCardString("plotFileName", "nuHistos.root");
  cards_->AddDataCardInt("numberOfRuns", 100);
  cards_->AddDataCardInt("numberOfEvents", -1);
  cards_->AddDataCardInt("plotterVerbose", 0);
  cards_->AddDataCardInt("pdg", 0);
  cards_->AddDataCardInt("leptonFlavour", 13);
  cards_->AddDataCardBool("justTPC",true);
  cards_->AddDataCardBool("ccOnly",false);
  cards_->AddDataCardInt("requiredHitsPerTrack",2);
  cards_->AddDataCardDouble("maxEnergy", 30.0*CLHEP::GeV);	//in GeV
  cards_->AddDataCardDouble("testMagField", 1.0);	//in Tesla
}

void myReconstructor::loadDataCards() {
  inputFileName_  	= cards_->fetchValueString("dataFileName");
  outputFileName_  	= cards_->fetchValueString("plotFileName");
  geomFileName_  	= cards_->fetchValueString("geomFileName");
  fiducialCut_ 		= cards_->fetchValueDouble("fiducialCut");
  runs_		  	= cards_->fetchValueInt("numberOfRuns");
  events_		= cards_->fetchValueInt("numberOfEvents");
  verbose_ 		= cards_->fetchValueInt("plotterVerbose");
  pdg_		  	= cards_->fetchValueInt("pdg");
  leptonPdg_		= cards_->fetchValueInt("leptonFlavour");
  justTPC_		= cards_->fetchValueBool("justTPC");
  ccOnly_		= cards_->fetchValueBool("ccOnly");
  hitCountCondition_  	= cards_->fetchValueInt("requiredHitsPerTrack");
  maxEnergy_ 		= cards_->fetchValueDouble("maxEnergy");
  testMagField_ 		= cards_->fetchValueDouble("testMagField");
}

void myReconstructor::initHistos(){

  //declare histograms
  muon_track_length_h1d.SetNameTitle("muon_track_length_h1d", "muon track length");
  muon_track_length_h1d.SetBins(800, 0., 4.0);
  muon_track_length_above1GeV_h1d.SetNameTitle("muon_track_length_above1GeV_h1d", "muon track length");
  muon_track_length_above1GeV_h1d.SetBins(800, 0., 4.0);
  muon_track_length_above2GeV_h1d.SetNameTitle("muon_track_length_above2GeV_h1d", "muon track length");
  muon_track_length_above2GeV_h1d.SetBins(800, 0., 4.0);
  muon_track_length_above3GeV_h1d.SetNameTitle("muon_track_length_above3GeV_h1d", "muon track length");
  muon_track_length_above3GeV_h1d.SetBins(800, 0., 4.0);
  muon_track_length_above4GeV_h1d.SetNameTitle("muon_track_length_above4GeV_h1d", "muon track length");
  muon_track_length_above4GeV_h1d.SetBins(800, 0., 4.0);
  muon_track_length_above5GeV_h1d.SetNameTitle("muon_track_length_above5GeV_h1d", "muon track length");
  muon_track_length_above5GeV_h1d.SetBins(800, 0., 4.0);
  
  muon_track_length_vs_mom_h2d.SetNameTitle("muon_track_length_vs_mom_h2d", "muon track length vs momentum");
  muon_track_length_vs_mom_h2d.SetBins(1000,0.0,10.0,800, 0., 4.0);
  
  truth_cc_muon_trans_spectrum_h1d.SetNameTitle("truth_cc_muon_trans_spectrum_h1d","truth_cc_muon_trans_spectrum_h1d");
  truth_cc_muon_trans_spectrum_h1d.SetBins(1000,0.0,10.0);
  lower_bound_cc_muon_trans_mom_spectrum_h1d.SetNameTitle("lower_bound_cc_muon_trans_mom_spectrum_h1d","lower_bound_cc_muon_trans_mom_spectrum_h1d");
  lower_bound_cc_muon_trans_mom_spectrum_h1d.SetBins(1000,0.0,10.0);
  upper_bound_cc_muon_trans_mom_spectrum_h1d.SetNameTitle("upper_bound_cc_muon_trans_mom_spectrum_h1d","upper_bound_cc_muon_trans_mom_spectrum_h1d");
  upper_bound_cc_muon_trans_mom_spectrum_h1d.SetBins(1000,0.0,10.0);
}

void myReconstructor::reconstruct() {

  //nodes
  int tpcFidCount = 0;
  int vesselCount = 0;
  int innerVesselCount = 0;
  int scintCount = 0;
  int magnetCount = 0;
  int cavityCount = 0;
  int rockCount = 0;
  int otherCount = 0;

  //QELCC count
  int qelccCount = 0;
    
  //exposure
  double exposure = 0;
  //define number of runs
  int runs = runs_;

  std::string filename = inputFileName_;
  std::stringstream ss;

  int totalNuEvents = 0;

  //create a progress bar for file loop and one for entries
  ProgressBar progFiles;
  ProgressBar progEntries;

  badFileCount_ = 0;

  //***loop over runs***
  for(int j = 0;j<runs_;j++){

    //read in the file
    ss << j;
    std::string filename_process = "";

    std::string filename_run;
    if(runs_==1)
    {
      filename_process = filename;
    }
    else
    {
    	filename_process = filename + ss.str() + ".root";
    }
    
    if(verbose_>0)
    {
      std::cout 	<< "\n************************************************************************************************"
                  << "\nFilename being read: " << filename_process
                  << "\n************************************************************************************************" <<std::endl;
    }

    //load the data and stepping files
    TFile datafile(filename_process.c_str());
    
    //if file is bad skip it
    if(datafile.GetNkeys()<1)
    {
      badFileCount_++;
      
      //clear stringstream
      ss.str("");
      continue;
    }

    //clear stringstream 
    ss.str("");

    //set data addresses
    if(!setupTrees(datafile)) return;

    //keep a count of neutrino events
    totalNuEvents +=nuEvents_nEntries;

    int i=0;
    if(verbose_>2 && verbose_<10)
    {
      std::cout << "\n----------------------------------"
	      << "\n    Looping over entries...   "
	      << "\n----------------------------------\n";
    }

    //make a vector to store the badly reconstructed events
    std::vector<ReconData> badEventID;

    // Detector Hits tree used for reconstruction
    if(detectorHitsTree)
    {
      this->setBackTracer(detectorHitsTree);
            
      ///create a tracker object to sort hits into "tracks"///
      Tracker * trackMgr = new Tracker;

      //loop over detector tree
      for(i=0;i<detHits_nEntries;i++)
      {
        detectorHitsTree->GetEntry(i);
      
        //set backtracers for neutrino event truth info
        if(nuHitTree)nuHitTree->GetEntry(backTracer_[0]);
        nuEventTree->GetEntry(backTracer_[1]);
        
        int nuPDG;
        if(nuHit_)nuPDG = nuHit_->getPDG();
        
        double nuE = nuEvent_->getNuEnergy() * (1./CLHEP::GeV);	//should be the same
        if(nuE > maxEnergy_)continue;
        //	std::cout << "\nBacktracer: " << backTracer_[0] << ", " << backTracer_[1] << ", " << backTracer_[2] <<std::endl;
	
        //node name
        std::string nodeName = nuEvent_->getNodeName();
        
        //only consider events in the TPC - skip to next event
        if(nodeName.compare( 0,6, "tpcFid" ) !=0)continue;

        //stop when reach number of required events - use -1 for all events
        if( (events_== tpcFidCount) && (events_>0) )break;

        //get the truth position of the neutrino event
        TLorentzVector eventVertex = nuEvent_->getPosition(); // in mm
        TLorentzVector nuHitPos;
        if(nuHit_)nuHitPos = nuHit_->getPosition();
        else nuHitPos.SetXYZT(0.,0.,0.,0.);

        //check if in fiducial cut - skip to next event if does not fit cut
        if( !inFiducialVolume(nuEvent_,fiducialCut_) )continue;
        
        //increment tpc count
        tpcFidCount++;
        
        //fss vector
        std::vector<ParticleDescrShortRecord> fssVector = nuEvent_->getFssVector();
	
        if( ccOnly_ && nuEvent_->getInteractionType()!=2)continue;	//CC only events

        if(verbose_>=0)
        {
          std::cout <<"\n *** Event Number: " << backTracer_[1];
          nuEvent_->printToStream(std::cout);
        }

        //loop over hits
        HitCollection tpcSdHits = simData_->getTpcFidHits();
        //size of the hit vectors per event
        int tpcHitCount = tpcSdHits.size();

        bool trackSet = false;
        bool goodEvent = false;	//use this to check if particles that leave hits in the event all have range>0
	
        //loop over and read tracks
        bool foundMuon = false;

        /// --- TPC HITS ---- ///
        TpcMomSmear smearer(magField_*(1./CLHEP::tesla),0.0);

        //create a vector of vectors for track information for each track
        std::vector<HitCollection> tpcTrackVector;
        int tpcTrackHitCount = 0;
	
        if(tpcHitCount>0)
        {
          //convert hits into vector of tracks to make it easier to reconstruct each one
          tpcTrackVector = trackMgr->MakeFromHits(tpcSdHits,true);

          //read from the TPC tracks
          for(std::vector< HitCollection >::iterator k = tpcTrackVector.begin();k != tpcTrackVector.end();++k)
          {
            int eventID = backTracer_[1];
            HitCollection tempVtr= *k;
		  
            //number of hits per track
            int numberOfHits = tempVtr.size();
		  
            //total number of hits from all tracks
            tpcTrackHitCount += numberOfHits;

            double initKinEnergy = tempVtr.front().getKinEnergy() * (1./CLHEP::GeV);
            double lastKinEnergy = tempVtr.back().getKinEnergy() * (1./CLHEP::GeV);
            
            double initEnergy = tempVtr.front().getP4().E()* (1./CLHEP::GeV);
            //std::cout <<"\ninitEnergy = " << tempVtr.front().getP4().E()* (1./CLHEP::GeV);
		
            //momentum reconstruction
            double trackTruthTransMomentum = smearer.getTruthTransMomentum(tempVtr);
		  
            //double trackTruthTransMomentum = smearer.getTruthTransMomentum(tempVtr);
            trackTruthTransMomentum *=  (1./CLHEP::GeV);
		  
            //get range of each particle
            double range = 0.;
		  
            //must have at least 2 points to get a range
            if(numberOfHits>1)
            {
              int fssPDG = tempVtr.front().getPDG();
              int fssTrackID = tempVtr.front().getTrackID();

              //range
              range = smearer.getTrackLength(tempVtr);
              range *= (1./CLHEP::m);

              //sagitta
              double truthSag = smearer.getTrackSagitta(tempVtr,trackTruthTransMomentum*CLHEP::GeV);
              double smearSag = smearer.getSmearSagitta(truthSag);
              truthSag *= (1./CLHEP::mm);
              smearSag *= (1./CLHEP::mm);

              TLorentzVector initPos = tempVtr.front().getPosition() * (1./CLHEP::m);
              TLorentzVector lastPos = tempVtr.back().getPosition() * (1./CLHEP::m);

              double eDep = tempVtr.front().getKinEnergy() * (1./CLHEP::MeV) - tempVtr.back().getKinEnergy() * (1./CLHEP::MeV);
			
              bool fullyContained = false;
              if( (tempVtr.back().getTrackLeftVolume()) == 0 )fullyContained = true;
			
              if(fssPDG == 13)foundMuon = true;

              if(verbose_>2 && verbose_<10)
              {
                std::cout << "\n Event: " << eventID
                          << "\n fss PDG " << fssPDG
				  << "\n fss Track Vector Size " << numberOfHits
				  << "\n fss Initial Pos (" <<  initPos.X() << ", " << initPos.Y() << ", " << initPos.Z() << ") m"
				  << "\n fss Final Pos (" <<  lastPos.X() << ", " << lastPos.Y() << ", " << lastPos.Z() << ") m"
				  << "\n fss truth mom " << trackTruthTransMomentum  << " GeV/c"
				  << "\n fss truth sagitta " << truthSag << " mm"
				  << "\n fss smear sagitta " << smearSag << " mm"
				  << "\n fss EDep " << eDep<< " MeV"
				  << "\n fss Initial Kinetic Energy " << initKinEnergy<< " GeV"
				  << "\n fss Final Kinetic Energy " << lastKinEnergy<< " GeV"
				  << "\n fss Range " << range  << " m"
				  << "\n fss eDep/Range " << eDep/(range*100.) << " MeV/cm"
				  << "\n first " << tempVtr.front().getTrackLeftVolume()
				  << "\n last " << tempVtr.back().getTrackLeftVolume()
				  << "\n ---Track fully contained in TPC? " <<  fullyContained << std::endl;
              }
	
              ///------------------- Fill histograms --------------------///
			
              //for particles that show in the TPC
              if(range>0. && (numberOfHits > 5) )
              {
                if(fssPDG==13 && nuPDG==14)
                {
                  muon_track_length_h1d.Fill(range);
                  const double deltaS_overS = TMath::Sqrt(720.0/(hitCountCondition_ + 4.0))*trackTruthTransMomentum*trackTruthTransMomentum*26.7*0.3e-3/(8.0*range*range*testMagField_);
                  const double deltaP_overP = TMath::Sqrt(deltaS_overS*deltaS_overS + 0.01*0.01); // 1% for magnetic field
                  const double deltaP = deltaP_overP*trackTruthTransMomentum;
                  
                  const double lowerP = trackTruthTransMomentum - deltaP;
                  const double upperP = trackTruthTransMomentum + deltaP;
                  
                  lower_bound_cc_muon_trans_mom_spectrum_h1d.Fill(lowerP);
                  upper_bound_cc_muon_trans_mom_spectrum_h1d.Fill(upperP);
                  truth_cc_muon_trans_spectrum_h1d.Fill(trackTruthTransMomentum);
                  
                  if(trackTruthTransMomentum >= 1.0)muon_track_length_above1GeV_h1d.Fill(range);
                  if(trackTruthTransMomentum >= 2.0)muon_track_length_above2GeV_h1d.Fill(range);
                  if(trackTruthTransMomentum >= 3.0)muon_track_length_above3GeV_h1d.Fill(range);
                  if(trackTruthTransMomentum >= 4.0)muon_track_length_above4GeV_h1d.Fill(range);
                  if(trackTruthTransMomentum >= 5.0)muon_track_length_above5GeV_h1d.Fill(range);
                  muon_track_length_vs_mom_h2d.Fill(trackTruthTransMomentum,range);
                }
              }
              goodEvent = true;
            }
            else goodEvent = false;
          
          } // end loop over tpc tracks

        } //end just TPC condition
        
        if(verbose_>1 && verbose_<10)
        {
          std::cout << "\n ================================================================ "
          << "\n Event: " << backTracer_[1]
          << "\n Tpc Hits: " << tpcHitCount
          << "\n Track Hits: " << tpcTrackHitCount
          << "\n Total Neutrino Truth Momentum [GeV/c]: " << nuE
          << "\n ================================================================ " << std::endl;
        }
        
        if(verbose_>2 && verbose_<10)std::cout << "\n*** ENTRY: " << i << ", tpcFidCount = " <<tpcFidCount;
      }//end event loop
      
      delete trackMgr;
    }//end detector tree

    this->deleteTrees();
    datafile.Close();

    //get and print statistics
    if(stats_){
      double tmp_exposure;
      if(verbose_>0){
        std::cout << "\n";
        stats_->printToStream(std::cout);
      }
      map<string,string,ci_less> statsInfo = stats_->getInfoMap();
      map<string,string>::const_iterator it;
	
      for (it = statsInfo.begin(); it != statsInfo.end(); it++)
      {
        std::string statName = it->first;
        std::size_t found = statName.find("Real Exposure");
	  
        if(found!=std::string::npos){
          std::stringstream ss_scientific(it->second);
          ss_scientific.precision(2);
          ss_scientific>> scientific >> tmp_exposure;
          exposure += tmp_exposure;
          //exposure += atoi(it->second.c_str());
        }
      }
    }

    //clear stringstream 
    ss.str("");
    if(verbose_<0)progFiles.showPercent(j,runs,runs,0);
  } //end loop over files

  std::cout << "\n------------------------------------------------------------------------"
		<< "\nTotal Exposure [p.o.t]:                            "
		<< exposure
		<< "\nAverage Exposure [p.o.t]:                          "
		<< exposure/(runs_-badFileCount_)
		<< "\nTotal number of neutrino interactions in target:   "
		<< totalNuEvents
	 	<< "\nTPC Interactions:                                  " 
		<< tpcFidCount
		<< "\n------------------------------------------------------------------------\n"
		<< "\n====================================================="
		<< "\n Files that could not be read: " << badFileCount_
		<< "\n=====================================================";

}

void myReconstructor::write()
{

  std::cout << "\n-----------------------------------------------------------------------------------------"
	    << "\n  Writing to file: " << outputFileName_.c_str()
	    << "\n-----------------------------------------------------------------------------------------\n";

  TFile histofile(outputFileName_.c_str(),"recreate");

  muon_track_length_h1d.SetMarkerColor(kBlue);
  muon_track_length_h1d.Write();
  
  muon_track_length_above1GeV_h1d.Write();
  muon_track_length_above2GeV_h1d.Write();
  muon_track_length_above3GeV_h1d.Write();
  muon_track_length_above4GeV_h1d.Write();
  muon_track_length_above5GeV_h1d.Write();

  muon_track_length_vs_mom_h2d.Write();
  
  lower_bound_cc_muon_trans_mom_spectrum_h1d.Write();
  upper_bound_cc_muon_trans_mom_spectrum_h1d.Write();
  truth_cc_muon_trans_spectrum_h1d.Write();

  std::cout << "\n----------------------------------"
	    << "\n             Closing...           "
	    << "\n----------------------------------\n";

  histofile.Close();
}

int main(int argc, char ** argv) {
  try {
    myReconstructor reconstructor;
    reconstructor.initDataCards();
    reconstructor.getRunTimeArguments(argc, argv);
    reconstructor.loadDataCards();
    if(!reconstructor.loadGeom())return 0;
    reconstructor.initHistos();
    reconstructor.reconstruct();
    reconstructor.write();

  } catch (LbnoException e) {
    std::cout << e.GetLocation() << std::endl;
    std::cout << e.GetDescription() << std::endl;
    return 1;
  }

  return 0;
}
