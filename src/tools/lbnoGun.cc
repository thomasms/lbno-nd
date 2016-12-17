//____________________________________________________________________________
/*!

\class    lbnoGun

\brief    Program to fire single particles using particle gun.
	  Only tracking, no flux or event generation

\author	  Tom Stainer <tstainer \at liv.ac.uk>
          University of Liverpool

\created  Oct 2014
\last update Oct 2014

*/
//____________________________________________________________________________
#include <iostream>

#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>
#include <TRandom3.h>

#include "LbnoProcessorManager.hh"
#include "GeometryLoader.hh"
#include "DataCards.hh"
#include "ParticleGunProcessor.hh"

int main(int argc, char ** argv) {
  try {
    // Create Processor Manage.
    LbnoProcessorManager manager;

    //Initiate geometry loader
    GeometryLoader * geomLoader = new GeometryLoader();

    // Initialize your processors here. The processors will define
    // default values for all runtime parameters.
    ParticleGunProcessor   *g4tracking = new ParticleGunProcessor(geomLoader);

    // Get command line arguments if any. 
    // This will modify the default values for the runtime parameters.
    manager.getRunTimeArguments(argc, argv);
    //manager.loadDataCards();

    // Now add your processors to the manages. 
    manager.addProcessor(g4tracking);

    //reinitiate the geometry loader
    //geomLoader->initialize(eventLoader->get);

    manager.printToStream(cout);

    // Run your processors.
    manager.go();

    // Store the result.
    manager.write();

   //close the geometry
   //geomLoader->close();

  } catch (LbnoException e) {
    std::cout << e.GetLocation() << std::endl;
    std::cout << e.GetDescription() << std::endl;
    return 1;
  }
  return 0;
}


