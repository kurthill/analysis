int run_macro( 
	      //std::string infile = "/phenix/hhj/kurthill/pFlowJets/data/pythia_pthatmin50_calo_dst_0.root",
	      std::string infile = "G4sPHENIX.root",
	      std::string outfile = "test.root"
	       )
{
  
  gSystem->Load("libg4dst.so");
  gSystem->Load("libfun4all.so");
  gSystem->Load("libphfield_io.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libphhepmc.so");
  gSystem->Load("libg4testbench.so");
  //gSystem->Load("libg4hough.so");
  //gSystem->Load("libcemc.so");
  gSystem->Load("libg4eval.so");

  gSystem->Load("libpFlowTreeMaker.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  recoConsts *rc = recoConsts::instance();

  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  hitsin->fileopen( infile );
  se->registerInputManager(hitsin);

  pFlowTreeMaker *tm = new pFlowTreeMaker( outfile );
  se->registerSubsystem( tm );

  se->run();

  se->End();
  std::cout << "All done" << std::endl;
  delete se;

  gSystem->Exit(0);
}
