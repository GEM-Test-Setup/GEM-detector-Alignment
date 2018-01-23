void Run_MC(int run_num =1000){

  TDatime *time1 = new TDatime();
  int random_seed = time1->Get();
  gROOT->ProcessLine(".L TrackGen.C");
  
  TRandom *random_gen = new TRandom(random_seed);

  for(int irun=0; irun<run_num; irun++){
    random_seed = random_gen->Uniform(10*run_num);R
    TrackGen(random_seed);
    gSystem->Exec("root -b -q -l RegressTrack.C");
  }
}
