# DoEasyTree

Sort code which takes standard AnalysisTrees created in GRSISort and creates higher level 'Easy' trees which can be analyzed much quicker and with more control than low-level stuff. Doppler corrected energies, excitation energies and theoretical energies are calculated for each reaction type so that in-beam calibrations can be carried out and data can be analyzed using simple draw commands. 

In order to run the code, TSharcAnalysis must be installed and compiled and an input file (.txt) must be included.

The input configuration file specifies;
  - Beam nucleus inc. energy
  - Calibration [.cal] file and coefficients
  - Energy cuts and particle multiplicities
  - Particle cuts [.root] file 
  - Configuration parameters such as target position and thickness
  - Other awesome stuff I can't think of right now

The output 'Easy' trees can be further analyzed using the TTigressAnalysis, MakeEasyMats, TSharcAnalysis, AngularDistribution and FrescoAnalysis packages.

grsisort -l
TChain *chain = new TChain("EasyTree")
chain->Add("easy*.root")
chain->Draw("exc:eadd>h(4000,0,4000,1000,-2000,8000)","type


#Notes
  -Requires TSharcAnalysis to run. 
  -Sharc timing is not included (was not working when this code was produced). To include, just add another branch into the trees e.g. sharctime.
  -Trifoil is also not included, but can easily be added.
