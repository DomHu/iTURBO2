
################################################
### iTURBO2 toolbox for MATLAB
################################################

* The original model by Trauth (2013) was adapted for Hülse et al. (2022)

Hülse et al. (2022). Assessing the impact of bioturbation on sedimentary isotopic records: a model study. Paleoceanography, XXX, XXX

Trauth, M. H. (2013). TURBO2: A MATLAB simulation to study the effects of bioturbation on paleoceanographic time series. Computers & Geosciences, 61:1–10.


################################################
### How to run experiments shown in the paper
################################################

** To run all experiments with one function, call the routine run_iturbo2_paper_experiments.m /
** To run just a single experiment select the respective function-calls in run_iturbo2_paper_experiments.m

* Run all Hülse et al. (2022) experiments but with fewer model iterations so the experiments are done faster:
run_iturbo2_paper_experiments('data/iTURBO2_input_data.xlsx', true)

* Run all experiments with the same model iterations as in Hülse et al. (2022):
run_iturbo2_paper_experiments('data/iTURBO2_input_data.xlsx', false)

* output figures: Are saved in folder ./output and have the general structure:
expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps'
with
expname	:	name as specified by the user (see e.g. in file run_iturbo2_paper_experiments.m)
abutxt	:	max number of carriers in the entire core
numbtxt	:	number of carriers to be measured
expstxt	:	number of different experiments

* The MATLAB results of the specific experiments are saved in the folder ./output/mat and have the same general structure as the figures


################################################
### Other / simple iTURBO2 example experiments
################################################

* run_iturbo2_simple_experiments.m
run_iturbo2_simple_experiments('data/iTURBO2_input_data.xlsx')
Performs two simple iTURBO2 experiments using abudance and isotope changes for only one model simulation.


################################################
### Description of the MATLAB files:
################################################

run_iturbo2_paper_experiments.m			- MATLAB script to run the experiments presented in the paper
run_iturbo2_exps_PETM690.m			- MATLAB script to run the PETM OPD 690 experiments presented in the paper (is also called from run_iturbo2_paper_experiments.)
run_iturbo2_simple_experiments.m		- MATLAB script to run simple iTURBO2 experiments where only one simulation is run (isotope & abundance changes simultaneously)

iturbo2_plus_TM.m       		  	- MATLAB code of iTURBO2, now including the transition matrix (TM) approach
iturbo2script_multipleSims.m		  	- MATLAB script to run iTURBO2 using multiple simulations (i.e. with different random mixing realisations) for one biotrbation depths
iturbo2script_multipleSims_3zbio.m 	  	- MATLAB script to run iTURBO2 using multiple simulations for three different biotrbation depths
iturbo2script_4zbio_ASH.m			- MATLAB script to run iTURBO2 to simulate the observed ash layers using four different biotrbation depths
iturbo2script_3zbio_ASH.m			- MATLAB script to run iTURBO2 to simulate the observed ash layers using three different biotrbation depths
iturbo2script_multipleSims_dissolution.m	- MATLAB script to run iTURBO2 using multiple simulations for three different dissolution scenarios

./PETM690_code/iturbo2_PETM690.m			- MATLAB code of iTURBO2 for the PETM 690 experiment (does not include the TM approach)
./PETM690_code/iturbo2script_PETM690.m		- MATLAB script to run iTURBO2 to simulate the observed PETM ODP690 record

################################################
### Description of the data files that contain information for the mixing experiments:
################################################

* in ./data:
iTURBO2_input_data.xlsx				- Information for experiments (artificial shapes, artificial sinus signal, Laskar solution, isotope & abundace change for cold and waRM species) 
iTURBO2_input_ash_experiment.xlsx		- Information for the ash experiments
iTURBO2_input_dissolution_experiment.xlsx	- Information for the dissolution experiments
iTURBO2_input_PETM_690_Exp.xlsx			- Information for the PETM 690 experiments

Excel files as in ./data just with fewer particles so a test run can be performed faster

* in ./data/PETM_690:
PETM690.mat								- PETM OPD 690 observations
zbio_10cm_thermo_0.1_and_1.0_percent_NoDeepMix_v2.mat			- iTURBO2 results shown in the paper for measuring n=4 carriers
zbio10cm_12carriers_Thermo_0.1_and_1.0_percent_IntermedDeepMix.mat	- iTURBO2 results shown in the paper for measuring n=12 carriers

################################################
### Description of other folders:
################################################

* transmtx:	- The transition matrices generated and used for the local (upward) mixing for bioturbation depths of 5, 7, 9, 10, 11, 13 and 20 cm
		- New, user-specified transition matrices have to be saved in this folder and need to be read and used by the iTURBO2 model in the file iturbo2_plus_TM.m. 
* utils:	- MATLAB utility functions for plotting purposes
* output: 	- Result figures are saved here as .eps files. 
* output/mat:	- MATLAB .mat files are saved here that record the results of the specific experiments so they can be recreated. 


################################################
### Licensing
################################################
Most of the routines included are MIT licensed. But please read details in individual files, as it includes some codes that are not authored by us.


