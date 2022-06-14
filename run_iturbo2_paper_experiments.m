function run_iturbo2_paper_experiments(datafile, test_run)
%%  MATLAB script to run multiple experiments of iTURBO2 that are included in the paper
%% to run: run_iturbo2_paper_experiments('data/iTURBO2_input_data.xlsx', false)
%
% datafile  :   INPUT FILE WITH REQUIRED DATA
% test_run  :   true = fewer simulations -> to run faster
%               false = run as for paper
%
% Example call: run_iturbo2_paper_experiments('data/iTURBO2_input_data.xlsx', false)
%
%   1) The six idealized shapes
%   2) Artificial sinusoidal isotopic change with 20, 40, and 100 kyr periode
%   3) Observed ash profiles
%   4) Isotopic change following Laskar solution
%   5) Dissolution experiments
%   6) Isotope & Abu-change for cold and warm species
%   7) Call the function to excecute the PETM experiment
%
%   below it is specified where the data of the isotopes shapes and ash observations can be found in the excel file
%   and iTURBO2 is called and the results are plotted
%   here 3 different bioturbation depths are compared (i.e. mxl cm (as specified in datafile) 10 cm and 20 cm))
%   or for the dissolution experiment: three different dissolution scenrios are simulated


% 
tmp_mutlab = version('-release');
str_mutlab = tmp_mutlab(1:4);
par_mutlab = str2num(str_mutlab);

if (par_mutlab > 2019)
    disp([' ']);
    disp(['!!!!!!!!!!!!!!!!!!!!!']);
    disp(['>>> Warning: The code might not work with your matlab version ...']);
    disp(['>>> The code has been created with version 2019b  ...']);
    disp(['!!!!!!!!!!!!!!!!!!!!!']);
    disp([' ']);

end    


% add path to plotting utils
addpath('./utils');

%%  Set default plotting as 'false' here
%  For each experiment one can specify which plots should be produced
settings.plot_iso_spec1 = false;            	% plot isotopes for species 1 only
settings.plot_isotope_both_species = false;  	% plot isotops for both species
settings.plot_abu_iso = false;                  % plot abundance and isotope figures for both species
settings.plot_just_mean = false;                % only plot the mean result, not every single experiment in grey
settings.plot_ash_expl = false;                 % plot ash example

settings.sprectral_ana = false;                 % signal needs to be flipped for my spectral analysis

disp([' ']);
if test_run
    disp(['>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>']);
    disp(['>>> Running all experiments of the paper but using fewer model simulations ...']);
    disp(['>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>']);
    disp([' ']);
else
    disp(['>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>']);
    disp(['>>> Running all experiments of the paper ...']);
    disp(['>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>']);
    disp([' ']);
end

if false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%               1) The six idealized shapes                   %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' ']);
disp(['>>> Running "The six idealized shapes" ...']);


carriers = 10000;     	% to be measured    [10000]
Exps = 100;             % Experiments to run; each one and the mean is plotted  [100]

if test_run
    datafile = 'data/test_run/iTURBO2_input_data_test_run.xlsx';
    carriers = 100;     	% to be measured    [10000]
    Exps = 2;         	% Experiments to run; each one and the mean is plotted  [100]
end

settings.SetXaxis_lim = true;   % Set lim for Xaxis
settings.Xaxis_lim = 200;       % value for Xaxis lim

% What do you want to plot?
settings.plot_iso_spec1 = true;              % plot isotopes for species 1 only

TM = [false, true];     % simulate homogeneous mixing and using a transition matrix (specify TM in iturbo2script_multipleSims_3zbio)

for i=1:length(TM)
    data=xlsread(datafile,'shapes','C4:F263');
    iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'Fig_IdealizedShapes_point_event',TM(i), settings)
    close all
    
    data=xlsread(datafile,'shapes','H4:K263');
    iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'Fig_IdealizedShapes_multiple_points',TM(i), settings)
	close all

    data=xlsread(datafile,'shapes','M4:P263');
    iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'Fig_IdealizedShapes_step_sequence',TM(i), settings)
	close all
    
    data=xlsread(datafile,'shapes','R4:U263');
    iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'Fig_IdealizedShapes_gradual_change', TM(i), settings)
	close all
    
    data=xlsread(datafile,'shapes','W4:Z263');
    iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'Fig_IdealizedShapes_stepwise_down', TM(i), settings)
	close all
    
    data=xlsread(datafile,'shapes','AB4:AE263');
    iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'Fig_IdealizedShapes_40kyrs',TM(i), settings)
	close all
end

settings.plot_iso_spec1 = false;              % plot isotopes for species 1 only

disp(['>>> Done with "The six idealized shapes" ...']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%         2) Artificial sinusoidal isotopic change            %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Artificial signal with 20, 40 and 100 kyr periode - 1.1 Mio years
disp([' ']);
disp(['>>> Running "Artificial sinusoidal isotopic change" ...']);

carriers = 10000;     	% to be measured [10000]
Exps = 100;              % Experiments to run; each one and the mean is plotted [100]

if test_run
    datafile = 'data/test_run/iTURBO2_input_data_test_run.xlsx';
    carriers = 50;     	% to be measured    [10000]
    Exps = 2;             % Experiments to run; each one and the mean is plotted  [100]
end

settings.plot_iso_spec1 = true;              % plot isotopes for species 1 only

settings.SetXaxis_lim = false;   % Set lim for Xaxis
settings.Xaxis_lim = 1000;       % value for Xaxis lim

% What do you want to plot?
settings.plot_iso_spec1 = true;              % plot isotopes for species 1 only

TM = false;

data=xlsread(datafile,'Artificial_SinusSignal','C4:F1103');
iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'Fig_ArtificialSinusSignal',TM, settings)
close all

settings.plot_iso_spec1 = false;              % plot isotopes for species 1 only

disp(['>>> Done with "Artificial sinusoidal isotopic change" ...']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%                 3) Observed ash profiles                    %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' ']);
disp(['>>> Running "Observed ash profiles" ...']);

carriers = 10000;     	% to be measured  [10000]
Exps = 100;        	% Experiments to run; each one and the mean is plotted [100]
%% ASH experiment run with input from iTURBO2_input_ash_experiment.xlsx
datafile2 = 'data/iTURBO2_input_ash_experiment.xlsx';

if test_run
    carriers = 100;     	% to be measured    [10000]
    Exps = 200;             % Experiments to run; each one and the mean is plotted  [100]
    datafile2 = 'data/test_run/iTURBO2_input_ash_experiment_test_run.xlsx';
end

data=xlsread(datafile2,'ash_data','C4:F73');

TM = [false, true];     % simulate homogeneous mixing and using a transition matrix

for i=1:length(TM)

    Observation = 1;
    iturbo2script_3zbio_ASH(data, carriers, Exps, 'Ash_experiment_0.5cmkyr', TM(i), Observation)
%    close all
    Observation = 2;
    iturbo2script_3zbio_ASH(data, carriers, Exps, 'Ash_experiment_2.0cmkyr', TM(i), Observation)
%    close all
end

disp(['>>> Done with "Observed ash profiles" ...']);

if false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%         4) Isotopic change following Laskar solution        %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the bioturbated signals is saved in a .mat file under ./data/mat
% file name : ['3zbio_',expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_homogen']
% and can then be used with a spectral analysis toolbox

% As the core record is very long it will take a pretty long time to measure all 10000 particles
% consider using fewer particles (e.g. carriers here, and also ABU in the excel file)

disp([' ']);
disp(['>>> Running "Isotopic change following Laskar solution" ...']);

carriers = 10000;     	% to be measured
Exps = 10;              % Experiments to run; each one and the mean is plotted

if test_run
    datafile = 'data/test_run/iTURBO2_input_data_test_run.xlsx';
    carriers = 10;     	% to be measured    [10000]
    Exps = 1;             % Experiments to run; each one and the mean is plotted  [100]
end

settings.SetXaxis_lim = false;      % Set lim for Xaxis
settings.Xaxis_lim = 4000;          % value for Xaxis lim

settings.sprectral_ana = true;      % just flip the bioturbated before it is saved - not sure if this is correct here

% What do you want to plot?
settings.plot_iso_spec1 = true;              % plot isotopes for species 1 only

TM = false;

% Get Laskar 4.1 Mio signal flipped - 19.02.2021
data=xlsread(datafile,'Laskar_solution','C4:F4103');
iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'Fig_Laskar_solution',TM, settings)
close all

settings.plot_iso_spec1 = false;              % plot isotopes for species 1 only

disp(['>>> Done with "Isotopic change following Laskar solution" ...']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%                 5) Dissolution experiments                  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' ']);
disp(['>>> Running "Dissolution experiments" ...']);

carriers = 1000;     	% to be measured [1000]
Exps = 100;        	% Experiments to run; each one and the mean is plotted [100]
%% dissolution experiment run with input from iTURBO2_input_dissolution_experiment.xlsx
datafile3 = 'data/iTURBO2_input_dissolution_experiment.xlsx';

if test_run
    datafile3 = 'data/test_run/iTURBO2_input_dissolution_experiment_test_run.xlsx';
    carriers = 1000;     	% to be measured    [10000]
    Exps = 2;             % Experiments to run; each one and the mean is plotted  [100]
end


TM = false;

% % % % % %% point event
plot_min = true; % just for plotting purposes -- to make sure the min point is plotted

data=xlsread(datafile3,'nochange','C4:F263');
data50=xlsread(datafile3,'50perc_diss','C4:F263');
data90=xlsread(datafile3,'90perc_diss','C4:F263');
data100=xlsread(datafile3,'100perc_diss','C4:F263');
iturbo2script_multipleSims_dissolution(data, data50, data90, data100, carriers, Exps, 'Dis_pointevent', TM, plot_min)
close all

% %%step sequence
plot_min = false; % set to false again as it is not a single min value
data=xlsread(datafile3,'nochange','H4:K263');
data50=xlsread(datafile3,'50perc_diss','H4:K263');
data90=xlsread(datafile3,'90perc_diss','H4:K263');
data100=xlsread(datafile3,'100perc_diss','H4:K263');
iturbo2script_multipleSims_dissolution(data, data50, data90, data100, carriers, Exps, 'Dis_step_sequence', TM, plot_min)
close all

% %% gradual change
plot_min = true; % set to true again
data=xlsread(datafile3,'nochange','O4:R263');
data50=xlsread(datafile3,'50perc_diss','O4:R263');
data90=xlsread(datafile3,'90perc_diss','O4:R263');
data100=xlsread(datafile3,'100perc_diss','O4:R263');
iturbo2script_multipleSims_dissolution(data, data50, data90, data100, carriers, Exps, 'Dis_gradual_change', TM, plot_min)
close all

disp(['>>> Done with "Dissolution experiments" ...']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%     6) Isotope & Abu-change for cold and warm species       %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' ']);
disp(['>>> Running "Isotope & Abu-change for cold and warm species" ...']);

carriers = 20;     	% to be measured [20]
Exps = 100;        	% Experiments to run; each one and the mean is plotted [100]

if test_run
    datafile = 'data/test_run/iTURBO2_input_data_test_run.xlsx';
    carriers = 20;     	% to be measured    [10000]
    Exps = 6;             % Experiments to run; each one and the mean is plotted  [100]
end

TM = false;

% % multistep change with 1000 total abundance
data=xlsread(datafile,'iso_abu_change_2Species','C4:F263');
iturbo2script_multipleSims(data, carriers, Exps, 'Cold_vs_warm_species', TM)
close all

disp(['>>> Done with "Isotope & Abu-change for cold and warm species" ...']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%   7) Call the function to plot the PETM ODP 690 results     %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' ']);
disp(['>>> Running "PETM OPD 690 experiments" ...']);

plot_paper_results = true;
numb_carriers = 4;
% pick and measure n = 4 carriers
run_iturbo2_exps_PETM690('data/iTURBO2_input_PETM_690_Exp.xlsx' , 'PETM_results', plot_paper_results, numb_carriers)
close all

% pick and measure n = 12 carriers
numb_carriers = 12;
run_iturbo2_exps_PETM690('data/iTURBO2_input_PETM_690_Exp.xlsx' , 'PETM_results', plot_paper_results, numb_carriers)
close all

% Note: Pick new carriers to measure (and decide how many) by setting
% plot_paper_results = false and selecting the desired number of carriers (numb_carriers)
% For example:
plot_paper_results = false;
numb_carriers = 20;
run_iturbo2_exps_PETM690('data/iTURBO2_input_PETM_690_Exp.xlsx' , 'PETM_results', plot_paper_results, numb_carriers)

disp(['>>> Done with "PETM OPD 690 experiments" ...']);

% final messages
disp([' ']);
disp(['------------------------------------------------------------']);
disp(['   Congratulations! All experiments are done!  ']);
disp(['   The .eps files are saved in the directory "output" ']);
disp(['------------------------------------------------------------']);
disp([' ']);

end