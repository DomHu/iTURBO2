function run_iturbo2_simple_experiments(datafile)
%%  MATLAB script to run some simple experiments of iTURBO2
%%  to run: run_iturbo2_paper_experiments('data/iTURBO2_input_data.xlsx')
%
% datafile  :   INPUT FILE WITH REQUIRED DATA

%%  Set default plotting as 'false' here
%  For each experiment one can specify which plots should be produced
settings.plot_iso_spec1 = false;            	% plot isotopes for species 1 only
settings.plot_isotope_both_species = false;  	% plot isotops for both species
settings.plot_abu_iso = false;                  % plot abundance and isotope figures for both species
settings.plot_just_mean = false;                % only plot the mean result, not every single experiment in grey
settings.plot_ash_expl = false;                 % plot ash example

settings.sprectral_ana = false;                 % signal needs to be flipped for my spectral analysis



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  1) Two idealized shapes with isotope and abundance change  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Just run a single simulation
Exps = 1;               
% Only measure 20 forams
carriers = 20;         


settings.SetXaxis_lim = true;   % Set lim for Xaxis
settings.Xaxis_lim = 200;       % value for Xaxis lim


% What do you want to plot?
settings.plot_abu_iso = true;                  % plot abundance and isotope figures for both species BUT HERE ONLY FOR ONE ZBIO (the one specified in the excel file!)

TM = true;

data=xlsread(datafile,'shapes_low_abu','C4:F263');
iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'point_event',TM, settings)

data=xlsread(datafile,'shapes_low_abu','H4:K263');
iturbo2script_multipleSims_3zbio(data, carriers, Exps, 'step_sequence',TM, settings)


