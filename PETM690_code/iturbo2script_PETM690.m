function [biopart1, bioabu] = iturbo2script_PETM690(data, nocarriers, noexps, expname, color_points)
%% MATLAB script to run iTURBO2 experiment with input data to reproduce PETM 690 observations

% data = matrix of required data (age, mxl, abu, iso) - explanation see below
% nocarriers = number of carriers to be measured
% noexps = number of experiments for SA
% expname = experiment name
% color_points = defines the color of the points and lines

% age = age of sediment layer down core
% mxl = series of mixed layer thicknesses (zbio) down core
% abu = series of abundances of carrier type 1 down core
% iso = original isotope signature of both carrier  types 1 and 2


age   = data(:,1);
age_flipud = flipud(age);
mxl1   = data(:,2);

abu   = data(:,3);
iso   = data(:,4);
lngth = length(data(:,1));

numb  = nocarriers;     % number of carriers to be measured
exps = noexps;          % number of different experiments

%%
for i = 1:exps
    [oriabu(i,:,:),bioabu(i,:,:),oriiso(i,:,:),bioiso(i,:,:),biopart1(i,:,:)] = iturbo2_PETM690(abu,iso,mxl1,numb);
end
% normalize bioabu:
bioabu_norm = bioabu./nocarriers;

% variable for mean results of mxl1
mean_bioabu1_mxl1 = zeros(1,lngth);
mean_bioabu2_mxl1 = zeros(1,lngth);
mean_bioiso1_mxl1 = zeros(1,lngth);
mean_bioiso2_mxl1 = zeros(1,lngth);

%%
mxltext = num2str(mean(mxl1));

numbtxt = num2str(numb,2);
expstxt = num2str(exps,2);
abutxt = num2str(max(abu),2);

set(0,'DefaultAxesFontSize',16)

%%  Plot isotope plots only and just for species 1 (no abundance change here)
figure, hold on
for i = 1:exps
    % mxl1
    plot(bioiso(i,:,1),age_flipud, 'Color', [1.0 0.8 0.8],'Linewidth',1.5)
    mean_bioiso1_mxl1 = mean_bioiso1_mxl1+bioiso(i,:,1);
    for j=1:numb
        plot(biopart1(i,:,j),age_flipud, 's', 'MarkerEdgeColor', color_points)
    end
    
end
plot(oriiso(1,:,1),age_flipud,'k','Linewidth',2.0) % plot original iso for species 1
% mxl1
mean_bioiso1_mxl1 = mean_bioiso1_mxl1/exps;
mean_bioiso2_mxl1 = mean_bioiso2_mxl1/exps;
plot(mean_bioiso1_mxl1,age_flipud, '-', 'Color', color_points ,'Linewidth',2.0)

set(gca,'YDir','Reverse','XGrid','On','YGrid','On','Box','On', 'XLim',[-2.0,4.0], 'YLim',[170.4,171.0]);%,'YTick',[1.0 1.5 2.0 2.5 3.0])
ylabel('Meters below seafloor (mbsf) ');
xlabel('\delta^{13}C');
titletxt = ['Isotopes of Carriers 1+2, ',mxltext,...
    ' cm Mixed Layer, ',numbtxt,' Carriers'];
title(titletxt)

%% Do not print these figures
% printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps'];
% save(['data/mat/',printfilename,'.mat'],'printfilename', 'lngth','bioiso','bioiso2','bioiso3', 'oriiso', 'mean_bioiso1_mxl1', 'mean_bioiso1_mxl2', 'mean_bioiso1_mxl3','expname', 'exps', 'bioabu','bioabu2','bioabu3', 'oriabu', 'oriabu2', 'oriabu3', 'mean_bioabu1_mxl1', 'mean_bioabu1_mxl2', 'mean_bioabu1_mxl3')
% print('-depsc', ['output/',printfilename]);   % save figure in extra output folder



