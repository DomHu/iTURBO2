function iturbo2script_multipleSims_dissolution(data, data50, data90, data100, nocarriers, noexps, expname, TM, plot_min)
%% MATLAB script to run TURBO2 experiment with input data
%% RUN three different dissolution experiments AND PLOT results in same Fig.

% data = matrix of required data (age, mxl, abu, iso) without dissoution (abundance change)
% data50 = matrix of required data with 50% dissoution (50% abundance change)
% data90 = matrix of required data with 90% dissoution (90% abundance change)
% data100 = matrix of required data with 100% dissoution (100% abundance change)
% nocarriers = number of carriers to be measured
% noexps = number of experiments for SA
% expname = experiment name
% TM = true or false for using trnasition matrix or homogeneous mixing
% plot_min = true/false; just for plotting purposes -- to make sure the min point is plotted

% age = age of sediment layer down core
% mxl = series of mixed layer thicknesses (zbio) down core
% abu = series of abundances of carrier type 1 down core
% iso = original isotope signature of both carrier  types 1 and 2

c = @cmu.colors; % shortcut function handle

age   = data(:,1);  % here asme for all data-files
mxl1   = data(:,2); % here asme for all data-files
abu   = data(:,3);
abu50   = data50(:,3);
abu90   = data90(:,3);
abu100   = data100(:,3);
iso   = data(:,4); % here asme for all data-files


lngth = length(data(:,1));

numb  = nocarriers;     % number of carriers to be measured
exps = noexps;          % number of different experiments

%%
for i = 1:exps
    i
    [oriabu(i,:,:),bioabu(i,:,:),oriiso(i,:,:),bioiso(i,:,:)] = iturbo2_plus_TM(abu,iso,mxl1,numb, TM);
    [oriabu2(i,:,:),bioabu2(i,:,:),oriiso2(i,:,:),bioiso2(i,:,:)] = iturbo2_plus_TM(abu50,iso,mxl1,numb, TM);
    [oriabu3(i,:,:),bioabu3(i,:,:),oriiso3(i,:,:),bioiso3(i,:,:)] = iturbo2_plus_TM(abu90,iso,mxl1,numb, TM);
    [oriabu4(i,:,:),bioabu4(i,:,:),oriiso4(i,:,:),bioiso4(i,:,:)] = iturbo2_plus_TM(abu100,iso,mxl1,numb, TM);
end
%%
% variable for mean results of no dissolution
mean_bioabu1_diss0 = zeros(1,lngth);
mean_bioabu2_diss0 = zeros(1,lngth);
mean_bioiso1_diss0 = zeros(1,lngth);
mean_bioiso2_diss0 = zeros(1,lngth);
% variable for mean results of 50% dissolution
mean_bioabu1_mxl2 = zeros(1,lngth);
mean_bioabu2_mxl2 = zeros(1,lngth);
mean_bioiso1_mxl2 = zeros(1,lngth);
mean_bioiso2_mxl2 = zeros(1,lngth);
% variable for mean results of 90% dissolution
mean_bioabu1_mxl3 = zeros(1,lngth);
mean_bioabu2_mxl3 = zeros(1,lngth);
mean_bioiso1_mxl3 = zeros(1,lngth);
mean_bioiso2_mxl3 = zeros(1,lngth);
% variable for mean results of 100% dissolution
mean_bioabu1_diss100 = zeros(1,lngth);
mean_bioabu2_diss100 = zeros(1,lngth);
mean_bioiso1_diss100 = zeros(1,lngth);
mean_bioiso2_diss100 = zeros(1,lngth);

%%
mxltext = num2str(mean(mxl1));
numbtxt = num2str(numb,2);
expstxt = num2str(exps,2);
abutxt = num2str(max(abu),2);


set(0,'DefaultAxesFontSize',16)

%%  Plot isotope plots only and just for species 1 (no abundance change here)
fig1 = figure;
hold on
for i = 1:exps
    % dissolution 90%
    plot(1:lngth,bioiso3(i,:,1), 'Color', [0.8 0.8 1.0],'Linewidth',1.5)
    mean_bioiso1_mxl3 = mean_bioiso1_mxl3+bioiso3(i,:,1);
    % dissolution 50%
    plot(1:lngth,bioiso2(i,:,1), 'Color', [0.8 1.0 0.8],'Linewidth',1.5)
    mean_bioiso1_mxl2 = mean_bioiso1_mxl2+bioiso2(i,:,1);
    % dissolution 0%
    plot(1:lngth,bioiso(i,:,1), 'Color', [1.0 0.8 0.8],'Linewidth',1.5)
    mean_bioiso1_diss0 = mean_bioiso1_diss0+bioiso(i,:,1);
    % dissolution 100%
    plot(1:lngth,bioiso4(i,:,1), 'Color',c('carrot orange'),'Linewidth',1.5)
    mean_bioiso1_diss100 = mean_bioiso1_diss100+bioiso4(i,:,1);
end
plot(1:lngth,oriiso(1,:,1),'k','Linewidth',2.0) % plot original iso for species 1
% mxl3
mean_bioiso1_mxl3 = mean_bioiso1_mxl3/exps;
plot(1:lngth,mean_bioiso1_mxl3, '-b','Linewidth',2.0)
% mxl2
mean_bioiso1_mxl2 = mean_bioiso1_mxl2/exps;
plot(1:lngth,mean_bioiso1_mxl2, '-g','Linewidth',2.0)
% mxl1
mean_bioiso1_diss0 = mean_bioiso1_diss0/exps;
plot(1:lngth,mean_bioiso1_diss0, '-r','Linewidth',2.0)
% dissolution 100%
mean_bioiso1_diss100 = mean_bioiso1_diss100/exps;
plot(1:lngth,mean_bioiso1_diss100, '-','Color',c('deep carrot orange'),'Linewidth',2.0)

set(gca,'YDir','Reverse','XGrid','On','YGrid','On','Box','On', 'XLim',[0,200], 'YLim',[1.0,3.0],'YTick',[1.0 1.5 2.0 2.5 3.0])
xlabel('Core depth (cm) ');
ylabel('\delta^{18}O');

printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps'];
save(['output/mat/',printfilename,'.mat'],'printfilename', 'lngth','bioiso','bioiso2','bioiso3','bioiso4', 'oriiso', 'mean_bioiso1_diss0', 'mean_bioiso1_mxl2', 'mean_bioiso1_mxl3', 'mean_bioiso1_diss100','expname', 'exps', 'bioabu','bioabu2','bioabu3','bioabu4', 'oriabu', 'oriabu2', 'oriabu3', 'oriabu4', 'mean_bioabu1_diss0', 'mean_bioabu1_mxl2', 'mean_bioabu1_mxl3', 'mean_bioabu1_diss100')
print(fig1, '-depsc', ['output/',printfilename]);   % save figure in extra output folder

%% Plot abundance of species 1
fig2 = figure;
hold on

for i = 1:exps
    % 0% dissolution
    plot(1:lngth,bioabu(i,:,1), 'Color', [1.0 0.8 0.8],'Linewidth',1.5)
    mean_bioabu1_diss0 = mean_bioabu1_diss0+bioabu(i,:,1);
    % 50% dissolution
    plot(1:lngth,bioabu2(i,:,1), 'Color', [0.8 1.0 0.8],'Linewidth',1.5)
    mean_bioabu1_mxl2 = mean_bioabu1_mxl2+bioabu2(i,:,1);
    % 90% dissolution
    plot(1:lngth,bioabu3(i,:,1), 'Color', [0.8 0.8 1.0],'Linewidth',1.5)
    mean_bioabu1_mxl3 = mean_bioabu1_mxl3+bioabu3(i,:,1);
    % 100% dissolution
    plot(1:lngth,bioabu4(i,:,1), 'Color',c('carrot orange'),'Linewidth',1.5)
    mean_bioabu1_diss100 = mean_bioabu1_diss100+bioabu4(i,:,1);
end

% 100% dissolution
mean_bioabu1_diss100 = mean_bioabu1_diss100/exps;
plot(1:lngth,mean_bioabu1_diss100, 'Color',c('deep carrot orange'),'Linewidth',2.0)
% 90% dissolution
mean_bioabu1_mxl3 = mean_bioabu1_mxl3/exps;
plot(1:lngth,mean_bioabu1_mxl3, '-b','Linewidth',2.0)
% 50% dissolution
mean_bioabu1_mxl2 = mean_bioabu1_mxl2/exps;
plot(1:lngth,mean_bioabu1_mxl2, '-g','Linewidth',2.0)
% 0% dissolution
mean_bioabu1_diss0 = mean_bioabu1_diss0/exps;
plot(1:lngth,mean_bioabu1_diss0, '-r','Linewidth',2.0)

% plot one of the original abu 0% dissolution
line_fewer_markers(1:lngth,oriabu(1,:,1), 20,'--ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',6);
% plot one of the original abu 50% dissolution
line_fewer_markers(1:lngth,oriabu2(1,:,1),20,'--ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',6);
% plot one of the original abu 90% dissolution
line_fewer_markers(1:lngth,oriabu3(1,:,1), 20,'--ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'MarkerSize',6);
% plot one of the original abu 100% dissolution
line_fewer_markers(1:lngth,oriabu4(1,:,1), 20,'--ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',c('deep carrot orange'),...
    'MarkerSize',6);

% plot minimum number of particles
if(plot_min)
    % find min
    idmin50perc = find(oriabu2(1,:,1) == min(oriabu2(1,:,1)));
    idmin90perc = find(oriabu3(1,:,1) == min(oriabu3(1,:,1)));
    idmin100perc = find(oriabu4(1,:,1) == min(oriabu4(1,:,1)));
    
    plot(idmin50perc, oriabu2(1,idmin50perc,1),'o-', 'LineWidth',2, 'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',6)
    plot(idmin90perc, oriabu3(1,idmin90perc,1),'o-', 'LineWidth',2, 'MarkerEdgeColor','k',...
        'MarkerFaceColor','b',...
        'MarkerSize',6)
    plot(idmin100perc, oriabu4(1,idmin100perc,1),'o-', 'LineWidth',2, 'MarkerEdgeColor','k',...
        'MarkerFaceColor',c('deep carrot orange'),...
        'MarkerSize',6)
end

set(gca,'XGrid','On','YGrid','On', 'XLim',[0,200], 'YLim',[0,1000], 'YTick',[0 200 400 600 800 1000])
xlabel('Core depth (cm) ');
ylabel('Number of Particles');
title('Abundance of Species 1')

printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_ABU'];
print(fig2, '-depsc', ['output/',printfilename]);   % save figure in extra output folder


