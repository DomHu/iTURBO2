function iturbo2script_4zbio_ASH(data, nocarriers, noexps, expname, TM, Observation)
%% run iTURBO2 experiments and plot results against ash observations
% use 4 different bioturbation depth depending on the observed ash-profiles

% data = matrix of required data (age, mxl, abu, iso) - explanation see below
% nocarriers = number of carriers to be measured
% noexps = number of experiments for SA
% expname = experiment name
% Observation:  1: cores with 0.5       cm kyr-1 sed-rate
%               2: cores with 2.0 - 2.5 cm kyr-1 sed-rate
% age = age of sediment layer down core
% mxl = series of mixed layer thicknesses (zbio) down core
% abu = series of abundances of carrier type 1 down core
% iso = original isotope signature of both carrier  types 1 and 2


age   = data(:,1);
abu   = data(:,3);
iso   = data(:,4);
lngth = length(data(:,1));

% distinguish between ash observation experiment
if(Observation ==1 )     % 0.5 cm kyr-1 rate
    mxl1   = data(:,2).*11;
    mxl2 = data(:,2).*13;
    mxl3 = data(:,2).*15;
    mxl4 = data(:,2).*7;
    Change_rel_depth = 53;  % to align max ash concentration observed with simulated using intermediate zbio
    txt = '0.5 cm kyr^{-1}';
    if(TM == true)
        Change_rel_depth = 38;  % to align max ash concentration observed with simulated using intermediate zbio
    end
    
    % load ash observations vs relative depth:
    data_01=xlsread('data/iTURBO2_input_ash_experiment.xlsx','ash_data','P30:Q43');
    data_02=xlsread('data/iTURBO2_input_ash_experiment.xlsx','ash_data','S30:T43');

else                    % 2.0 - 2.5 cm kyr-1 rate
    mxl1   = data(:,2).*5;
    mxl2 = data(:,2).*7;
    mxl3 = data(:,2).*9;
    mxl4 = data(:,2).*13;
    Change_rel_depth = 47;  % to align max ash concentration observed with simulated using intermediate zbio
    txt = '2.0 - 2.5 cm kyr^{-1}';
    if(TM == true)
        Change_rel_depth = 40;  % to align max ash concentration observed with simulated using intermediate zbio
    end    
    data_01=xlsread('data/iTURBO2_input_ash_experiment.xlsx','ash_data','V30:W45');
    data_02=xlsread('data/iTURBO2_input_ash_experiment.xlsx','ash_data','Y30:Z47');
end

numb  = nocarriers;     % number of carriers to be measured
exps = noexps;          % number of different experiments

%%
for i = 1:exps
    [oriabu(i,:,:),bioabu(i,:,:),oriiso(i,:,:),bioiso(i,:,:)] = iturbo2_plus_TM(abu,iso,mxl1,numb, TM);
    [oriabu2(i,:,:),bioabu2(i,:,:),oriiso2(i,:,:),bioiso2(i,:,:)] = iturbo2_plus_TM(abu,iso,mxl2,numb, TM);
    [oriabu3(i,:,:),bioabu3(i,:,:),oriiso3(i,:,:),bioiso3(i,:,:)] = iturbo2_plus_TM(abu,iso,mxl3,numb, TM);
    [oriabu4(i,:,:),bioabu4(i,:,:),oriiso4(i,:,:),bioiso4(i,:,:)] = iturbo2_plus_TM(abu,iso,mxl4,numb, TM);
end
% normalize bioabu:
bioabu_norm = bioabu./nocarriers;
bioabu2_norm = bioabu2./nocarriers;
bioabu3_norm = bioabu3./nocarriers;
bioabu4_norm = bioabu4./nocarriers;

% variable for mean results of mxl1
mean_bioabu1_mxl1 = zeros(1,lngth);
mean_bioabu2_mxl1 = zeros(1,lngth);
mean_bioiso1_mxl1 = zeros(1,lngth);
mean_bioiso2_mxl1 = zeros(1,lngth);
% variable for mean results of mxl2
mean_bioabu1_mxl2 = zeros(1,lngth);
mean_bioabu2_mxl2 = zeros(1,lngth);
mean_bioiso1_mxl2 = zeros(1,lngth);
mean_bioiso2_mxl2 = zeros(1,lngth);
% variable for mean results of mxl3
mean_bioabu1_mxl3 = zeros(1,lngth);
mean_bioabu2_mxl3 = zeros(1,lngth);
mean_bioiso1_mxl3 = zeros(1,lngth);
mean_bioiso2_mxl3 = zeros(1,lngth);
% variable for mean results of mxl4
mean_bioabu1_mxl4 = zeros(1,lngth);
mean_bioabu2_mxl4 = zeros(1,lngth);
mean_bioiso1_mxl4 = zeros(1,lngth);
mean_bioiso2_mxl4 = zeros(1,lngth);


%%
mxltext = num2str(mean(mxl1));
mxltext2 = num2str(mean(mxl2));
mxltext3 = num2str(mean(mxl3));
mxltext4 = num2str(mean(mxl4));
numbtxt = num2str(numb,2);
expstxt = num2str(exps,2);
abutxt = num2str(max(abu),2);


rel_depth = (1:lngth)-Change_rel_depth;


%%  Plot normalized abundance only and just for species 1 (no isotope change here)
set(0,'DefaultAxesFontSize',16)

figure
hold on
for i = 1:exps
    % mxl1
    mean_bioabu1_mxl1 = mean_bioabu1_mxl1+bioabu_norm(i,:,1);
    % mxl2
    mean_bioabu1_mxl2 = mean_bioabu1_mxl2+bioabu2_norm(i,:,1);
    % mxl3
    mean_bioabu1_mxl3 = mean_bioabu1_mxl3+bioabu3_norm(i,:,1);
    % mxl4
    mean_bioabu1_mxl4 = mean_bioabu1_mxl4+bioabu4_norm(i,:,1);
end
mean_bioabu1_mxl1 = mean_bioabu1_mxl1/exps;
mean_bioabu1_mxl2 = mean_bioabu1_mxl2/exps;
mean_bioabu1_mxl3 = mean_bioabu1_mxl3/exps;
mean_bioabu1_mxl4 = mean_bioabu1_mxl4/exps;


%% calculate residual sum of squares
% save relative bioturbated-abundance of respective depths in the
% observations in columns 3, 4, 5, 6
for i=1:size(data_01,1)
    index = find(rel_depth==data_01(i,1));
    if(isempty(index)) % depth not in our bioturbated core --> abu = 0
        data_01(i,3) = 0;
        data_01(i,4) = 0;
        data_01(i,5) = 0;
        data_01(i,6) = 0;
    else
    data_01(i,3) = mean_bioabu1_mxl1(index);
    data_01(i,4) = mean_bioabu1_mxl2(index);
    data_01(i,5) = mean_bioabu1_mxl3(index);
    data_01(i,6) = mean_bioabu1_mxl4(index);
    end
end

for i=1:size(data_02,1)
    index = find(rel_depth==data_02(i,1));
    if(isempty(index)) % depth not in our bioturbated core --> abu = 0
        data_02(i,3) = 0;
        data_02(i,4) = 0;
        data_02(i,5) = 0;
        data_02(i,6) = 0;
    else
    data_02(i,3) = mean_bioabu1_mxl1(index);
    data_02(i,4) = mean_bioabu1_mxl2(index);
    data_02(i,5) = mean_bioabu1_mxl3(index);
    data_02(i,6) = mean_bioabu1_mxl4(index);
    end
end

% mxl1
RSS_mxl1 = sum((data_01(:,2)-data_01(:,3)).^2) + sum((data_02(:,2)-data_02(:,3)).^2);
% mxl2
RSS_mxl2 = sum((data_01(:,2)-data_01(:,4)).^2) + sum((data_02(:,2)-data_02(:,4)).^2);
% mxl3
RSS_mxl3 = sum((data_01(:,2)-data_01(:,5)).^2) + sum((data_02(:,2)-data_02(:,5)).^2);
% mxl4
RSS_mxl4 = sum((data_01(:,2)-data_01(:,6)).^2) + sum((data_02(:,2)-data_02(:,6)).^2);

RSS_mxl1_text = num2str(round(1000*RSS_mxl1)/1000);
RSS_mxl2_text = num2str(round(1000*RSS_mxl2)/1000);
RSS_mxl3_text = num2str(round(1000*RSS_mxl3)/1000);
RSS_mxl4_text = num2str(round(1000*RSS_mxl4)/1000);

% plot relative depths:
plot(rel_depth,oriabu(1,:,1)./nocarriers,'--k','Linewidth',2.0) % plot one of the original abu

% mxl1
plot(rel_depth,mean_bioabu1_mxl1, '-r','Linewidth',2.0)
% mxl2
plot(rel_depth,mean_bioabu1_mxl2, '-g','Linewidth',2.0)
% mxl3
plot(rel_depth,mean_bioabu1_mxl3, '-b','Linewidth',2.0)
% mxl4
plot(rel_depth,mean_bioabu1_mxl4, '-','Color',[0.9290 0.6940 0.1250],'Linewidth',2.0)

% plot observations
plot(data_01(:,1),data_01(:,2),'ko','MarkerFaceColor','k')
plot(data_02(:,1),data_02(:,2),'k^','MarkerFaceColor','k')
hleg=legend(['z_{bio}= 0 cm'], ['z_{bio}= ' , mxltext, ' cm, RSS=', RSS_mxl1_text],['z_{bio}= ' , mxltext2,' cm, RSS=', RSS_mxl2_text],['z_{bio}= ', mxltext3,' cm, RSS=', RSS_mxl3_text],['z_{bio}= ', mxltext4,' cm, RSS=', RSS_mxl4_text]);

set(gca,'XGrid','On','YGrid','On', 'YLim',[0, 0.2],'YTick',[0.0 0.05 0.1 0.15 0.2])
set(hleg,'FontSize',8);
set(hleg,'Location','NorthEast');
xlim([-30 20])
xlabel('Core depth (cm) ');
ylabel('Normalized ash concentration');
text(0.04, 0.90, txt, 'FontSize', 14, 'Units', 'normalized');
hold off
printfilename = ['4zbio_',expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_',num2str(TM),'TM'];

% check if output directory exists -- if not create it:
if ~(exist('output/mat','dir') == 7), mkdir('output/mat'); end

    
save(['output/mat/',printfilename,'.mat'],'printfilename', 'lngth','bioiso','bioiso2','bioiso3','bioiso4', 'oriiso', 'mean_bioiso1_mxl1', 'mean_bioiso1_mxl2', 'mean_bioiso1_mxl3', 'mean_bioiso1_mxl4','expname', 'exps', 'bioabu','bioabu2','bioabu3','bioabu4', 'oriabu', 'oriabu2', 'oriabu3', 'oriabu4', 'mean_bioabu1_mxl1', 'mean_bioabu1_mxl2', 'mean_bioabu1_mxl3', 'mean_bioabu1_mxl4')
print('-depsc', ['output/',printfilename]);   % save figure in extra output folder

