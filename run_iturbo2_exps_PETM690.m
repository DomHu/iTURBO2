function run_iturbo2_exps_PETM690(datafile , outfilename, load_plot_existing_result, carriers)
%% MATLAB script to reproduce PETM690 observations with iTURBO2 using foram-data in an excel file & 690 observations
%
% datafile:         input excel file with my foram info (age, depth, zbio,abu, iso) to reproduce the observations
% filename:         output filename
% load_plot_existing_result:    load/plot existing results, e.g. the results as shown in the maunscript
%                               measuring n = 4 and n = 12 carriers, specify .mat file further down (~line 20 ff.)
% carriers: how many forams shall be picked and measured
%           Note: existing results only for carriers - 4 and 12!
%
% Example call  run_iturbo2_exps_PETM690('data/3_turbo2input_PETM_690_Exp.xlsx' , 'PETM_results', true)
%               run_iturbo2_exps_PETM690('data/3_turbo2input_PETM_690_Exp.xlsx' , 'PETM_results', false)

% add path to PETM 690 code changes
addpath('./PETM690_code');


% load PETM690 observations
load('data/PETM_690/PETM690.mat')

Exps = 1;

if(load_plot_existing_result)
    if(carriers == 4)
        load('data/PETM_690/zbio_10cm_thermo_0.1_and_1.0_percent_NoDeepMix_v2.mat');     % n=4 carriers
    elseif(carriers == 12)
        load('data/PETM_690/zbio10cm_12carriers_Thermo_0.1_and_1.0_percent_IntermedDeepMix.mat');     % n=12 carriers
    else
        error('Only results for n = 4 and 12 carriers exist.');
    end
else
    % Surface foraminifera (mixed layer)
    data1=xlsread(datafile,'Sheet1','C6:F56');
    %    data1=xlsread(datafile,'Sheet1','P6:S56');
    [biopart11,bioabu1] = iturbo2script_PETM690(data1, carriers, Exps, 'PETM_690_MixedLayer_2.5NewLayer', 'r');
    oriabu1   = data1(:,3);
    oriabu_flipud1 = flipud(oriabu1);
    % % Thermocline foraminifera
    data2=xlsread(datafile,'Sheet1','I6:L56');
    %    data2=xlsread(datafile,'Sheet1','V6:Y56');
    [biopart12,bioabu2] = iturbo2script_PETM690(data2, carriers, Exps, 'PETM_690_Thermocline_2.5NewLayer' , 'b');
    oriabu2   = data2(:,3);
    oriabu_flipud2 = flipud(oriabu2);
    
    age   = data1(:,1);
    age_flipud = flipud(age);
end

set(0,'DefaultAxesFontSize',16)

fig01 = figure('Renderer', 'painters', 'Position', [10 10 900 900]);
% plot 690 single foram observations
subplot(1,3,1)
plot(d13C_690a(:,2),d13C_690a(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
set(gca,'YDir','reverse')
hold on
plot(d13C_690b(:,2),d13C_690b(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690c(:,2),d13C_690c(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690d(:,2),d13C_690d(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690e(:,2),d13C_690e(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690f(:,2),d13C_690f(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690g(:,2),d13C_690g(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690h(:,2),d13C_690h(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690i(:,2),d13C_690i(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690j(:,2),d13C_690j(:,1),'s','MarkerFaceColor',red,'MarkerEdgeColor',red)
plot(d13C_690sa(:,2),d13C_690sa(:,1),'x','MarkerEdgeColor',dblue)
plot(d13C_690sb(:,2),d13C_690sb(:,1),'x','MarkerEdgeColor',dblue)
plot(d13C_690sc(:,2),d13C_690sc(:,1),'x','MarkerEdgeColor',dblue)
plot(d13C_690sd(:,2),d13C_690sd(:,1),'x','MarkerEdgeColor',dblue)
plot(d13C_690se(:,2),d13C_690se(:,1),'x','MarkerEdgeColor',dblue)
plot(d13C_690sf(:,2),d13C_690sf(:,1),'x','MarkerEdgeColor',dblue)
plot(d13C_690sg(:,2),d13C_690sg(:,1),'x','MarkerEdgeColor',dblue)
plot(d13C_690sh(:,2),d13C_690sh(:,1),'x','MarkerEdgeColor',dblue)

plot([-2,5],[170.78,170.78],'-','Color',red)    % example CIE in mixed layer forams
plot([-2,5],[170.7,170.7],'-','Color',dblue)    % example CIE in thermocline forams
ylim([170.4 171])
xlim([-2 5])
ylabel('Meters below seafloor (mbsf) ');
xlabel('\delta^{13}C');
txt = 'ODP 690 obs.';
title(txt)

% plot d13C of iTURBO2 experiment
subplot(1,3,2)
hold on
for i = 1:Exps
    for j=1:carriers
        plot(biopart11(i,:,j),age_flipud, 's','MarkerFaceColor',red,'MarkerEdgeColor',red)
        plot(biopart12(i,:,j),age_flipud, 'x','MarkerEdgeColor',dblue)
    end
end
% create and original d13C
ori_d13C_mixed=oriabu_flipud2;
ori_d13C_mixed(oriabu_flipud2==10000)=3.1;
ori_d13C_mixed(oriabu_flipud2<10000)=0.0;
ori_d13C_thermo=oriabu_flipud2;
ori_d13C_thermo(oriabu_flipud2==10000)=1.6;
ori_d13C_thermo(oriabu_flipud2<10000)=-0.4;
plot(ori_d13C_mixed,age_flipud, '--', 'Color', red,'LineWidth',2)
plot(ori_d13C_thermo,age_flipud, '--','Color', dblue,'LineWidth',2)

plot([-2,5],[170.78,170.78],'-','Color',red)    % example CIE in mixed layer forams
plot([-2,5],[170.7,170.7],'-','Color',dblue)    % example CIE in thermocline forams
set(gca,'YDir','Reverse','XGrid','On','YGrid','On','Box','On', 'XLim',[-2.0,5.0], 'YLim',[170.4,171.0]);%,'YTick',[1.0 1.5 2.0 2.5 3.0])
xlabel('\delta^{13}C');
txt = ['n = ', num2str(carriers) ];
title(txt)

subplot(1,3,3)
hold on
for i = 1:Exps
    plot(100*bioabu1(i,:,1)/10000,age_flipud, 's','MarkerFaceColor',red,'MarkerEdgeColor',red)
    plot(100*bioabu2(i,:,1)/10000,age_flipud, 'x','MarkerEdgeColor',dblue)
end
% original abundance
plot(100*oriabu_flipud2/10000,age_flipud, 'Color', dblue,'LineWidth',2)
plot(100*oriabu_flipud1/10000,age_flipud, 'Color', red,'LineWidth',2)
plot(100*oriabu_flipud2/10000,age_flipud, '--','Color', dblue,'LineWidth',2)
set(gca,'YDir','Reverse','XGrid','On','YGrid','On','Box','On', 'XLim',[0,105.0], 'YLim',[170.4,171.0]);%,'YTick',[1.0 1.5 2.0 2.5 3.0])
xlabel('Rel. pop. size (%)');

if(load_plot_existing_result)
    print(fig01,'-depsc', ['output/Combined_',outfilename , '_as_paper_n=' , num2str(carriers) , 'carriers']);   % save figure in extra output folder
else
    save(['output/mat/',outfilename,'_new_n=' , num2str(carriers) , 'carriers.mat'],'outfilename', 'biopart11','biopart12','age_flipud','bioabu1','bioabu2','oriabu_flipud1','oriabu_flipud2')
    print(fig01, '-depsc', ['output/Combined_',outfilename, '_new_resultsr_n=' , num2str(carriers) , 'carriers']);   % save figure in extra output folder
end

