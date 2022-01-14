function iturbo2script_multipleSims(data, nocarriers, noexps, expname, TM)
%% MATLAB script to run TURBO2 experiment with input data
%% ALSO plots the results - either just the isotope signal
% or if abundance change is wanted set plot_abu_iso = true;
% or if just the mean result of all experiments is wanted set plot_just_mean = true;

% data = matrix of required data (age, mxl, abu, iso) - explanation see below
% nocarriers = number of carriers to be measured
% noexps = number of experiments for SA
% expname = experiment name

% age = age of sediment layer down core
% mxl = series of mixed layer thicknesses (zbio) down core
% abu = series of abundances of carrier type 1 down core
% iso = original isotope signature of both carrier  types 1 and 2

%%
plot_results = true;       % plot results or just input signal?
plot_iso_in_2figs = true;
plot_abu = true;            % plot the abundance figures?

age   = data(:,1);
mxl   = data(:,2);
abu   = data(:,3);
iso   = data(:,4);
lngth = length(data(:,1));

numb  = nocarriers;     %50;     % number of carriers to be measured
exps = noexps;     % number of different experiments

%%
for i = 1:exps
    i
    [oriabu(i,:,:),bioabu(i,:,:),oriiso(i,:,:),bioiso(i,:,:)] = iturbo2_plus_TM(abu,iso,mxl,numb, TM);
end
%%
mean_bioabu1 = zeros(1,lngth);
mean_bioabu2 = zeros(1,lngth);
mean_bioiso1 = zeros(1,lngth);
mean_bioiso2 = zeros(1,lngth);

%%
mxltext = num2str(mean(mxl));
numbtxt = num2str(numb,2);
expstxt = num2str(exps,2);
abutxt = num2str(max(abu),2);

set(0,'DefaultAxesFontSize',16)

mean_bioiso1 = zeros(1,lngth);
mean_bioiso2 = zeros(1,lngth);

%%  Plot isotope plots only and in one figure
fig1 = figure;
hold on
for i = 1:exps
    if(plot_results)
        plot(1:lngth,bioiso(i,:,1), 'Color', [0.8 0.8 1.0],'Linewidth',1.5)
        plot(1:lngth,bioiso(i,:,2), 'Color', [1.0 0.8 0.8],'Linewidth',1.5)
    end
    mean_bioiso1 = mean_bioiso1+bioiso(i,:,1);
    mean_bioiso2 = mean_bioiso2+bioiso(i,:,2);
end
plot(1:lngth,oriiso(1,:,1),'k','Linewidth',2.0) % plot one of the original iso
mean_bioiso1 = mean_bioiso1/exps;
mean_bioiso2 = mean_bioiso2/exps;
if(plot_results)
    plot(1:lngth,mean_bioiso1, '-b','Linewidth',2.0)
    plot(1:lngth,mean_bioiso2, '--r','Linewidth',2.0)
end
set(gca,'YDir','Reverse','XGrid','On','YGrid','On','Box','On', 'XLim',[0,200], 'YLim',[1.0,3.0],'YTick',[1.0 1.5 2.0 2.5 3.0])
xlabel('Core depth (cm) ');
ylabel('\delta^{18}O');
titletxt = ['Isotopes of Carriers 1+2, ',mxltext,...
    ' cm Mixed Layer, ',numbtxt,' Carriers'];
title('Isotope signals')


if(plot_results)
    printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_ISO'];
else
    printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_ISO_JUST_INPUT'];
end
print(fig1, '-depsc', ['output/',printfilename]);   % save figure in extra output folder


%% plot isotope results in 2 figures
if(plot_iso_in_2figs)
    plot_results = true;
    % Species 1
    fig2 = figure; 
    hold on
    for i = 1:exps
        if(plot_results)
            plot(1:lngth,bioiso(i,:,1), 'Color', [0.8 0.8 1.0],'Linewidth',1.5)
        end
    end
    plot(1:lngth,oriiso(1,:,1),'k','Linewidth',2.0) % plot one of the original iso
    if(plot_results)
        plot(1:lngth,mean_bioiso1, '-b','Linewidth',2.0)
    end
    set(gca,'YDir','Reverse','XGrid','On','YGrid','On','Box','On', 'XLim',[0,200], 'YLim',[1.0,3.0],'YTick',[1.0 1.5 2.0 2.5 3.0])
    xlabel('Core depth (cm) ');
    ylabel('\delta^{18}O');
    titletxt = ['Isotopes of Carriers 1, ',mxltext,...
        ' cm Mixed Layer, ',numbtxt,' Carriers'];
    title('Isotope signals cold species')
    
    if(plot_results)
        printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_ISO_Species1'];
    else
        printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_ISO_Species1_JUST_INPUT'];
    end
    print(fig2, '-depsc', ['output/',printfilename]);   % save figure in extra output folder
    
    % Species 2
    fig3 = figure;
    hold on
    for i = 1:exps
        if(plot_results)
            plot(1:lngth,bioiso(i,:,2), 'Color', [1.0 0.8 0.8],'Linewidth',1.5)
        end
    end
    plot(1:lngth,oriiso(1,:,1),'k','Linewidth',2.0) % plot one of the original iso
    if(plot_results)
        plot(1:lngth,mean_bioiso2, '-r','Linewidth',2.0)
    end
    set(gca,'YDir','Reverse','XGrid','On','YGrid','On','Box','On', 'XLim',[0,200], 'YLim',[1.0,3.0],'YTick',[1.0 1.5 2.0 2.5 3.0])
    xlabel('Core depth (cm) ');
    ylabel('\delta^{18}O');
    titletxt = ['Isotopes of Carriers 2, ',mxltext,...
        ' cm Mixed Layer, ',numbtxt,' Carriers'];
    title('Isotope signals warm species')
    
    if(plot_results)
        printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_ISO_Species2'];
    else
        printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_ISO_Species2_JUST_INPUT'];
    end
    print(fig3, '-depsc', ['output/',printfilename]);   % save figure in extra output folder
    
end

%% Plot abundance plots
if(plot_abu)
    
    fig4 = figure;
    hold on
    for i = 1:exps
        if(plot_results)
            plot(1:lngth,bioabu(i,:,1), 'Color', [0.8 0.8 1.0],'Linewidth',2.0)
        end
        mean_bioabu1 = mean_bioabu1+bioabu(i,:,1);
    end
    plot(1:lngth,oriabu(1,:,1),'k','Linewidth',3.0) % plot one of the original abu
    mean_bioabu1 = mean_bioabu1/exps;
    if(plot_results)
        plot(1:lngth,mean_bioabu1, '-b','Linewidth',3.0)
    end
    plot(1:lngth,numb*ones(lngth),'g')
    set(gca,'XGrid','On','YGrid','On','Box','On', 'XLim',[0,200], 'YLim',[0,1000],'YTick',(0:200:1000))
    xlabel('Core depth (cm) ');
    ylabel('Number of Particles');
    % legend('Original abundance','Bioturbated abundance Species 1')
    title('Abundance of cold species')
    
    if(plot_results)
        printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_SP1'];
    else
        printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_SP1_JUST_INPUT'];
    end
    print(fig4, '-depsc', ['output/', printfilename]);
    
    
    fig5 = figure;
    hold on
    for i = 1:exps
        if(plot_results)
            plot(1:lngth,bioabu(i,:,2), 'Color', [1.0 0.8 0.8],'Linewidth',2.0)
        end
        mean_bioabu2 = mean_bioabu2+bioabu(i,:,2);
    end
    plot(1:lngth,oriabu(1,:,2),'k','Linewidth',3.0) % plot one of the original abu
    mean_bioabu2 = mean_bioabu2/exps;
    if(plot_results)
        plot(1:lngth,mean_bioabu2, '-r','Linewidth',3.0)
    end
    plot(1:lngth,numb*ones(lngth),'g')
    set(gca,'XGrid','On','YGrid','On','Box','On', 'XLim',[0,200], 'YLim',[0,1000],'YTick',(0:200:1000))
    
    xlabel('Core depth (cm) ');
    ylabel('Number of Particles');
    % legend('Original abundance','Bioturbated abundance Species 2')
    title('Abundance of warm species')
    
    if(plot_results)
        printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_SP2'];
    else
        printfilename = [expname,'_',abutxt,'abu_',numbtxt,'carriers_',expstxt,'Exps_SP2_JUST_INPUT'];
    end
    print(fig5, '-depsc', ['output/', printfilename]);
    
    % save variables necessary to plot in variable
    save(['output/mat/',printfilename,'.mat'],'bioiso','oriiso','mean_bioiso1','mean_bioiso2','lngth','exps','bioabu', 'oriabu', 'mean_bioabu1', 'mean_bioabu2')
    
    
end

