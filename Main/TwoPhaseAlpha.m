%% Housekeeping
clear;
close all;
clc;

% Select subjects to exclude from group analysis
exSub = {'Sub_151','Sub_152','126'}; %,'Sub_102','Sub_104','Sub_126','Sub_138'};

%% Analysis parameters
% LabView streaming command
params.strcmd = 'str,m,dt,eeg1,eeg2,out,iphs,ps_pulse_count,ps_pulse_id,ps_onphs,ps_offphs,ps_onontime,ps_durbyphs,ps_durbytime,ps_isi';

% Define pre-processing for raw EEG data
params = echtparams('filterImplmentation',0);

CH = {'Fpz','Fz'}; % Electrode channel labels

% Anonymous function to calculate instantaneous frequency from phase
iPhaseToiFreq = @(fs,phi)(fs/(2*pi)*diff(unwrap(phi)));

syncToECHT = 1; % Use instantaneous phase values from ECHT device? [0:no, 1:yes]

%% Load summarized data table
% Load current data table
datafilename = '/Applications/Toolbox/MATLAB/TwoPhaseExp_table.mat';    %% Changed for Windows
load(datafilename);
dataTable = sortrows(dataTable);

%% Load data
pname = genpath('/Applications/Toolbox/SubjectData');  %% Changed for Windows
pname = strsplit(pname,':')'; %% Switched

% Locate Two-Phase ERP datasets
pname = pname(contains(pname,'TwoPhaseERP85dB') & ~contains(pname,'(Low'));

% Exclude subjects?
if(~isempty(exSub))  
    pname = pname(~contains(pname,exSub)); % exclude subject(s) from analysis
end

% Whats my sample?
n = length(pname);

% Get all data
for k = 1:n

    % Get subject info
    tmp = strsplit(pname{k},filesep);
    Sub_ID = tmp{contains(tmp,'Sub_')};
    subIdx = contains(dataTable.SubID,Sub_ID([5:end]));

    % Setup alpha center frequency
    params.alphaCF = dataTable(subIdx,:).IAF(1);

    % Run
    [ERP,params] = batchERP(pname{k},params);

    % Check order of optimal/pessimal phase list in ERP structure against the
    % calculated phases in 'dataTable.' Fix wrap-around error if necessary
    if(dataTable(subIdx,:).OptPhase > dataTable(subIdx,:).PesPhase)
        ERP.setPhases = ERP.setPhases([2 1]);
        fprintf('Switcharoo!\n')
    end

    % Store IAF 
    iaf(k) = params.alphaCF;

    % Store time
    t = ERP.t0;

    % Separate by Phase - all subsequent created matrices have dimensions:
    %
    %           [participant, time, electrode, condition]
    %
    for ph = 1:2
        pIdx{ph} = intersect(ERP.goodTrials,find(round(ERP.echtPhaseLabel)==ERP.setPhases(ph)));

        % Get EROs
        ero(k,:,:,ph) = squeeze(mean(ERP.walpha(:,:,pIdx{ph}),3));

        % Get envelopes
        env(k,:,:,ph) = abs(ero(k,:,:,ph));

        % Get evoked Power
        evo(k,:,:,ph) = env(k,:,:,ph).^2;

        % Get ITPC
        itc(k,:,:,ph) = squeeze(abs(mean(exp(1i*angle(ERP.walpha(:,:,pIdx{ph}))),3)));

        % Z-transform ITPC
        itc(k,:,:,ph) = atan(itc(k,:,:,ph));

        % Get Total power
        tot(k,:,:,ph) = squeeze(mean(abs(ERP.walpha(:,:,pIdx{ph})).^2,3));

        % Get Induced activity
        erp = mean(ERP.erp(:,:,pIdx{ph}),3);
        noerp = ERP.erp(:,:,pIdx{ph}) - erp;

        % Filter alpha band
        [bA,aA] = butter(4,params.fpassA/(params.fs/2),'bandpass');
        indalpha = filtfilt(bA,aA,noerp);
        indalpha = hilbert(indalpha);

        % Get induced power
        ind(k,:,:,ph) = squeeze(mean(abs(indalpha).^2,3));

        clearvars noerp indalpha
    end

end

%% Post Process

% Baseline Period
idxBL = find(t< -50 & t>-250);

% Get time for later fill use
time = [t ; flipud(t)];

% Baseline normalize power measures
evo = evo - mean(evo(:,idxBL,:,:),2);
tot = tot - mean(tot(:,idxBL,:,:),2);
ind = ind - mean(ind(:,idxBL,:,:),2);

% Mean and SEM for all measures - note: SEM are adapted for future use of 
% fill functions
evomu = squeeze(mean(evo)); % Evoked
evose = squeeze(std(evo))./sqrt(n);
evose = cat(1, evomu + evose, flip(evomu - evose,1));

totmu = squeeze(mean(tot)); % Total
totse = squeeze(std(tot))./sqrt(n);
totse = cat(1, totmu + totse, flip(totmu - totse,1));

indmu = squeeze(mean(ind)); % Induced
indse = squeeze(std(ind))./sqrt(n);
indse = cat(1, indmu + indse, flip(indmu - indse,1));

itcmu = squeeze(mean(itc)); % ITPC
itcse = squeeze(std(itc))./sqrt(n);
itcse = cat(1, itcmu + itcse, flip(itcmu - itcse,1));

eromu = squeeze(mean(ero)); % ERO
erose = squeeze(std(ero))./sqrt(n);
erose = cat(1, eromu + erose, flip(eromu - erose,1));

envmu = squeeze(mean(env)); % Envelope
envse = squeeze(std(env))./sqrt(n);
envse = cat(1, envmu + envse, flip(envmu - envse,1));

% Grab the gas for figures
co2 = [0.0000 0.4470 0.7410
       0.8500 0.10 0.0980];

%% Plot EROs

% Trough FPZ
f0 = figure;
plot(t, eromu(:,1,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time, erose(:,1,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
pbaspect([3 1 1])
xlim([min(t) max(t)])
box off
title('Trough FPZ')

% Peak FPZ
f1 = figure;
plot(t, eromu(:,1,2),'color',co2(2,:),'linewidth',2)
hold on
fill(time, erose(:,1,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
pbaspect([3 1 1])
xlim([min(t) max(t)])
box off
title('Peak FPZ')

% Trough FZ
f2 = figure;
plot(t, eromu(:,2,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time, erose(:,2,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
pbaspect([3 1 1])
xlim([min(t) max(t)])
box off
title('Trough FZ')

% Peak FZ
f3 = figure;
plot(t, eromu(:,2,2),'color',co2(2,:),'linewidth',2)
hold on
fill(time, erose(:,2,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
pbaspect([3 1 1])
xlim([min(t) max(t)])
box off
title('Peak FZ')


%% Plot Envelope Comparisons across electrodes

% Statistics
[penv] = npperm(env(:,:,1,:), env(:,:,2,:),1000,"cluster",0,0);
[m,idx] = min(penv);

% Figure Panel
f4 = figure;
for ch = 1:2

    % Raw plots
    ax(3*(ch-1)+1) = subplot(2,3,3*(ch-1)+1);
    hold on;
    patch(ax(3*(ch-1)+1),time,envse(:,1,ch),co2(ch,:),'LineStyle','none','Facealpha',0.1);
    patch(ax(3*(ch-1)+1),time,envse(:,2,ch),co2(ch,:),'LineStyle','none','Facealpha',0.1);
    axis square
    p2a(ch,:) = plot(ax(3*(ch-1)+1),t,squeeze(envmu(:,:,ch)),'LineWidth',2);
    line(ax(3*(ch-1)+1),[0 0],[0 6],'Color','k','LineStyle','--');
    p2a(ch,1).Color = co2(ch,:);
    p2a(ch,2).Color = co2(ch,:);
    ax(3*(ch-1)+1).XGrid = 'on'; ax(3*(ch-1)+1).YGrid = 'on';
    ax(3*(ch-1)+1).XLabel.String = 'Time (ms)';
    ax(3*(ch-1)+1).YLabel.String = 'Envelope Amplitude (µV)';
    ax(3*(ch-1)+1).Title.String = sprintf('%s: Grand Average Alpha Envelopes',CH{ch});
    ax(3*(ch-1)+1).Title.FontSize = 14;
    ax(3*(ch-1)+1).XLim = [-250 500];
    ax(3*(ch-1)+1).YLim = [0 8];
    ax(3*(ch-1)+1).YTick = [0 2 4 6 8];
    ax(3*(ch-1)+1).XTick = [-200 -100 0 100 200 300 400 500];
    title(ax(3*(ch-1)+1).Legend,'Onset Phase');
    pbaspect(ax(3*(ch-1)+1), [3 1 1])

    % Mean differences
    ax(3*(ch-1)+2) = subplot(2,3,3*(ch-1)+2);
    hold on;
    axis square
    mY(:,ch) = envmu(:,2,ch) - envmu(:,1,ch);
    semY(:,ch) = squeeze(std(env(:,:,2,ch) - env(:,:,1,ch)))/sqrt(size(env,1));
    patch(ax(3*(ch-1)+2),time,cat(1,mY(:,ch)+semY(:,ch),flipud(mY(:,ch)-semY(:,ch))),...
        'k','LineStyle','none','FaceAlpha',0.1);
    p2b(ch) = plot(ax(3*(ch-1)+2),t,mY(:,ch),'Color','k','LineWidth',2);
    text(ax(3*(ch-1)+2),t(idx(ch)),-.9,sprintf('* p = %f',m(ch)));
    line(ax(3*(ch-1)+2),xlim,[0 0],'Color','k','LineStyle','--');
    line(ax(3*(ch-1)+2),[0 0],[1 1],'Color','k','LineStyle','--');
    ax(3*(ch-1)+2).XGrid = 'on'; ax(3*(ch-1)+2).YGrid = 'on';
    ax(3*(ch-1)+2).Title.String = sprintf('Ave Within-Subject Diff: %s (n=%d)',CH{ch},size(env,1));
    ax(3*(ch-1)+2).Title.FontSize = 14;
    ax(3*(ch-1)+2).XLabel.String = 'Time (ms)';
    ax(3*(ch-1)+2).YLabel.String = 'Voltage (µV)';
    ax(3*(ch-1)+2).YLim = [-1 1];
    ax(3*(ch-1)+2).XLim = [-250 500];
    ax(3*(ch-1)+2).XTick = [-200 -100 0 100 200 300 400 500];
    pbaspect(ax(3*(ch-1)+2), [3 1 1])

    % Signifigance
    ax(3*(ch-1)+3) = subplot(2,3,3*(ch-1)+3);
    hold on;
    axis square
    plot(ax(3*(ch-1)+3), t,penv(:,ch),'k-')
    fill(ax(3*(ch-1)+3),time, [penv(:,ch); ones(1,length(penv(:,ch)))'],'k')
    ax(3*(ch-1)+3).Title.String = sprintf('Cluster Permutation Results: %s (n=%d)',CH{ch},size(env,1));
    ax(3*(ch-1)+3).YLim = [0 0.05] ;
    ax(3*(ch-1)+3).XLim = [-250 500] ;
    ax(3*(ch-1)+3).XTick = [-200 -100 0 100 200 300 400 500];
    ax(3*(ch-1)+3).YLabel.String = 'P-Value';
    ax(3*(ch-1)+3).XLabel.String = 'Time (ms)';
    ax(3*(ch-1)+3).YDir = 'reverse';
    pbaspect(ax(3*(ch-1)+3), [3 1 1])
    
end


%% Plot Envelope Comparisons across conditions

% Statistics
[penv] = npperm(env(:,:,:,1), env(:,:,:,2),1000,"cluster",0,0);
[pks, loc] = findpeaks(-penv(:,1),t);
[m,idx] = min(penv);

% Figure Panel
f5 = figure;
for ch = 1:2

    % Raw plots
    ax1(3*(ch-1)+1) = subplot(2,3,3*(ch-1)+1);
    hold on;
    patch(ax1(3*(ch-1)+1),time,envse(:,ch,1),co2(1,:),'LineStyle','none','Facealpha',0.1);
    patch(ax1(3*(ch-1)+1),time,envse(:,ch,2),co2(2,:),'LineStyle','none','Facealpha',0.1);
    axis square
    p2a(ch,:) = plot(ax1(3*(ch-1)+1),t,squeeze(envmu(:,ch,:)),'LineWidth',2);
    line(ax1(3*(ch-1)+1),[0 0],[0 6],'Color','k','LineStyle','--');
    p2a(ch,1).Color = co2(1,:);
    p2a(ch,2).Color = co2(2,:);
    ax1(3*(ch-1)+1).XGrid = 'on'; ax1(3*(ch-1)+1).YGrid = 'on';
    ax1(3*(ch-1)+1).XLabel.String = 'Time (ms)';
    ax1(3*(ch-1)+1).YLabel.String = 'Envelope Amplitude (µV)';
    ax1(3*(ch-1)+1).Title.String = sprintf('%s: Grand Average Alpha Envelopes',CH{ch});
    ax1(3*(ch-1)+1).Title.FontSize = 14;
    ax1(3*(ch-1)+1).XLim = [-250 500];
    ax1(3*(ch-1)+1).YLim = [0 8];
    ax1(3*(ch-1)+1).YTick = [0 2 4 6 8];
    ax1(3*(ch-1)+1).XTick = [-200 -100 0 100 200 300 400 500];
    title(ax1(3*(ch-1)+1).Legend,'Onset Phase');
    pbaspect(ax1(3*(ch-1)+1), [3 1 1])

    % Mean differences
    ax1(3*(ch-1)+2) = subplot(2,3,3*(ch-1)+2);
    hold on;
    axis square
    mY(:,ch) = envmu(:,ch,2) - envmu(:,ch,1);
    semY(:,ch) = squeeze(std(env(:,:,ch,2) - env(:,:,ch,1)))/sqrt(size(env,1));
    patch(ax1(3*(ch-1)+2),time,cat(1,mY(:,ch)+semY(:,ch),flipud(mY(:,ch)-semY(:,ch))),...
        'k','LineStyle','none','FaceAlpha',0.1);
    p2b(ch) = plot(ax1(3*(ch-1)+2),t,mY(:,ch),'Color','k','LineWidth',2);
    text(ax1(3*(ch-1)+2),t(idx(ch)),-.9,sprintf('* p = %f',m(ch)));
    line(ax1(3*(ch-1)+2),xlim,[0 0],'Color','k','LineStyle','--');
    line(ax1(3*(ch-1)+2),[0 0],[-2 2],'Color','k','LineStyle','--');
    ax1(3*(ch-1)+2).XGrid = 'on'; ax1(3*(ch-1)+2).YGrid = 'on';
    ax1(3*(ch-1)+2).Title.String = sprintf('Ave Within-Subject Diff: %s (n=%d)',CH{ch},size(env,1));
    ax1(3*(ch-1)+2).Title.FontSize = 14;
    ax1(3*(ch-1)+2).XLabel.String = 'Time (ms)';
    ax1(3*(ch-1)+2).YLabel.String = 'Voltage (µV)';
    ax1(3*(ch-1)+2).YLim = [-2 2];
    ax1(3*(ch-1)+2).XLim = [-250 500];
    ax1(3*(ch-1)+2).XTick = [-200 -100 0 100 200 300 400 500];
    pbaspect(ax1(3*(ch-1)+2), [3 1 1])

    % Signifigance
    ax1(3*(ch-1)+3) = subplot(2,3,3*(ch-1)+3);
    hold on;
    axis square
    plot(ax1(3*(ch-1)+3), t,penv(:,ch),'k-')
    fill(ax1(3*(ch-1)+3),time, [penv(:,ch); ones(1,length(penv(:,ch)))'],'k')
    ax1(3*(ch-1)+3).Title.String = sprintf('Cluster Permutation Results: %s (n=%d)',CH{ch},size(env,1));
    ax1(3*(ch-1)+3).YLim = [0 0.05] ;
    ax1(3*(ch-1)+3).XLim = [-250 500] ;
    ax1(3*(ch-1)+3).XTick = [-200 -100 0 100 200 300 400 500];
    ax1(3*(ch-1)+3).YLabel.String = 'P-Value';
    ax1(3*(ch-1)+3).XLabel.String = 'Time (ms)';
    ax1(3*(ch-1)+3).YDir = 'reverse';
    pbaspect(ax1(3*(ch-1)+3), [3 1 1])
    
end


%% Plot power & ITPC analyses

% Statistics
[ptot] = npperm(tot(:,:,:,1), tot(:,:,:,2),1000,"cluster",0,0);
[pevo] = npperm(evo(:,:,:,1), evo(:,:,:,2),1000,"cluster",0,0);
[pind] = npperm(ind(:,:,:,1), ind(:,:,:,2),1000,"cluster",0,0);
[pitc] = npperm(itc(:,:,:,1), itc(:,:,:,2),1000,"cluster",0,0);

% ITPC Figures
f6 = figure;
subplot(2,1,1)
plot(t, itcmu(:,1,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time,itcse(:,1,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, itcmu(:,1,2),'color',co2(2,:),'linewidth',2)
fill(time,itcse(:,1,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([0 max(t)])
xlabel('Time (ms)')
ylabel('ITPC')
axis square
box off
title('ITPC FPZ')

subplot(2,1,2)
plot(t,pitc(:,ch),'k-')
hold on
fill(time, [pitc(:,ch); ones(1,length(pitc(:,ch)))'],'k')
xlim([0 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
axis square
box off

f7 = figure;
subplot(2,1,1)
plot(t, itcmu(:,2,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time,itcse(:,2,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, itcmu(:,2,2),'color',co2(2,:),'linewidth',2)
fill(time,itcse(:,2,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([0 max(t)])
xlabel('Time (ms)')
ylabel('ITPC')
axis square
box off
title('ITPC FZ')

subplot(2,1,2)
plot(t,pitc(:,ch),'k-')
hold on
fill(time, [pitc(:,ch); ones(1,length(pitc(:,ch)))'],'k')
xlim([0 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
axis square
box off

% Evoked Power Figures
f8 = figure;
subplot(2,1,1)
plot(t, evomu(:,1,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time,evose(:,1,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, evomu(:,1,2),'color',co2(2,:),'linewidth',2)
fill(time,evose(:,1,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([0 max(t)])
xlabel('Time (ms)')
ylabel('Evoked Power')
axis square
box off
title('Evoked FPZ')

subplot(2,1,2)
plot(t,pevo(:,ch),'k-')
hold on
fill(time, [pevo(:,ch); ones(1,length(pevo(:,ch)))'],'k')
xlim([0 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
axis square
box off

f9 = figure;
subplot(2,1,1)
plot(t, evomu(:,2,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time,evose(:,2,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, evomu(:,2,2),'color',co2(2,:),'linewidth',2)
fill(time,evose(:,2,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([0 max(t)])
xlabel('Time (ms)')
ylabel('Evoked Power')
axis square
box off
title('Evoked FZ')

subplot(2,1,2)
plot(t,pevo(:,ch),'k-')
hold on
fill(time, [pevo(:,ch); ones(1,length(pevo(:,ch)))'],'k')
xlim([0 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
axis square
box off

% Total Power Figures
f10 = figure;
subplot(2,1,1)
plot(t, totmu(:,1,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time,totse(:,1,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, totmu(:,1,2),'color',co2(2,:),'linewidth',2)
fill(time,totse(:,1,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([0 max(t)])
xlabel('Time (ms)')
ylabel('Total Power')
axis square
box off
title('Total FPZ')

subplot(2,1,2)
plot(t,ptot(:,ch),'k-')
hold on
fill(time, [ptot(:,ch); ones(1,length(ptot(:,ch)))'],'k')
xlim([0 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
axis square
box off

f11 = figure;
subplot(2,1,1)
plot(t, totmu(:,2,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time,totse(:,2,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, totmu(:,2,2),'color',co2(2,:),'linewidth',2)
fill(time,totse(:,2,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([0 max(t)])
xlabel('Time (ms)')
ylabel('Total Power')
axis square
box off
title('Total FZ')

subplot(2,1,2)
plot(t,ptot(:,ch),'k-')
hold on
fill(time, [ptot(:,ch); ones(1,length(ptot(:,ch)))'],'k')
xlim([0 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
axis square
box off

% Induced Power Figures
f12 = figure;
subplot(2,1,1)
plot(t, indmu(:,1,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time,indse(:,1,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, indmu(:,1,2),'color',co2(2,:),'linewidth',2)
fill(time,indse(:,1,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([0 max(t)])
xlabel('Time (ms)')
ylabel('Induced Power')
axis square
box off
title('Induced FPZ')

subplot(2,1,2)
plot(t,pind(:,ch),'k-')
hold on
fill(time, [pind(:,ch); ones(1,length(pind(:,ch)))'],'k')
xlim([0 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
axis square
box off

f13 = figure;
subplot(2,1,1)
plot(t, indmu(:,2,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time,indse(:,2,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, indmu(:,2,2),'color',co2(2,:),'linewidth',2)
fill(time,indse(:,2,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([0 max(t)])
xlabel('Time (ms)')
ylabel('Induced Power')
axis square
box off
title('Induced FZ')

subplot(2,1,2)
plot(t,pind(:,ch),'k-')
hold on
fill(time, [pind(:,ch); ones(1,length(pind(:,ch)))'],'k')
xlim([0 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
axis square
box off

%% Save Figures

svpath = '/Applications/Toolbox/MATLAB/Figures/Alpha HT/';

%%%%%%%%%%% PNG
% Save EROs
print(f0,fullfile(svpath,'Trough ERO FPZ'),'-dpng');
print(f1,fullfile(svpath,'Peak ERO FPZ'),'-dpng');
print(f2,fullfile(svpath,'Trough ERO FZ'),'-dpng');
print(f3,fullfile(svpath,'Peak ERO FZ'),'-dpng');

% Save Envelopes
print(f4,fullfile(svpath,'Envelopes by Electrode'),'-dpng');
print(f5, fullfile(svpath,'Envelopes by Condition'),'-dpng');

% Save 
print(f6, fullfile(svpath, 'ITPC FPZ'),'-dpng');
print(f7, fullfile(svpath, 'ITPC FZ'),'-dpng');
print(f8, fullfile(svpath, 'Evoked Power FPZ'),'-dpng');
print(f9, fullfile(svpath, 'Evoked Power FZ'),'-dpng');
print(f10, fullfile(svpath,'Total Power FPZ'),'-dpng');
print(f11, fullfile(svpath,'Total Power FZ'),'-dpng');
print(f12, fullfile(svpath,'Induced Power FPZ'),'-dpng');
print(f13, fullfile(svpath,'Induced Power FZ'),'-dpng');

%%%%%%%%%%% SVG
% Save EROs
print(f0,fullfile(svpath,'Trough ERO FPZ'),'-dsvg');
print(f1,fullfile(svpath,'Peak ERO FPZ'),'-dsvg');
print(f2,fullfile(svpath,'Trough ERO FZ'),'-dsvg');
print(f3,fullfile(svpath,'Peak ERO FZ'),'-dsvg');

% Save Envelopes
print(f4,fullfile(svpath,'Envelopes by Electrode'),'-dsvg');
print(f5, fullfile(svpath,'Envelopes by Condition'),'-dsvg');

% Save 
print(f6, fullfile(svpath, 'ITPC FPZ'),'-dsvg');
print(f7, fullfile(svpath, 'ITPC FZ'),'-dsvg');
print(f8, fullfile(svpath, 'Evoked Power FPZ'),'-dsvg');
print(f9, fullfile(svpath, 'Evoked Power FZ'),'-dsvg');
print(f10, fullfile(svpath,'Total Power FPZ'),'-dsvg');
print(f11, fullfile(svpath,'Total Power FZ'),'-dsvg');
print(f12, fullfile(svpath,'Induced Power FPZ'),'-dsvg');
print(f13, fullfile(svpath,'Induced Power FZ'),'-dsvg');


