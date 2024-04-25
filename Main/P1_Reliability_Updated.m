%% Housekeeping
clear;
close all;
clc;

%% Analysis parameters

% Define pre-processing for raw EEG data
params = echtparams('zeroPhase',1,'filterImplementation',1,'fslide',0,'tWin',[-0.250 0.500]);

CH = {'Fpz','Fz'}; % Electrode channel labels

% Select subjects to exclude from group analysis
exSub = {'Sub_151','Sub_152','126'}; %,'Sub_102','Sub_104','Sub_126','Sub_138'};

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
pname1 = pname(contains(pname,'ISI750ms85dB') & ~contains(pname,'(Low'));
pname2 = pname(contains(pname,'TwoPhaseERP85dB') & ~contains(pname,'(Low'));

% Exclude subjects?
if(~isempty(exSub))  
    pname1 = pname1(~contains(pname1,exSub)); % exclude subject(s) from analysis
    pname2 = pname2(~contains(pname2,exSub)); % exclude subject(s) from analysis
end

% Whats my sample?
n = length(pname1);

% Create figure directory
if ~isdir('./P1 Analyses')
    mkdir('./P1 Analyses');
end
if ~isdir('./ParticipantEROs')
    mkdir('./ParticipantEROs');
end

% Grab nice colors
co2 = [0.0000 0.4470 0.7410
       0.8500 0.10 0.0980];


% Get all data
for k = 1:n

    % Get subject info
    tmp = strsplit(pname1{k},filesep);
    Sub_ID = tmp{contains(tmp,'Sub_')};
    subIdx = contains(dataTable.SubID,Sub_ID([5:end]));

    % Setup alpha center frequency
    params.alphaCF = dataTable(subIdx,:).IAF(1);
    iaf(k) = params.alphaCF;

    % Other parameters
    params.fslide = 0;
    params.useZeroPhaseFiltering = 0;
    params.filterImplementation = 1;

    % Run random
    [ERP,params] = batchERP(pname1{k},params);

    % Run 2Phase
    [ERP2,params] = batchERP(pname2{k},params);

    % Store time
    t = ERP.t0;

    % Separate out recordings sessions
    erp = ERP.erp;

    % Get evoked responses
    erp = squeeze(mean(erp,3));

    % Check order of optimal/pessimal phase list in ERP structure against the
    % calculated phases in 'dataTable.' Fix wrap-around error if necessary
    if(dataTable(subIdx,:).OptPhase > dataTable(subIdx,:).PesPhase)
        ERP2.setPhases = ERP2.setPhases([2 1]);
    end

    % Separate by phase (& loop through electrodes for fft)
    for ph = 1:2
        pIdx{ph} = intersect(ERP2.goodTrials,find(round(ERP2.echtPhaseLabel)==ERP2.setPhases(ph)));
    end

    % Get peak and trough ERPs
    terp = squeeze(mean(ERP2.erp(:,:,pIdx{1}),3));
    perp = squeeze(mean(ERP2.erp(:,:,pIdx{2}),3));

    % Get baseline period
    idxBL = find(t < -50 & t > -250);

    % Baseline normalize ERPs
    erp = erp - mean(erp(idxBL,:));
    terp = terp - mean(terp(idxBL,:));
    perp = perp - mean(perp(idxBL,:));

    % Find Latencies
    [rl,ra] = findERPcomp(erp,t);
    [tl,ta] = findERPcomp(terp,t);
    [pl,pa] = findERPcomp(perp,t);

    % Save "P1" latencies
    rlat(k,:) = rl.P50;
    tlat(k,:) = tl.P50;
    plat(k,:) = pl.P50;

    % Grab latency difference
    dtrLat(k,:) = (rlat(k,:) - tlat(k,:));
    dprLat(k,:) = (rlat(k,:) - plat(k,:));

    % Get echt Alphas
    echtR = ERP.ealpha;
    echtT = ERP2.ealpha(:,:,pIdx{1});
    echtP = ERP2.ealpha(:,:,pIdx{2});

    % Get standard hilbert alphas
    sthtR = ERP.halpha;
    sthtT = ERP2.halpha(:,:,pIdx{1});
    sthtP = ERP2.halpha(:,:,pIdx{2});

    % Define phase lag formula
    plag = @(wave,dim) mean(exp(1i * (squeeze(angle(wave(:,1,:))) ...
                                    - squeeze(angle(wave(:,2,:))))),dim);

    % Get Phase lag between electrodes
    echtRlag(k,:) = plag(echtR,2);
    sthtRlag(k,:) = plag(sthtR,2);
    echtTlag(k,:) = plag(echtT,2);
    sthtTlag(k,:) = plag(sthtT,2);
    echtPlag(k,:) = plag(echtP,2);
    sthtPlag(k,:) = plag(sthtP,2);

    % Define phase locking value formula
    plv = @(w1,w2) abs(mean(exp(1i * (w1 - w2))));

    % Get ECHT PLV
    echtplvR(k) = plv(ERP.echtPhase(ERP.goodTrials),deg2rad(ERP.echtPhaseLabel(ERP.goodTrials)));
    echtplvT(k) = plv(ERP2.echtPhase(pIdx{1}),deg2rad(ERP2.echtPhaseLabel(pIdx{1})));
    echtplvP(k) = plv(ERP2.echtPhase(pIdx{2}),deg2rad(ERP2.echtPhaseLabel(pIdx{2})));
    echtplv2(k) = plv(ERP2.echtPhase(ERP2.goodTrials),deg2rad(ERP2.echtPhaseLabel(ERP2.goodTrials)));
    
    % Get Standard hilbert transform PLV
    sthtplvR(k) = plv(ERP.aPhaseOn(ERP.goodTrials,1),deg2rad(ERP.echtPhaseLabel(ERP.goodTrials)));
    sthtplvT(k) = plv(ERP2.aPhaseOn(pIdx{1},1),deg2rad(ERP2.echtPhaseLabel(pIdx{1})));
    sthtplvP(k) = plv(ERP2.aPhaseOn(pIdx{2},1),deg2rad(ERP2.echtPhaseLabel(pIdx{2})));
    sthtplv2(k) = plv(ERP2.aPhaseOn(ERP2.goodTrials,1),deg2rad(ERP2.echtPhaseLabel(ERP2.goodTrials)));

    % Get trough & peak alphas
    eAt = squeeze(mean(ERP2.ealpha(:,:,pIdx{1}),3));
    eAp = squeeze(mean(ERP2.ealpha(:,:,pIdx{2}),3));
    hAt = squeeze(mean(ERP2.halpha(:,:,pIdx{1}),3));
    hAp = squeeze(mean(ERP2.halpha(:,:,pIdx{2}),3));

    % Find time for P50
    for j = 1:2
        if isempty(find(t == tl.P50(j)))
            idxT = [];
            eAphaseT(k,j) = nan;
            hAphaseT(k,j) = nan;
        else
            idxT = find(t == tl.P50(j));
            eAphaseT(k,j) = angle(eAt(idxT));
            hAphaseT(k,j) = angle(hAt(idxT));
        end
        if isempty(find(t == pl.P50(j)))
            idxP = [];
            eAphaseP(k,j) = nan;
            hAphaseP(k,j) = nan;
        else
            idxP = find(t == pl.P50(j));
            eAphaseP(k,j) = angle(eAp(idxP));
            hAphaseP(k,j) = angle(hAp(idxP));
        end
    end

    % Generate figure
    path = './ParticipantEROs';
    tmp = figure;
    str = sprintf('Participant %s',Sub_ID(5:end));
    sgtitle(str)
    subplot(2,1,1)
    plot(t,eAt(:,1),'linewidth',2)
    hold on
    plot([tlat(k,1) tlat(k,1)],[-20 20],...
        '--','color',[0.5 0.5 0.5],'linewidth',2)
    ylim([-15 15])
    yticks([-10 0 10])
    box off 
    px1 = polaraxes('Position',[.7 .7 .2 .2])
    polarplot(px1,[0 wrapTo2Pi(eAphaseT(k,1))],[0 1],...
        '-o','markerfacecolor',co2(1,:),'LineWidth',...
        2,'Color',co2(1,:),'markersize',6);
    rticks(px1,[])
    thetaticks(px1,[])
   
    subplot(2,1,2)
    plot(t,eAp(:,1),'r','linewidth',2)
    hold on
    plot([plat(k,1) plat(k,1)],[-20 20],...
        '--','color',[0.5 0.5 0.5],'linewidth',2)
    ylim([-15 15])
    yticks([-10 0 10])
    box off
    px2 = polaraxes('Position',[.7 .22 .2 .2])
    polarplot(px2,[0 wrapTo2Pi(eAphaseP(k,1))],[0 1],...
        '-o','markerfacecolor',co2(2,:),'LineWidth',...
        2,'Color',co2(2,:),'markersize',6);
    rticks(px2,[])
    thetaticks(px2,[])
    set(gca,'fontsize',14)

    % Save and close
    print(tmp,fullfile(path,sprintf('%s.png',str)),'-dpng')
    % pause(2)
    close all


end

%%



%% Make reliability plots
path = './P1 Analyses';


% Initiate figure
f0 = figure;
scatterDiagHist(rlat(:,1),tlat(:,1),12,co2(1,:),[45 80])
title('FPZ Random-to-Trough Consistency')
xlabel('Random Latency (ms)')
ylabel('Trough Latency (ms)')
set(gca,'fontsize',14)

text(90, 50,sprintf("r = %0.2f",corr(rlat(:,1),tlat(:,1),'rows','complete')),'fontsize',14)

% Fz latency test-retest
f1 = figure;
scatterDiagHist(rlat(:,2),tlat(:,2),12,co2(1,:),[45 80])
title('FZ Random-to-Trough Consistency')
xlabel('Random Latency (ms)')
ylabel('Trough Latency (ms)')
set(gca,'fontsize',14)

text(90, 50,sprintf("r = %0.2f",corr(rlat(:,2),tlat(:,2),'rows','complete')),'fontsize',14)

% Fpz latency peak
f2 = figure;
scatterDiagHist(rlat(:,1),plat(:,1),12,co2(2,:),[45 80])
title('FPZ Random-to-Peak Consistency')
xlabel('Random Latency (ms)')
ylabel('Peak Latency (ms)')
set(gca,'fontsize',14)

text(90, 50,sprintf("r = %0.2f",corr(rlat(:,1),plat(:,1),'rows','complete')),'fontsize',14)

% Fz latency test-retes
f3 = figure;
scatterDiagHist(rlat(:,2),plat(:,2),12,co2(2,:),[45 80])
title('FZ Random-to-Peak Consistency')
xlabel('Random Latency (ms)')
ylabel('Peak Latency (ms)')
set(gca,'fontsize',14)

text(90, 50,sprintf("r = %0.2f",corr(rlat(:,2),plat(:,2),'rows','complete')),'fontsize',14)


print(f0,fullfile(path,'FPZ Trough to Random'),'-dsvg')
print(f0,fullfile(path,'FPZ Trough to Random'),'-dpng')

print(f1,fullfile(path,'FZ Trough to Random'),'-dsvg')
print(f1,fullfile(path,'FZ Trough to Random'),'-dpng')

print(f2,fullfile(path,'FPZ Peak to Random'),'-dsvg')
print(f2,fullfile(path,'FPZ Peak to Random'),'-dpng')

print(f3,fullfile(path,'FZ Peak to Random'),'-dsvg')
print(f3,fullfile(path,'FZ Peak to Random'),'-dpng')


%% Translate to proportion of cycle & plot

% Latency differences as a proportion of alpha cycle
dtrProp = round(dtrLat ./ (1000 ./ iaf'),2);
dprProp = round(dprLat ./ (1000 ./ iaf)',2);

% Vectorize and remove NaN
dtrProp = dtrProp(:,1); % dtrProp = dtrProp(find(~isnan(dtrProp)));
dprProp = dprProp(:,1); % dprProp = dprProp(find(~isnan(dprProp)));

% Generate simulated alpha waves
basis = -0.5:0.001:0.5;
peak = real(cos(2 * pi * basis));
trough = real(cos(2 * pi * basis - pi));

% Trough and Random
f4 = figure;
for j = 1:length(dtrProp)
    [~,idx(j)] = min(abs(basis - dtrProp(j)));
end
plot(basis,trough,'color',co2(1,:),'linewidth',3)
hold on
scatter(basis(idx),trough(idx),120,'markerfacecolor',co2(1,:),'markeredgecolor','none');
plot([0 0],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',2)
ylim([-1.2 1.2])
yticks({})
xlim([-0.5 0.5])
xticks(-0.5:0.1:0.5)
xlabel('Proportions of IAF Cycle')
title('FPZ Latency Difference as Proportion of Alpha Cycle')
set(gca,'fontsize',14)
legend({'','','Target'})
box off
axis square


% Pek and Random
f5 = figure;
for j = 1:length(dprProp)
    [~,idx(j)] = min(abs(basis - dprProp(j)));
end
plot(basis,peak,'color',co2(2,:),'linewidth',3)
hold on
scatter(basis(idx),peak(idx),120,'markerfacecolor',co2(2,:),'markeredgecolor','none');
plot([0 0],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',2)
ylim([-1.2 1.2])
yticks({})
xlim([-0.5 0.5])
xticks(-0.5:0.1:0.5)
xlabel('Proportions of IAF Cycle')
title('FPZ Latency Difference as Proportion of Alpha Cycle')
set(gca,'fontsize',14)
box off
axis square


% Latency differences as a proportion of alpha cycle
dtrProp = round(dtrLat ./ (1000 ./ iaf'),2);
dprProp = round(dprLat ./ (1000 ./ iaf)',2);

% Vectorize and remove NaN
dtrProp = dtrProp(:,2); % dtrProp = dtrProp(find(~isnan(dtrProp)));
dprProp = dprProp(:,2); % dprProp = dprProp(find(~isnan(dprProp)));

% Generate simulated alpha waves
basis = -0.5:0.001:0.5;
peak = real(cos(2 * pi * basis));
trough = real(cos(2 * pi * basis - pi));

% Trough and Random
f6 = figure;
for j = 1:length(dtrProp)
    [~,idx(j)] = min(abs(basis - dtrProp(j)));
end
plot(basis,trough,'color',co2(1,:),'linewidth',3)
hold on
scatter(basis(idx),trough(idx),120,'markerfacecolor',co2(1,:),'markeredgecolor','none');
plot([0 0],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',2)
ylim([-1.2 1.2])
yticks({})
xlim([-0.5 0.5])
xticks(-0.5:0.1:0.5)
xlabel('Proportions of IAF Cycle')
title('FZ Latency Difference as Proportion of Alpha Cycle')
set(gca,'fontsize',14)
legend({'','','Target'})
box off
axis square


% Pek and Random
f7 = figure;
for j = 1:length(dprProp)
    [~,idx(j)] = min(abs(basis - dprProp(j)));
end
plot(basis,peak,'color',co2(2,:),'linewidth',3)
hold on
scatter(basis(idx),peak(idx),120,'markerfacecolor',co2(2,:),'markeredgecolor','none');
plot([0 0],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',2)
ylim([-1.2 1.2])
yticks({})
xlim([-0.5 0.5])
xticks(-0.5:0.1:0.5)
xlabel('Proportions of IAF Cycle')
title('FZ Latency Difference as Proportion of Alpha Cycle')
set(gca,'fontsize',14)
box off
axis square


print(f4,fullfile(path,'FPZ Trough Proportion of Cycle'),'-dsvg')
print(f4,fullfile(path,'FPZ Trough Proportion of Cycle'),'-dpng')

print(f5,fullfile(path,'FPZ Peak Proportion of Cycle'),'-dsvg')
print(f5,fullfile(path,'FPZ Peak Proportion of Cycle'),'-dpng')

print(f6,fullfile(path,'FZ Trough Proportion of Cycle'),'-dsvg')
print(f6,fullfile(path,'FZ Trough Proportion of Cycle'),'-dpng')

print(f7,fullfile(path,'FZ Peak Proportion of Cycle'),'-dsvg')
print(f7,fullfile(path,'FZ Peak Proportion of Cycle'),'-dpng')



%% Just histogram plots

% Setup parameters
edges = [-30:3:30];

% Histogram Trough-Random differences FPZ
f8 = figure;
h0 = histogram(dtrLat(:,1),edges,'FaceColor',co2(1,:),'EdgeColor',co2(1,:));
xlim([-25 25])
ylim([0 6])
yticks([ 1 2 3 4 5 6])
hold off
box off
axis square
title('FPZ Random - Trough Latency')
xlabel('Number of Participants')
xlabel('Latency Difference (ms)')
set(gca,'fontsize',14)

f9 = figure;
h1 = histogram(dtrLat(:,2),edges,'FaceColor',co2(1,:),'EdgeColor',co2(1,:));
xlim([-25 25])
ylim([0 6])
yticks([ 1 2 3 4 5 6])
box off
axis square
title('FZ Random - Trough Latency')
xlabel('Number of Participants')
xlabel('Latency Difference (ms)')
set(gca,'fontsize',14)

% Histogram Peak-Random differences FPZ
f10 = figure;
h3 = histogram(dprLat(:,1),edges,'FaceColor',co2(2,:),'EdgeColor',co2(2,:));
xlim([-25 25])
ylim([0 6])
yticks([ 1 2 3 4 5 6])
box off
axis square
title('FPZ Random - Peak Latency')
xlabel('Number of Participants')
xlabel('Latency Difference (ms)')
set(gca,'fontsize',14)

f11 = figure;
h3 = histogram(dprLat(:,2),edges,'FaceColor',co2(2,:),'EdgeColor',co2(2,:));
xlim([-25 25])
ylim([0 6])
yticks([ 1 2 3 4 5 6])
box off
axis square
title('FZ Random - Peak Latency')
xlabel('Number of Participants')
xlabel('Latency Difference (ms)')
set(gca,'fontsize',14)


print(f8,fullfile(path,'FPZ Random - Trough Hist'),'-dsvg')
print(f8,fullfile(path,'FPZ Random - Trough Hist'),'-dpng')

print(f9,fullfile(path,'FZ Random - Trough Hist'),'-dsvg')
print(f9,fullfile(path,'FZ Random - Trough Hist'),'-dpng')

print(f10,fullfile(path,'FPZ Random - Peak Hist'),'-dsvg')
print(f10,fullfile(path,'FPZ Random - Peak Hist'),'-dpng')

print(f11,fullfile(path,'FZ Random - Peak Hist'),'-dsvg')
print(f11,fullfile(path,'FZ Random - Peak Hist'),'-dpng')


%% Plot Phase lags

% Grab reformated time for fill
time = [t ; flipud(t)];

% ECHT Random Phase lag
f12 = figure;

mu = squeeze(circ_mean(angle(echtRlag)));
se = squeeze(circ_std(angle(echtRlag)))./sqrt(n);
se = [mu + se , fliplr(mu - se)];

plot(t,mu,'color',[0.5 0.5 0.5],'linewidth',2)
hold on
fill(time,se,[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
plot(t,zeros(length(t),1),'--','color',[0.5 0.5 0.5],'linewidth',2)
pbaspect([3 1 1])
xlim([-250 500])
ylim([-0.1 0.1])
xticks([-200 -100 0 100 200 300 400 500])
xlabel('Time (ms)')
ylabel('Phase Lag (radians)')
title('ECHT Random Phase Lag')
set(gca,'fontsize',14)
box off

% Standard Random Phase lag
f13 = figure;

mu = squeeze(circ_mean(angle(sthtRlag)));
se = squeeze(circ_std(angle(sthtRlag)))./sqrt(n);
se = [mu + se , fliplr(mu - se)];

plot(t,mu,'color',[0.5 0.5 0.5],'linewidth',2)
hold on
fill(time,se,[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
plot(t,zeros(length(t),1),'--','color',[0.5 0.5 0.5],'linewidth',2)
pbaspect([3 1 1])
xlim([-250 500])
ylim([-0.1 0.1])
xticks([-200 -100 0 100 200 300 400 500])
xlabel('Time (ms)')
ylabel('Phase Lag (radians)')
title('Standard HT Random Phase Lag')
set(gca,'fontsize',14)
box off

% ECHT Trough Phase lag
f14 = figure;

mu = squeeze(circ_mean(angle(echtTlag)));
se = squeeze(circ_std(angle(echtTlag)))./sqrt(n);
se = [mu + se , fliplr(mu - se)];

plot(t,mu,'color',co2(1,:),'linewidth',2)
hold on
fill(time,se,co2(1,:),'FaceAlpha',0.3,'EdgeColor','none')
plot(t,zeros(length(t),1),'--','color',[0.5 0.5 0.5],'linewidth',2)
pbaspect([3 1 1])
xlim([-250 500])
ylim([-0.1 0.1])
xticks([-200 -100 0 100 200 300 400 500])
xlabel('Time (ms)')
ylabel('Phase Lag (radians)')
title('ECHT Trough Phase Lag')
set(gca,'fontsize',14)
box off


% Standard Trough Phase lag
f15 = figure;

mu = squeeze(circ_mean(angle(sthtTlag)));
se = squeeze(circ_std(angle(sthtTlag)))./sqrt(n);
se = [mu + se , fliplr(mu - se)];

plot(t,mu,'color',co2(1,:),'linewidth',2)
hold on
fill(time,se,co2(1,:),'FaceAlpha',0.3,'EdgeColor','none')
plot(t,zeros(length(t),1),'--','color',[0.5 0.5 0.5],'linewidth',2)
pbaspect([3 1 1])
xlim([-250 500])
ylim([-0.1 0.1])
xticks([-200 -100 0 100 200 300 400 500])
xlabel('Time (ms)')
ylabel('Phase Lag (radians)')
title('Standard HT Trough Phase Lag')
set(gca,'fontsize',14)
box off


% ECHT Peak Phase lag
f16 = figure;

mu = squeeze(circ_mean(angle(echtPlag)));
se = squeeze(circ_std(angle(echtPlag)))./sqrt(n);
se = [mu + se , fliplr(mu - se)];

plot(t,mu,'color',co2(2,:),'linewidth',2)
hold on
fill(time,se,co2(2,:),'FaceAlpha',0.3,'EdgeColor','none')
plot(t,zeros(length(t),1),'--','color',[0.5 0.5 0.5],'linewidth',2)
pbaspect([3 1 1])
xlim([-250 500])
ylim([-0.1 0.1])
xticks([-200 -100 0 100 200 300 400 500])
xlabel('Time (ms)')
ylabel('Phase Lag (radians)')
title('ECHT Peak Phase Lag')
set(gca,'fontsize',14)
box off

% Standard Peak Phase lag
f17 = figure;

mu = squeeze(circ_mean(angle(sthtPlag)));
se = squeeze(circ_std(angle(sthtPlag)))./sqrt(n);
se = [mu + se , fliplr(mu - se)];

plot(t,mu,'color',co2(2,:),'linewidth',2)
hold on
fill(time,se,co2(2,:),'FaceAlpha',0.3,'EdgeColor','none')
plot(t,zeros(length(t),1),'--','color',[0.5 0.5 0.5],'linewidth',2)
pbaspect([3 1 1])
xlim([-250 500])
ylim([-0.1 0.1])
xticks([-200 -100 0 100 200 300 400 500])
xlabel('Time (ms)')
ylabel('Phase Lag (radians)')
title('Standard HT Peak Phase Lag')
set(gca,'fontsize',14)
box off


print(f12,fullfile(path,'ECHT Random Phase Lag'),'-dsvg')
print(f12,fullfile(path,'ECHT Random Phase Lag'),'-dpng')

print(f13,fullfile(path,'Standard HT Random Phase Lag'),'-dsvg')
print(f13,fullfile(path,'Standard HT Random Phase Lag'),'-dpng')

print(f14,fullfile(path,'ECHT Trough Phase Lag'),'-dsvg')
print(f14,fullfile(path,'ECHT Trough Phase Lag'),'-dpng')

print(f15,fullfile(path,'Standard HT Trough Phase Lag'),'-dsvg')
print(f15,fullfile(path,'Standard HT Trough Phase Lag'),'-dpng')

print(f16,fullfile(path,'ECHT Peak Phase Lag'),'-dsvg')
print(f16,fullfile(path,'ECHT Peak Phase Lag'),'-dpng')

print(f17,fullfile(path,'Standard HT Peak Phase Lag'),'-dsvg')
print(f17,fullfile(path,'Standard HT Peak Phase Lag'),'-dpng')


%% PLV reliability

% Initiate figure
f18 = figure;
scatterDiagHist(echtplvR,echtplvT,12,co2(1,:),[0.6 1])
title('ECHT Random-to-Trough PLV Consistency')
xlabel('Random PLV')
ylabel('Trough PLV')
set(gca,'fontsize',14)

[~,pechtRT] = corr(echtplvR',echtplvT','rows','complete')
text(1,0.9,sprintf("r = %0.2f",corr(echtplvR',echtplvT','rows','complete')),'fontsize',14)

% Fz latency test-retest
f19 = figure;
scatterDiagHist(echtplvR,echtplvP,12,co2(2,:),[0.6 1])
title('ECHT Random-to-Peak PLV Consistency')
xlabel('Random PLV')
ylabel('Peak PLV')
set(gca,'fontsize',14)

[~,pechtRP] = corr(echtplvR',echtplvP','rows','complete')
text(1,0.9,sprintf("r = %0.2f",corr(echtplvR',echtplvP','rows','complete')),'fontsize',14)

% Fpz latency peak
f20 = figure;
scatterDiagHist(sthtplvR,sthtplvT,12,co2(1,:),[0.6 1])
title('Standard HT Random-to-Trough Consistency')
xlabel('Random PLV')
ylabel('Trough PLV')
set(gca,'fontsize',14)

[~,psthtRT] = corr(sthtplvR',sthtplvT','rows','complete')
text(1,0.9,sprintf("r = %0.2f",corr(sthtplvR',sthtplvT','rows','complete')),'fontsize',14)

% Fz latency test-retes
f21 = figure;
scatterDiagHist(sthtplvR,sthtplvP,12,co2(2,:),[0.6 1])
title('Standard HT Random-to-Peak Consistency')
xlabel('Random PLV')
ylabel('Peak PLV')
set(gca,'fontsize',14)

[~,psthtRP] = corr(sthtplvR',sthtplvP','rows','complete')
text(1,0.9,sprintf("r = %0.2f",corr(sthtplvR',sthtplvP','rows','complete')),'fontsize',14)


print(f18,fullfile(path,'ECHT Random-Trough PLV Reliability'),'-dsvg')
print(f18,fullfile(path,'ECHT Random-Trough PLV Reliability'),'-dpng')

print(f19,fullfile(path,'ECHT Random-Peak PLV Reliability'),'-dsvg')
print(f19,fullfile(path,'ECHT Random-Peak PLV Reliability'),'-dpng')

print(f20,fullfile(path,'Standard HT Random-Trough PLV Reliability'),'-dsvg')
print(f20,fullfile(path,'Standard HT Random-Trough PLV Reliability'),'-dpng')

print(f21,fullfile(path,'Standard HT Random-Peak PLV Reliability'),'-dsvg')
print(f21,fullfile(path,'Standard HT Random-Peak PLV Reliability'),'-dpng')


%% Random to overall 2Phase


% Initiate figure
f22 = figure;
scatterDiagHist(echtplvR,echtplv2,12,[0.5 0.5 0.5],[0.6 1])
title('ECHT Random-to-2Phase PLV Consistency')
xlabel('Random PLV')
ylabel('Trough PLV')
set(gca,'fontsize',14)

[~,pechtR2] = corr(echtplvR',sthtplv2','rows','complete')
text(1,0.9,sprintf("r = %0.2f",corr(echtplvR',echtplv2','rows','complete')),'fontsize',14)

% Initiate figure
f23 = figure;
scatterDiagHist(sthtplvR,sthtplv2,12,[0.5 0.5 0.5],[0.6 1])
title('Standard HT Random-to-2Phase PLV Consistency')
xlabel('Random PLV')
ylabel('Trough PLV')
set(gca,'fontsize',14)

[~,psthtR2] = corr(sthtplvR',sthtplv2','rows','complete')
text(1,0.9,sprintf("r = %0.2f",corr(sthtplvR',sthtplv2','rows','complete')),'fontsize',14)

print(f22,fullfile(path,'ECHT Random-2Phase PLV Reliability'),'-dsvg')
print(f22,fullfile(path,'ECHT Random-2Phase PLV Reliability'),'-dpng')

print(f23,fullfile(path,'Standard HT Random-2Phase PLV Reliability'),'-dsvg')
print(f23,fullfile(path,'Standard HT Random-2Phase PLV Reliability'),'-dpng')

%% P1 Latency over alpha Phase

[ptrough,utrough] = watsons_U2_approx_p(eAphaseT(~isnan(eAphaseT(:,1)),1),eAphaseT(~isnan(eAphaseT(:,1)),2));
[ppeak,upeak] = watsons_U2_approx_p(eAphaseT(~isnan(eAphaseP(:,1)),1),eAphaseT(~isnan(eAphaseP(:,1)),2));
p = [ptrough ppeak];
h = benhoch(p);

% ECHT
% Trough  FPZ
f24 = figure;
ax1 = polaraxes;
hold on
for i = 1:n
polarplot(ax1,[0 wrapTo2Pi(eAphaseT(i,1))],[0 echtplvT(i)],...
    'o-','LineWidth',2,'Color',[co2(1,:) 0.5],'MarkerFaceColor','w');
end

% get non nan values
idx = eAphaseT(~isnan(eAphaseT(:,1)),1);

% Mean angle and PLV
polarplot(ax1,[0 circ_mean(idx)],[0 1],...
    'LineWidth',4,'Color',co2(1,:));

% SD angle arc
polarplot(ax1,circ_mean(idx) + linspace(-circ_std(idx),circ_std(idx),100),...
    ones(size(linspace(-circ_std(idx),circ_std(idx),100))),...
    'LineWidth',2,'Color',co2(1,:));

ax1.Title.String = 'ECHT Trough Phase of P1 FPZ';
ax1.RAxisLocation = 70;
ax1.FontSize = 14;
rlim(ax1, [0 1]);
rticks(ax1,[0 0.2 0.4 0.6 0.8 1])

% Peak FPZ
f25 = figure;
ax2 = polaraxes;
hold on
for i = 1:n
polarplot(ax2,[0 wrapTo2Pi(eAphaseP(i,1))],[0 echtplvP(i)],...
    'o-','LineWidth',2,'Color',[co2(2,:) 0.5],'MarkerFaceColor','w');
end

% get non nan values
idx = eAphaseP(~isnan(eAphaseP(:,1)),1);

% Mean angle and PLV
polarplot(ax2,[0 circ_mean(idx)],[0 1],...
    'LineWidth',4,'Color',co2(2,:));

% SD angle arc
polarplot(ax2,circ_mean(idx) + linspace(-circ_std(idx),circ_std(idx),100),...
    ones(size(linspace(-circ_std(idx),circ_std(idx),100))),...
    'LineWidth',2,'Color',co2(2,:));

ax2.Title.String = 'ECHT Peak Phase of P1 FPZ';
ax2.RAxisLocation = 70;
ax2.FontSize = 14;
rlim(ax2, [0 1]);
rticks(ax2,[0 0.2 0.4 0.6 0.8 1])

% Trough  FZ
f26 = figure;
ax3 = polaraxes;
hold on
for i = 1:n
polarplot(ax3,[0 wrapTo2Pi(eAphaseT(i,2))],[0 echtplvT(i)],...
    'o-','LineWidth',2,'Color',[co2(1,:) 0.5],'MarkerFaceColor','w');
end

% get non nan values
idx = eAphaseT(~isnan(eAphaseT(:,2)),2);

% Mean angle and PLV
polarplot(ax3,[0 circ_mean(idx)],[0 1],...
    'LineWidth',4,'Color',co2(1,:));

% SD angle arc
polarplot(ax3,circ_mean(idx) + linspace(-circ_std(idx),circ_std(idx),100),...
    ones(size(linspace(-circ_std(idx),circ_std(idx),100))),...
    'LineWidth',2,'Color',co2(1,:));

ax3.Title.String = 'ECHT Trough Phase of P1 FZ';
ax3.RAxisLocation = 70;
ax3.FontSize = 14;
rlim(ax3, [0 1]);
rticks(ax3,[0 0.2 0.4 0.6 0.8 1])


% Peak FZ
f27 = figure;
ax4 = polaraxes;
hold on
for i = 1:n
polarplot(ax4,[0 wrapTo2Pi(eAphaseP(i,2))],[0 echtplvP(i)],...
    'o-','LineWidth',2,'Color',[co2(2,:) 0.5],'MarkerFaceColor','w');
end

% get non nan values
idx = eAphaseP(~isnan(eAphaseP(:,2)),2);

% Mean angle and PLV
polarplot(ax4,[0 circ_mean(idx)],[0 1],...
    'LineWidth',4,'Color',co2(2,:));

% SD angle arc
polarplot(ax4,circ_mean(idx) + linspace(-circ_std(idx),circ_std(idx),100),...
    ones(size(linspace(-circ_std(idx),circ_std(idx),100))),...
    'LineWidth',2,'Color',co2(2,:));

ax4.Title.String = 'ECHT Peak Phase of P1 FZ';
ax4.RAxisLocation = 70;
ax4.FontSize = 14;
rlim(ax4, [0 1]);
rticks(ax4,[0 0.2 0.4 0.6 0.8 1])

% Standard HT
% Trough  FPZ
f28 = figure;
ax5 = polaraxes;
hold on
for i = 1:n
polarplot(ax5,[0 wrapTo2Pi(hAphaseT(i,1))],[0 sthtplvT(i)],...
    'o-','LineWidth',2,'Color',[co2(1,:) 0.5],'MarkerFaceColor','w');
end

% get non nan values
idx = hAphaseT(~isnan(hAphaseT(:,1)),1);

% Mean angle and PLV
polarplot(ax5,[0 circ_mean(idx)],[0 1],...
    'LineWidth',4,'Color',co2(1,:));

% SD angle arc
polarplot(ax5,circ_mean(idx) + linspace(-circ_std(idx),circ_std(idx),100),...
    ones(size(linspace(-circ_std(idx),circ_std(idx),100))),...
    'LineWidth',2,'Color',co2(1,:));

ax5.Title.String = 'Standard HT Trough Phase of P1 FPZ';
ax5.RAxisLocation = 70;
ax5.FontSize = 14;
rlim(ax5, [0 1]);
rticks(ax5,[0 0.2 0.4 0.6 0.8 1])

% Peak FPZ
f29 = figure;
ax6 = polaraxes;
hold on
for i = 1:n
polarplot(ax6,[0 wrapTo2Pi(hAphaseP(i,1))],[0 sthtplvP(i)],...
    'o-','LineWidth',2,'Color',[co2(2,:) 0.5],'MarkerFaceColor','w');
end

% get non nan values
idx = hAphaseP(~isnan(hAphaseP(:,1)),1);

% Mean angle and PLV
polarplot(ax6,[0 circ_mean(idx)],[0 1],...
    'LineWidth',4,'Color',co2(2,:));

% SD angle arc
polarplot(ax6,circ_mean(idx) + linspace(-circ_std(idx),circ_std(idx),100),...
    ones(size(linspace(-circ_std(idx),circ_std(idx),100))),...
    'LineWidth',2,'Color',co2(2,:));

ax6.Title.String = 'Standard HT Peak Phase of P1 FPZ';
ax6.RAxisLocation = 70;
ax6.FontSize = 14;
rlim(ax6, [0 1]);
rticks(ax6,[0 0.2 0.4 0.6 0.8 1])

% Trough  FZ
f30 = figure;
ax7 = polaraxes;
hold on
for i = 1:n
polarplot(ax7,[0 wrapTo2Pi(hAphaseT(i,2))],[0 sthtplvT(i)],...
    'o-','LineWidth',2,'Color',[co2(1,:) 0.5],'MarkerFaceColor','w');
end

% get non nan values
idx = hAphaseT(~isnan(hAphaseT(:,2)),2);

% Mean angle and PLV
polarplot(ax7,[0 circ_mean(idx)],[0 1],...
    'LineWidth',4,'Color',co2(1,:));


% SD angle arc
polarplot(ax7,circ_mean(idx) + linspace(-circ_std(idx),circ_std(idx),100),...
    ones(size(linspace(-circ_std(idx),circ_std(idx),100))),...
    'LineWidth',2,'Color',co2(1,:));

ax7.Title.String = 'Standard HT Trough Phase of P1 FZ';
ax7.RAxisLocation = 70;
ax7.FontSize = 14;
rlim(ax7, [0 1]);
rticks(ax7,[0 0.2 0.4 0.6 0.8 1])


% Peak FZ
f31 = figure;
ax8 = polaraxes;
hold on
for i = 1:n
polarplot(ax8,[0 wrapTo2Pi(hAphaseP(i,2))],[0 sthtplvP(i)],...
    'o-','LineWidth',2,'Color',[co2(2,:) 0.5],'MarkerFaceColor','w');
end

% get non nan values
idx = hAphaseP(~isnan(hAphaseP(:,2)),2);

% Mean angle and PLV
polarplot(ax8,[0 circ_mean(idx)],[0 1],...
    'LineWidth',4,'Color',co2(2,:));

% SD angle arc
polarplot(ax8,circ_mean(idx) + linspace(-circ_std(idx),circ_std(idx),100),...
    ones(size(linspace(-circ_std(idx),circ_std(idx),100))),...
    'LineWidth',2,'Color',co2(2,:));

ax8.Title.String = 'Standard HT Peak Phase of P1 FZ';
ax8.RAxisLocation = 70;
ax8.FontSize = 14;
rlim(ax8, [0 1]);
rticks(ax8,[0 0.2 0.4 0.6 0.8 1])


print(f24,fullfile(path,'ECHT Alpha Phase at P1 Trough FPZ'),'-dsvg')
print(f24,fullfile(path,'ECHT Alpha Phase at P1 Trough FPZ'),'-dpng')

print(f25,fullfile(path,'ECHT Alpha Phase at P1 Peak FPZ'),'-dsvg')
print(f25,fullfile(path,'ECHT Alpha Phase at P1 Peak FPZ'),'-dpng')

print(f26,fullfile(path,'ECHT Alpha Phase at P1 Trough FZ'),'-dsvg')
print(f26,fullfile(path,'ECHT Alpha Phase at P1 Trough FZ'),'-dpng')

print(f27,fullfile(path,'ECHT Alpha Phase at P1 Peak FZ'),'-dsvg')
print(f27,fullfile(path,'ECHT Alpha Phase at P1 Peak FZ'),'-dpng')

print(f28,fullfile(path,'Standard HT Alpha Phase at P1 Trough FPZ'),'-dsvg')
print(f28,fullfile(path,'Standard HT Alpha Phase at P1 Trough FPZ'),'-dpng')

print(f29,fullfile(path,'Standard HT Alpha Phase at P1 Peak FPZ'),'-dsvg')
print(f29,fullfile(path,'Standard HT Alpha Phase at P1 Peak FPZ'),'-dpng')

print(f30,fullfile(path,'Standard HT Alpha Phase at P1 Trough FZ'),'-dsvg')
print(f30,fullfile(path,'Standard HT Alpha Phase at P1 Trough FZ'),'-dpng')

print(f31,fullfile(path,'Standard HT Alpha Phase at P1 Peak FZ'),'-dsvg')
print(f31,fullfile(path,'Standard HT Alpha Phase at P1 Peak FZ'),'-dpng')

%% P1 latency relative to target

% Generate simulated alpha waves
basis = -0.5:0.001:0.5;
trough = real(cos(2 * pi * basis));
peak = real(cos(2 * pi * basis - pi));

% New basis
basis = linspace(-pi,pi,length(basis));

% Trough and Random
f32 = figure;
for j = 1:length(dtrProp)
    [~,idx(j)] = min(abs(basis - wrapToPi(eAphaseT(j,1) - 0)));
end
plot(basis,trough,'color',co2(1,:),'linewidth',3)
hold on
scatter(basis(idx),trough(idx),120,'markerfacecolor',co2(1,:),'markeredgecolor','none');
plot([0 0],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',2)
ylim([-1.2 1.2])
yticks({})
xlim([-pi pi])
xticks(linspace(-pi,pi,11))
xticklabels({'-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5'})
xlabel('Distance from Target (proportion of cycle)')
title('FPZ Trough P1 Latency Relative to Target Phase')
set(gca,'fontsize',14)
legend({'','','Target'})
box off
axis square


% Peak and Random
f33 = figure;
for j = 1:length(dprProp)
    [~,idx(j)] = min(abs(basis - wrapToPi(eAphaseP(j,1) - pi)));
end
plot(basis,peak,'color',co2(2,:),'linewidth',3)
hold on
scatter(basis(idx),peak(idx),120,'markerfacecolor',co2(2,:),'markeredgecolor','none');
plot([0 0],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',2)
ylim([-1.2 1.2])
yticks({})
xlim([-pi pi])
xticks(linspace(-pi,pi,11))
xticklabels({'-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5'})
xlabel('Distance from Target (proportion of cycle)')
title('FPZ Peak P1 Latency Relative to Target Phase')
set(gca,'fontsize',14)
box off
axis square

% Trough and Random
f34 = figure;
for j = 1:length(dtrProp)
    [~,idx(j)] = min(abs(basis - wrapToPi(eAphaseT(j,2) - 0)));
end
plot(basis,trough,'color',co2(1,:),'linewidth',3)
hold on
scatter(basis(idx),trough(idx),120,'markerfacecolor',co2(1,:),'markeredgecolor','none');
plot([0 0],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',2)
ylim([-1.2 1.2])
yticks({})
xlim([-pi pi])
xticks(linspace(-pi,pi,11))
xticklabels({'-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5'})
xlabel('Distance from Target (proportion of cycle)')
title('FZ Trough P1 Latency Relative to Target Phase')
set(gca,'fontsize',14)
legend({'','','Target'})
box off
axis square


% Peak and Random
f35 = figure;
for j = 1:length(dprProp)
    [~,idx(j)] = min(abs(basis - wrapToPi(eAphaseP(j,2) - pi)));
end
plot(basis,peak,'color',co2(2,:),'linewidth',3)
hold on
scatter(basis(idx),peak(idx),120,'markerfacecolor',co2(2,:),'markeredgecolor','none');
plot([0 0],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',2)
ylim([-1.2 1.2])
yticks({})
xlim([-pi pi])
xticks(linspace(-pi,pi,11))
xticklabels({'-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5'})
xlabel('Distance from Target (proportion of cycle)')
title('FZ Peak P1 Latency Relative to Target Phase')
set(gca,'fontsize',14)
box off
axis square


print(f32,fullfile(path,'FPZ Trough P1 Latency Relative to Target Phase'),'-dsvg')
print(f32,fullfile(path,'FPZ Trough P1 Latency Relative to Target Phase'),'-dpng')

print(f33,fullfile(path,'FPZ Peak P1 Latency Relative to Target Phase'),'-dsvg')
print(f33,fullfile(path,'FPZ Peak P1 Latency Relative to Target Phase'),'-dpng')

print(f34,fullfile(path,'FZ Trough P1 Latency Relative to Target Phase'),'-dsvg')
print(f34,fullfile(path,'FZ Trough P1 Latency Relative to Target Phase'),'-dpng')

print(f35,fullfile(path,'FZ Peak P1 Latency Relative to Target Phase'),'-dsvg')
print(f35,fullfile(path,'FZ Peak P1 Latency Relative to Target Phase'),'-dpng')

