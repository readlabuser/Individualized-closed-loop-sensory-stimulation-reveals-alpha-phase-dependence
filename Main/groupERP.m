% groupERP Analysis
%
% Group level analysis of Random & Two-Phase ERP data
%
% 2021-04-08 scott.bressler@elemindtech.com
% Updated 9-15-2023 tylor.harlow@uconn.edu

%% Housekeeping
clear;
close all;
clc;

%% Load Data

% Select subjects to exclude from group analysis
exSub = {'Sub_151','Sub_152','126'}; %,'Sub_102','Sub_104','Sub_126','Sub_138'};

%% Analysis parameters

% Define pre-processing for raw EEG data
params = echtparams('filterImplementation',1,'fslide',0,'tWin',[-0.250 0.500]);

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

% Get Data Paths
pname = genpath('/Applications/Toolbox/SubjectData'); 
pname = strsplit(pname,':')';

% Separate to Random & TwoPhase
pname1 = pname(contains(pname,'ISI750ms85dB') & ~contains(pname,'(Low'));
pname2 = pname(contains(pname,'TwoPhaseERP85dB') & ~contains(pname,'(Low'));

% Exclude partiicpants with incomplete data
if(~isempty(exSub))  
    pname1 = pname1(~contains(pname1,exSub)); % exclude subject(s) from analysis
    pname2 = pname2(~contains(pname2,exSub)); % exclude subject(s) from analysis
end

% Whats my sample?
n = length(pname1);

% Loop through directories
for k = 1:n

    % Get subject info
    tmp = strsplit(pname1{k},filesep);
    Sub_ID = tmp{contains(tmp,'Sub_')};
    subIdx = contains(dataTable.SubID,Sub_ID([5:end]));

    % Setup alpha center frequency
    params.alphaCF = dataTable(subIdx,:).IAF(1);

    % Process Random
    [ERP1(k),params1{k}] = batchERP(pname1{k},params);

    % Process TwoPhase
    [ERP2(k),params2{k}] = batchERP(pname2{k},params);

    % Store time
    t = ERP1(1).t0;

    % Check order of optimal/pessimal phase list in ERP structure against the
    % calculated phases in 'dataTable.' Fix wrap-around error if necessary
    if(dataTable(subIdx,:).OptPhase > dataTable(subIdx,:).PesPhase)
        ERP2(k).setPhases = ERP2(k).setPhases([2 1]);
    end

    % Sort Trial Type
    for ph = 1:2
        pIdx{ph} = intersect(ERP2(k).goodTrials,find(round(ERP2(k).echtPhaseLabel)==ERP2(k).setPhases(ph)));
    end

    % Grab ERPs
    ERP = squeeze(mean(permute(ERP1(k).erp(:,:,ERP1(k).goodTrials),[3,1,2])));  
    ERPtrough = squeeze(mean(permute(ERP2(k).erp(:,:,pIdx{1}),[3,1,2])));
    ERPpeak = squeeze(mean(permute(ERP2(k).erp(:,:,pIdx{2}),[3,1,2])));
    
    % Get Baseline period
    idxBL = find(t< -50 & t>-250);

    % Baseline Noramalize
    erp(k,:,:) = ERP - mean(ERP(idxBL,:)); 
    trough(k,:,:) = ERPtrough - mean(ERPtrough(idxBL,:)); 
    peak(k,:,:) = ERPpeak - mean(ERPpeak(idxBL,:)); 
    
    % Extract ERP component latency & amplitude:
    % Random
    [lat,amp] = findERPcomp(ERP,t);

    % Latency            % Amplitudes
    P0r(k,:) = lat.P30;  P0ra(k,:) = amp.P30;   
    P1r(k,:) = lat.P50;  P1ra(k,:) = amp.P50;
    N1r(k,:) = lat.N100; N1ra(k,:) = amp.N100;
    P2r(k,:) = lat.P2;   P2ra(k,:) = amp.P2;
    N2r(k,:) = lat.N2;   N2ra(k,:) = amp.N2;
    clearvars lat amp

    % Trough
    [lat,amp] = findERPcomp(ERPtrough,t);

    % Latency            % Amplitudes
    P0t(k,:) = lat.P30;  P0ta(k,:) = amp.P30;   
    P1t(k,:) = lat.P50;  P1ta(k,:) = amp.P50;
    N1t(k,:) = lat.N100; N1ta(k,:) = amp.N100;
    P2t(k,:) = lat.P2;   P2ta(k,:) = amp.P2;
    N2t(k,:) = lat.N2;   N2ta(k,:) = amp.N2;
    clearvars lat amp

    % Peak
    [lat,amp] = findERPcomp(ERPpeak,t); 

    % Latency            % Amplitudes
    P0p(k,:) = lat.P30;  P0pa(k,:) = amp.P30;   
    P1p(k,:) = lat.P50;  P1pa(k,:) = amp.P50;
    N1p(k,:) = lat.N100; N1pa(k,:) = amp.N100;
    P2p(k,:) = lat.P2;   P2pa(k,:) = amp.P2;
    N2p(k,:) = lat.N2;   N2pa(k,:) = amp.N2;
    clearvars lat amp

end


%% Post-process

% Means
troughmu = squeeze(mean(trough));
peakmu = squeeze(mean(peak));
erpmu = squeeze(mean(erp));

% Permutation and SEM
[ptvp,troughsem,peaksem] = npperm(trough,peak,1000,"cluster", 0, 0);
[~,~,erpsem] = npperm(trough,erp,1000,"cluster", 0, 0);

%% Figures

% Define Colors
gray = [0.5 0.5 0.5];
co2 = [0.0000 0.4470 0.7410
       0.8500 0.10 0.0980];

% ERP Component Anlaysis windows
p0win = t(t >= 0 & t <= 38);
p1win = t(t >= 39 & t <= 80);
n1win = t(t >= 80 & t <= 120);
p2win = t(t >= 120 & t <= 220);
n2win = t(t >= 220 & t <= 400);

% Wrap-around time
time = [t',fliplr(t')];

% Peak vs Trough
f1 = figure;

% Fpz
subplot(6,2,[1 3 5 7 9])
p1 = plot(t,troughmu(:,1),'color',co2(1,:),'linewidth',2);
hold on
p2 = plot(t,peakmu(:,1),'color',co2(2,:),'linewidth',2);
fill(time,peaksem(:,1),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
fill(time,troughsem(:,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([min(t) max(t)])
ylim([-6 9])
ylabel('Voltage (uV)')
box off
axis square
ax.FontSize = 9; 
ax = gca;
title('Fpz Trough v. Peak ERP')
legend('Trough','Peak')

% Fpz signifigance
subplot(6,2,11)
plot(t,ptvp(:,1),'k-')
hold on
fill([t, fliplr(t)],[ptvp(:,1), ones(length(ptvp(:,1)),1)],'k')
fill([p0win; flipud(p0win)], [zeros(1,length(p0win)), ones(1,length(p0win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([p1win; flipud(p1win)], [zeros(1,length(p1win)), ones(1,length(p1win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([n1win; flipud(n1win)], [zeros(1,length(n1win)), ones(1,length(n1win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([p2win; flipud(p2win)], [zeros(1,length(p2win)), ones(1,length(p2win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([n2win; flipud(n2win)], [zeros(1,length(n2win)), ones(1,length(n2win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
set(gca,'Ydir','reverse')
ylabel('p-value')
ylim([0 0.05])
yticks([0 0.02 0.04])
box off
xlim([min(t) max(t)])
xlabel('Time (ms)')
ax = gca;
ax.FontSize = 9; 

% Fz
subplot(6,2,[2 4 6 8 10])
p3 = plot(t,troughmu(:,2),'color',co2(1,:),'linewidth',2);
hold on
p4 = plot(t,peakmu(:,2),'color',co2(2,:),'linewidth',2);
fill(time,peaksem(:,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
fill(time,troughsem(:,2),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([min(t) max(t)])
ylim([-6 9])
ax = gca;
ax.FontSize = 9; 
ylabel('Voltage (uV)')
title('Fz Trough v. Peak ERP')
box off
axis square
legend('Trough','Peak')

% Fz signifigance
subplot(6,2,12)
plot(t,ptvp(:,2),'k-')
hold on
fill([p0win; flipud(p0win)], [zeros(1,length(p0win)), ones(1,length(p0win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([p1win; flipud(p1win)], [zeros(1,length(p1win)), ones(1,length(p1win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([n1win; flipud(n1win)], [zeros(1,length(n1win)), ones(1,length(n1win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([p2win; flipud(p2win)], [zeros(1,length(p2win)), ones(1,length(p2win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([n2win; flipud(n2win)], [zeros(1,length(n2win)), ones(1,length(n2win))],...
    [0.5 0.5 0.5], 'edgecolor','k','facealpha', 0.1)
fill([t, fliplr(t)],[ptvp(:,2), ones(length(ptvp(:,2)),1)],'k')
set(gca,'Ydir','reverse')
ylabel('p-value')
ylim([0 0.05])
yticks([0 0.02 0.04])
ax = gca;
ax.FontSize = 9; 
box off
xlim([min(t) max(t)])
xlabel('Time (ms)')


% Random ERP Figures
f2 = figure;

% Fpz
subplot(6,2,[1 3 5 7 9 11])
plot(t,erpmu(:,1),'k-','linewidth',2)
hold on
fill(time,erpsem(:,1),gray, ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([min(t) max(t)])
ylim([-6 8])
ax = gca;
ax.FontSize = 9; 
ylabel('Voltage (uV)')
title('Fpz ERP - baseline corrected t<0 all trials')
box off
axis square

% Fz
subplot(6,2,[2 4 6 8 10 12])
plot(t,erpmu(:,2),'k-','linewidth',2)
hold on
fill(time,erpsem(:,2),gray, ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlim([min(t) max(t)])
ylim([-6 8])
ax = gca;
ax.FontSize = 9; 
ylabel('Voltage (uV)')
title('Fz ERP - baseline corrected t<0 all trials')
box off
axis square


%% Beewarm/Violin Plots

% Grab nice colors
co2 = [0.0000 0.4470 0.7410
       0.8500 0.10 0.0980];

% Definite Color & Postitions
posl =    [.7  1.8  2.9   4   5.1]*4;
posr =    [1   2.1  3.2  4.3  5.4]*4;
posrand = [1.3 2.4  3.5  4.6  5.7]*4;
colr = [co2(2,:); co2(2,:); co2(2,:); co2(2,:); co2(2,:)];
coll = [co2(1,:); co2(1,:); co2(1,:); co2(1,:); co2(1,:)];
colrand = [0.5 0.5 0.5];

% Structure arrays for figures
fpzlatenciesl = [P0t(:,1) , P1t(:,1) , N1t(:,1) , P2t(:,1) , N2t(:,1)];
fzlatenciesl =  [P0t(:,2) , P1t(:,2) , N1t(:,2) , P2t(:,2) , N2t(:,2)];
fpzampl =   abs([P0ta(:,1), P1ta(:,1), N1ta(:,1), P2ta(:,1), N2ta(:,1)]);
fzampl =    abs([P0ta(:,2), P1ta(:,2), N1ta(:,2), P2ta(:,2), N2ta(:,2)]);

fpzlatenciesr = [P0p(:,1) , P1p(:,1) , N1p(:,1) , P2p(:,1) , N2p(:,1)];
fzlatenciesr =  [P0p(:,2) , P1p(:,2) , N1p(:,2) , P2p(:,2) , N2p(:,2)];
fpzampr =   abs([P0pa(:,1), P1pa(:,1), N1pa(:,1), P2pa(:,1), N2pa(:,1)]);
fzampr =    abs([P0pa(:,2), P1pa(:,2), N1pa(:,2), P2pa(:,2), N2pa(:,2)]);

fpzlatenciesrand = [P0r(:,1) , P1r(:,1) , N1r(:,1) , P2r(:,1) , N2r(:,1)];
fzlatenciesrand =  [P0r(:,2) , P1r(:,2) , N1r(:,2) , P2r(:,2) , N2r(:,2)];
fpzamprand =   abs([P0ra(:,1), P1ra(:,1), N1ra(:,1), P2ra(:,1), N2ra(:,1)]);
fzamprand =    abs([P0ra(:,2), P1ra(:,2), N1ra(:,2), P2ra(:,2), N2ra(:,2)]);


% Latency Beeswarm plot
f3 = figure;

% FPZ
subplot(1,2,1)
[~, ~, ~, ~, ~, ~, xl, ~] = al_goodplot(fpzlatenciesl,posl,0.25,coll,'left');
hold on
[~, ~, ~, ~, ~,~, xr, ~] = al_goodplot(fpzlatenciesr,posr,0.25,colr,'left');
[~, ~, ~, ~, ~, ~,xrand, ~] = al_goodplot(fpzlatenciesrand,posrand,0.25,colrand,'left');
for i = 1:5
    plot([xl(:,i), xr(:,i)]', [fpzlatenciesl(:,i), fpzlatenciesr(:,i)]','color',[0.5 0.5 0.5 0.1],'linewidth',2)
    plot([xrand(:,i), xr(:,i)]', [fpzlatenciesrand(:,i), fpzlatenciesr(:,i)]','color',[0.5 0.5 0.5 0.1],'linewidth',2)
end
box off
xticks([posr])
xlim([1 24])
xticklabels({'P0','P1','N1','P2','N2'})
title('FPZ')
ylabel('ERP Latency (ms)')
pbaspect([1.68 1 1])
ylim([-10 400])
grid off

% FZ
subplot(1,2,2)
[~, ~, ~, ~, ~, ~, xl, ~] = al_goodplot(fzlatenciesl,posl,0.25,coll,'left');
hold on
[~, ~, ~, ~, ~,~, xr, ~] = al_goodplot(fzlatenciesr,posr,0.25,colr,'left');
[~, ~, ~, ~, ~, ~,xrand, ~] = al_goodplot(fzlatenciesrand,posrand,0.25,colrand,'left');
for i = 1:5
    plot([xl(:,i), xr(:,i)]', [fzlatenciesl(:,i), fzlatenciesr(:,i)]','color',[0.5 0.5 0.5 0.1],'linewidth',2)
    plot([xrand(:,i), xr(:,i)]', [fzlatenciesrand(:,i), fzlatenciesr(:,i)]','color',[0.5 0.5 0.5 0.1],'linewidth',2)
end
box off
xticks([posr])
xlim([1 24])
xticklabels({'P0','P1','N1','P2','N2'})
title('FZ')
pbaspect([1.68 1 1])
ylim([-10 400])
grid off

% Amplitude swarmolin plot
f4 = figure;

% FPZ
subplot(1,2,1)
[~, ~, ~, ~, ~, ~, xl, ~] = al_goodplot(fpzampl,posl,0.25,coll,'left');
hold on
[~, ~, ~, ~, ~,~, xr, ~] = al_goodplot(fpzampr,posr,0.25,colr,'left');
[~, ~, ~, ~, ~, ~,xrand, ~] = al_goodplot(fpzamprand,posrand,0.25,colrand,'left');
for i = 1:5
    plot([xl(:,i), xr(:,i)]', [fpzampl(:,i), fpzampr(:,i)]','color',[0.5 0.5 0.5 0.1],'linewidth',2)
    plot([xrand(:,i), xr(:,i)]', [fpzamprand(:,i), fpzampr(:,i)]','color',[0.5 0.5 0.5 0.1],'linewidth',2)
end
box off
xticks(([1 2.1 3.2 4.3 5.4] - 0.15)*4)
xticklabels({'P0','P1','N1','P2','N2'})
title('FPZ')
ylabel('ERP Amplitude (uV)')
ylim([-5 15])
xlim([1 24])
grid off
pbaspect([1.68 1 1])

% FZ
subplot(1,2,2)
[~, ~, ~, ~, ~, ~, xl, ~] = al_goodplot(fzampl,posl,0.25,coll,'left');
hold on
[~, ~, ~, ~, ~,~, xr, ~] = al_goodplot(fzampr,posr,0.25,colr,'left');
[~, ~, ~, ~, ~, ~,xrand, ~] = al_goodplot(fzamprand,posrand,0.25,colrand,'left');
for i = 1:5
    plot([xl(:,i), xr(:,i)]', [fzampl(:,i), fzampr(:,i)]','color',[0.5 0.5 0.5 0.1],'linewidth',2)
    plot([xrand(:,i), xr(:,i)]', [fpzamprand(:,i), fzampr(:,i)]','color',[0.5 0.5 0.5 0.1],'linewidth',2)
end
box off
xticks(([1 2.1 3.2 4.3 5.4] - 0.15)*4)
xticklabels({'P0','P1','N1','P2','N2'})
ylim([-5 15])
xlim([1 24])
title('FZ')
pbaspect([1.68 1 1])
grid off


%% ANOVA for anlyses

% Aggregate
fpzamp = abs([[P0ta(:,1); P0pa(:,1); P0ra(:,1)] ; [P1ta(:,1); P1pa(:,1); P1ra(:,1)] ; [N1ta(:,1); N1pa(:,1); N1ra(:,1)] ; [P2ta(:,1); P2pa(:,1); P2ra(:,1)] ; [N2ta(:,1); N2pa(:,1); N2ra(:,1)]]);
fzamp =  abs([[P0ta(:,2); P0pa(:,2); P0ra(:,2)] ; [P1ta(:,2); P1pa(:,2); P1ra(:,2)] ; [N1ta(:,2); N1pa(:,2); N1ra(:,2)] ; [P2ta(:,2); P2pa(:,2); P2ra(:,2)] ; [N2ta(:,2); N2pa(:,2); N2ra(:,2)]]);
amplitude = [fpzamp ; fzamp];

% Setup Amplitude ANOVA
electrode = repelem(["fpz", "fz"], [285 285])';
condition = repmat(repelem(["Trough", "Peak", "Random"], [19 19 19]), [1 5])';
condition = [condition ; condition];
component = repelem(["P0", "P1", "N1", "P2", "N2"], [57 57 57 57 57])';
component = [component ; component];

% Amplitude ANOVA
[~,amptbl] = anovan(amplitude, {condition, electrode, component},...
    'varnames',{'Condition', 'Electrode','Component'},'display','off');

% Aggregate
fpzlatencies = [[P0t(:,1); P0p(:,1); P0r(:,1)]; [P1t(:,1);P1p(:,1); P1r(:,1)];[N1t(:,1);N1p(:,1); N1r(:,1)];[P2t(:,1); P2p(:,1); P2r(:,1)];[ N2t(:,1); N2p(:,1); N2r(:,1)]];
fzlatencies =  [[P0t(:,2); P0p(:,2); P0r(:,2)]; [P1t(:,2);P1p(:,2); P1r(:,2)];[N1t(:,2);N1p(:,2); N1r(:,2)];[P2t(:,2); P2p(:,2); P2r(:,2)];[N2t(:,2); N2p(:,2); N2r(:,2)]];
latency =      [ fpzlatencies ; fzlatencies ];

% Multiway anova
[~,lattbl] = anovan(latency, {condition, electrode, component},...
    'varnames',{'Condition', 'Electrode','Component'},'display','off');

% Ttest amplitudes
[~,p0] = ttest(P0ta - P0pa);
[~,p1] = ttest(P1ta - P1pa);
[~,n1] = ttest(N1ta - N1pa);
[~,p2] = ttest(P2ta - P2pa);
[~,n2] = ttest(N2ta - N2pa);
varnames = {'Component', 'fpz t-test p-value','fz t-test p-value'};
comp = ["P0", "P1", "N1", "P2", "N2"]';
ps = [p0; p1; n1; p2; n2];
[~,~,~,p]=fdr_bh(reshape(ps,[],1));
p = reshape(p,5,2);
pamptbl = table(comp, ps(:,1), ps(:,2), 'VariableNames', varnames);

% Ttest latencies
[~,p0] = ttest(P0t - P0p);
[~,p1] = ttest(P1t - P1p);
[~,n1] = ttest(N1t - N1p);
[~,p2] = ttest(P2t - P2p);
[~,n2] = ttest(N2t - N2p);
varnames = {'Component', 'fpz t-test p-value','fz t-test p-value'};
comp = ["P0", "P1", "N1", "P2", "N2"]';
ps = [p0; p1; n1; p2; n2];
[~,~,~,p]=fdr_bh(reshape(ps,[],1));
p = reshape(p,5,2);
plattbl = table(comp, ps(:,1), ps(:,2), 'VariableNames', varnames);


%% Save

% Directory
svpath =  '/Applications/Toolbox/MATLAB/Figures/';

% Save
print(f1,fullfile(svpath,'Peak v Trough'),'-dsvg');
print(f1,fullfile(svpath,'Peak v Trough'),'-dpng');
print(f2,fullfile(svpath,'Random'),'-dsvg');
print(f2,fullfile(svpath,'Random'),'-dpng');
print(f3,fullfile(svpath,'ERP Latency'),'-dsvg','-painters');
print(f3,fullfile(svpath,'ERP Latency'),'-dpng','-painters');
print(f4,fullfile(svpath,'ERP Amplitude'),'-dsvg','-painters');
print(f4,fullfile(svpath,'ERP Amplitude'),'-dpng','-painters');
writecell(lattbl,fullfile(svpath,'ANOVA Latency Comparison Table.csv'));
writecell(amptbl,fullfile(svpath,'ANOVA Amplitude Comparison Table.csv'));
writetable(plattbl,fullfile(svpath,'Latency t-test Table.csv'));
writetable(pamptbl,fullfile(svpath,'Amplitude t-test Table.csv'));
