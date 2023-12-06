%% Housekeeping
clear;
close all;
clc;

%% Load Data

% Select subjects to exclude from group analysis
exSub = {'Sub_151','Sub_152','126'}; %,'Sub_102','Sub_104','Sub_126','Sub_138'};

%% Analysis parameters

% Define pre-processing for raw EEG data
params = echtparams('filterImplementation',0,'fslide',1,'tWin',[-0.5 0.5]);

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
    params.causal = 1;
    params.filter = 'echt';

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

    % Store sampling rate
    fs = ERP.fs(1);

    % Reshape ERP
    erp = permute(ERP.erp,[1 3 2]);

    % Separate by phase (& loop through electrodes for fft)
    for ph = 1:2
        pIdx{ph} = intersect(ERP.goodTrials,find(round(ERP.echtPhaseLabel)==ERP.setPhases(ph)));

        % Baseline Period
        idxBL = find(t < 0 & t >-500);
        idxPS = find(t > 0 & t < 500);

        % Get prestimulus power spectrums
        [pre,f]  = pspectrum(erp(idxBL,:,ph),fs);
        [post,~] = pspectrum(erp(idxPS,:,ph),fs);

        % Detrend pre with 3rd order polynomial
        p = polyfit(f,10*log10(median(pre,2)),3); % 3rd-order fit
        pfit = polyval(p,f);
        pre(:,:,ph) = 10*log10(pre)-pfit;

        % Detrend post with 3rd order polynomial
        p = polyfit(f,10*log10(median(post,2)),3); % 3rd-order fit
        pfit = polyval(p,f);
        post(:,:,ph) = 10*log10(post)-pfit;
    end

    % Get mean peak frequency in the alpha range
    freqs = find(f > 8 & f < 14);

    % Loop through electrodes
    for ph = 1:2
        
        % Get average spectrums
        prepeak = median(pre(:,pIdx{1},ph),2);
        pretrough = median(pre(:,pIdx{2},ph),2);
        postpeak = median(post(:,pIdx{1},ph),2);
        posttrough = median(post(:,pIdx{2},ph),2);

        % Find peak alpha in this range trough condition
        [prepk,alphapre]   = findpeaks(pretrough(freqs),f(freqs));
        [postpk,alphapost] = findpeaks(posttrough(freqs),f(freqs));

        % Get frequency greatest power
        [~,idx] = max(prepk);
        prestimtr(k,ph) = alphapre(idx);
        [~,idx] = max(postpk);
        poststimtr(k,ph) = alphapost(idx);

        % Find peak alpha in this range peak condition
        [prepk,alphapre]   = findpeaks(prepeak(freqs),f(freqs));
        [postpk,alphapost] = findpeaks(postpeak(freqs),f(freqs));

        % Get frequency greatest power
        [~,idx] = max(prepk);
        prestimpk(k,ph) = alphapre(idx);
        [~,idx] = max(postpk);
        poststimpk(k,ph) = alphapost(idx);

    end

    % Get frequency sliding
    tslide(k,:,:) = squeeze(mean(ERP.slide(:,:,pIdx{1}),3));
    pslide(k,:,:) = squeeze(mean(ERP.slide(:,:,pIdx{2}),3));

    % Get average phase
    tphase(k,:,:) = angle(squeeze(mean(ERP.phase(:,:,pIdx{1}),3)));
    pphase(k,:,:) = angle(squeeze(mean(ERP.phase(:,:,pIdx{2}),3)));

    % Get Phase Lag 
    lag(k,:,:) = unwrap(tphase(k,:,:)) - unwrap(pphase(k,:,:));
    lag(k,:,:) = lag(k,:,:) - mean(lag(k,idxBL,:),2);

end

%% Post process

% Get change in freqeuency
tslid = tslide - mean(tslide,2);
pslid = pslide - mean(pslide,2);

% Get time for later fill use
time = [t ; flipud(t)];

% Phase means
tphasemu = squeeze(angle(sum(exp(1i*tphase))));
pphasemu = squeeze(angle(sum(exp(1i*pphase))));

% Lag, means, and standard errors
mulag = squeeze(mean(lag));
selag = squeeze(std(lag))./sqrt(20);
selag = [mulag + selag ; flipud(mulag - selag)];

% Mean and SEM - note: SEM are adapted for future use of fill functions
tmu = movmean(squeeze(mean(tslid)),25); % Trough frequency sliding
tse = squeeze(std(tslid))./sqrt(n);
tse = cat(1, tmu + tse, flip(tmu - tse, 1));

pmu = movmean(squeeze(mean(pslid)),25); %  Peak frequency sliding
pse = squeeze(std(pslid))./sqrt(n);
pse = cat(1, pmu + pse, flip(pmu - pse,1));

% Get signifigance
p = npperm(tslid, pslid, 1000,'cluster',0,0);
[~,plag] = rayleigh(lag,1); plag = squeeze(plag);
[~,CV] = benhoch(plag(:,1));
factor = 0.05./CV;
plag = plag * factor;

% Get Area under the curve
tpreauc = squeeze(trapz(t(idxBL),tslid(:,idxBL,:),2));
ppreauc = squeeze(trapz(p(idxBL),pslid(:,idxBL,:),2));
tpostauc = squeeze(trapz(t(idxPS),tslid(:,idxPS,:),2));
ppostauc = squeeze(trapz(p(idxPS),pslid(:,idxPS,:),2));

%% Plots

% Grab the colors figures
co2 = [0.0000 0.4470 0.7410
       0.8500 0.10 0.0980];

% Phase figure
f0 = figure;
sgtitle('Phases by Condition')

subplot(2,2,1)
plot(t, tphasemu(:,1),'color', co2(1,:), 'linewidth',2)
hold on
plot(t, pphasemu(:,1),'color', co2(2,:), 'linewidth',2)
title("Phase FPZ")
xlim([-250 500])
ylabel("Phase (degrees)")
pbaspect([3 1 1])
yticks([-pi pi])
yticklabels({'-180', '180'})
box off
grid on

subplot(2,2,2)
plot(t, tphasemu(:,2),'color', co2(1,:), 'linewidth',2)
hold on
plot(t, pphasemu(:,2),'color', co2(2,:), 'linewidth',2)
title("Phase FZ")
xlim([-250 500])
pbaspect([3 1 1])
yticks([-pi pi])
yticklabels({'-180', '180'})
box off
grid on

subplot(2,2,3)
plot(t, mulag(:,1),'color', [0.5 0.5 0.5], 'linewidth',2)
hold on
fill([t ; flipud(t)],selag(:,1),[0.5 0.5 0.5], 'facealpha', 0.3, 'edgecolor','none');
title("Phase Lag FPZ")
xlim([-250 500])
ylim([-0.5 4])
ylabel("Phase Lag (degrees)")
pbaspect([3 1 1])
yticks([0 pi/2 pi])
yticklabels({'180', '90','0'})
box off
grid on

subplot(2,2,4)
plot(t, mulag(:,2),'color', [0.5 0.5 0.5], 'linewidth',2)
hold on
fill([t ; flipud(t)],selag(:,2),[0.5 0.5 0.5], 'facealpha', 0.3, 'edgecolor','none');
title("Phase Lag FZ")
xlim([-250 500])
ylim([-0.5 4])
pbaspect([3 1 1])
yticks([0 pi/2 pi])
yticklabels({'180', '90','0'})
box off
grid on


% Start figure of frequency sliding
f1 = figure;

% FPZ
subplot(2,2,1)
plot(t, tmu(:,1),'color',co2(1,:),'linewidth',2)
hold on
fill(time, tse(:,1),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, pmu(:,1),'color',co2(2,:),'linewidth',2)
fill(time, pse(:,1),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlabel('Time (ms)')
ylabel('Frequency Change (Hz)')
pbaspect([3 1 1])
ylim([-0.4 0.4])
xlim([-250 500])
box off
title('Frequency Sliding FPZ')

% FZ
subplot(2,2,2)
plot(t, tmu(:,2),'color',co2(1,:),'linewidth',2)
hold on
fill(time, tse(:,2),co2(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
plot(t, pmu(:,2),'color',co2(2,:),'linewidth',2)
fill(time, pse(:,2),co2(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
xlabel('Time (ms)')
ylabel('Frequency Change (Hz)')
pbaspect([3 1 1])
ylim([-0.4 0.4])
xlim([-250 500])
box off
title('Frequency Sliding FZ')

% Fpz signifigance
subplot(2,2,3)
plot(t,p(:,1),'k-')
hold on
fill(time, [p(:,1); ones(1,length(p(:,1)))'],'k')
xlim([-250 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
pbaspect([3 1 1])
box off

% Fz signifigance
subplot(2,2,4)
plot(t,p(:,2),'k-')
hold on
fill(time, [p(:,2); ones(1,length(p(:,2)))'],'k')
xlim([-250 max(t)])
ylim([0 0.05])
set(gca,'ydir','reverse')
xlabel('Time (ms)')
pbaspect([3 1 1])
box off

min(p(find(t > 0 & t < 100),:))
min(p(find(t > 100 & t < 200),:))

% Save
svpath = '/Applications/Toolbox/MATLAB';
print(f0,fullfile(svpath,'Phase Figures.svg'),'-dsvg','-painters');
print(f0,fullfile(svpath,'Phase Figures.png'),'-dpng','-painters');
print(f1,fullfile(svpath,'Frequency Sliding.svg'),'-dsvg','-painters');
print(f1,fullfile(svpath,'Frequency Sliding.png'),'-dpng','-painters');