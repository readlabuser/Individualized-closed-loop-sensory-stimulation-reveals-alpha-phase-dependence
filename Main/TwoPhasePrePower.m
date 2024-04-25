%% Housekeeping
% Clear all figures, data, and set a timer.

close all
clear
clc

% Set timer
h = tic;

%% Prep Work

% This is where we load all of our participant files directories', as well as 
% parameters for processing the output of the ECHT device and exclude participant's 
% with artifacts.

% Load current data table
datafilename = '/Applications/Toolbox/MATLAB/TwoPhaseExp_table.mat';    %% Changed for Windows
load(datafilename);
dataTable = sortrows(dataTable);

% Prep ECHT Parameters - speciify wavelet frequency range, resolution + 
% low cycles creates wide frequency ranges
params = echtparams('filterImplementation',1,'freqs', [2 70], 'steps', 70,...
                    'cycles', [5 12],'tWin', [-0.7 0]);
params.causal = 1;
params.filter = 'echt';

% Load subject directories
pname = genpath('/Applications/Toolbox/SubjectData');  %% Changed for Windows
pname = strsplit(pname,':')'; %% Switched

% Locate Random Phase ERP or TwoPhase datasets
pname1 = pname(contains(pname,'ISI750ms85dB') & ~contains(pname,'(Low'));
pname2 = pname(contains(pname,'TwoPhaseERP85dB') & ~contains(pname,'(Low'));

% Select subjects to exclude from group analysis
exSub = {'Sub_151','Sub_152','126'}; 
if(~isempty(exSub))  
    pname1 = pname1(~contains(pname1,exSub)); % exclude subject(s)
    pname2 = pname2(~contains(pname2,exSub));
end

% Setup pyversion for fooof
% py = pyversion('/Users/tylor/anaconda3/bin/python');

%% Gather Data

% This is where we loop through each subject, and grab their prestimulus power 
% for each condition, showing that they are in fact equal, as well as consistent 
% with prior work (Iemi et al., 2019; Cabral-Calderin et al., 2020).

% If compute calderin ways
calderin = 0;
pspctrm = 0;

for k = 1:(length(pname1))

    % Get subject index
    tmp = strsplit(pname1{k},filesep);
    Sub_ID = tmp{contains(tmp,'Sub_')};
    subIdx = contains(dataTable.SubID,Sub_ID([5:end]));

    % Setup alpha center frequency
    params.alphaCF = dataTable(subIdx,:).IAF(1);

    % Process 2Phase data
    [ERP,params2] = batchERP(pname2{k},params);
    fs = ERP.fs(1);
    t = ERP.t0;

    % Reformat ERP
    erp = permute(ERP.erp,[3 1 2]);

    % Check order of optimal/pessimal phase list in ERP structure against the
    % calculated phases in 'dataTable.' Fix wrap-around error if necessary
    if(dataTable(subIdx,:).OptPhase > dataTable(subIdx,:).PesPhase)
        ERP.setPhases = ERP.setPhases([2 1]);
    end

    % Sort trial types & good trials
    for i = 1:2
        pIdx{i} = intersect(ERP.goodTrials,...
                  find(round(ERP.echtPhaseLabel)==ERP.setPhases(i)));

        % Get trial type
        tmp = squeeze(median(erp(pIdx{i},:,:)));

        % Define hanning window
        win = .54 - .46*cos(2*pi*(0:length(ERP.t0)-1)'/(length(ERP.t0)-1));
        
        % Power spectrum
        frange = [1 50];
        settings = struct('aperiodic_mode', 'knee');

        % if unfooofed
        if calderin
            f = linspace(0,fs,fs)';
            pxx = abs(squeeze(fft(squeeze(win .* tmp),fs)));
            pyy(k,i,:,1:2) = 10*log10(pxx.^2./mean(pxx.^2));
        else
            % Smoothed power spectrum
            if pspctrm
                [pxx, f] = pspectrum(tmp, fs, 'FrequencyLimits', frange);
            else
                f = linspace(0,fs,fs)';
                pxx = abs(squeeze(fft(squeeze(win .* tmp),fs)));
            end

            % each electrode
            for j = 1:2
                mdl = fooof(f, pxx(:,j), frange, settings, 1);

                % Detrend
                pyy(k,i,:,j) = mdl.fooofed_spectrum - mdl.ap_fit;

            end
        end

    end
    
    % Process Random data
    [ERP,params2] = batchERP(pname1{k},params);

    % Reformat ERP
    erp = permute(ERP.erp,[3 1 2]);

    % Get AEP
    tmp = zscore(squeeze(median(erp(ERP.goodTrials,:,:))));

    % Power spectrum
    frange = [1 50];
    settings = struct('aperiodic_mode', 'knee');

    % if unfooofed
    if calderin
        f = linspace(0,fs,fs)';
        pxx = abs(squeeze(fft(squeeze(win .* tmp),fs)));
        pyy(k,3,:,1:2) = 10*log10(pxx.^2./mean(pxx.^2));
    else
        % Smoothed power spectrum
        if pspctrm
            [pxx, f] = pspectrum(tmp, fs, 'FrequencyLimits', frange);
        else
            f = linspace(0,fs,fs)';
            pxx = abs(squeeze(fft(squeeze(win .* tmp),fs)));
        end

        % each electrode
        for j = 1:2
            mdl = fooof(f, pxx(:,j), frange, settings, 1);
    
            % Detrend
            pyy(k,3,:,j) = mdl.fooofed_spectrum - mdl.ap_fit;
    
        end
    end

end
%% Post Process

% Permutation test parameters
iter = 1000;
correction = 'cluster';
circ = 0;
trial_level = 0;

% Permute
trough = squeeze(pyy(:,1,:,:));
peak = squeeze(pyy(:,2,:,:));
random = squeeze(pyy(:,3,:,:));
[ppwr,~,~] = npperm(trough, peak, iter, correction, circ, trial_level);

% Resize frequency vector
if calderin
    f = f';
    ylims = [-10 20];
    ytick = [-10 0 10 20];
    xlims = [0 15];
    ylbl = '(dB)';
else
    f = mdl.freqs;
    ylims = [-0.2 2];
    ytick = [0 0.5 1 1.5 2];
    xlims = [1 50];
    ylbl = '(uV^2/Hz)';
end
freqs = [f , fliplr(f)];

% Get colors
color = [0.0000 0.4470 0.7410
         0.8500 0.1000 0.0980];

% First set
f0 = figure;

% Get mean and standard error
mu = squeeze(mean(pyy(:,:,:,1)));
se = squeeze(std(pyy(:,:,:,1))./sqrt(19));

% Resize standard errors
se = [mu - se, fliplr(mu + se)];

% Plot Trough Spectrum
plot(f, mu(1,:),'color', color(1,:), 'linewidth', 2)
hold on
fill(freqs, se(1,:), color(1,:),'edgecolor','none','facealpha', 0.1)
xlim(xlims)
xlabel('Frequency (Hz)')
ylabel(sprintf('PSD %s',ylbl))
ylim(ylims)
yticks(ytick)
box off
axis square
pbaspect([2 1 1])
hold on

% Plot Peak Spectrum 
plot(f, mu(2,:),'color', color(2,:), 'linewidth', 2)
hold on
fill(freqs, se(2,:), color(2,:),'edgecolor','none','facealpha', 0.1)
xlim(xlims)
xlabel('Frequency (Hz)')
ylabel(sprintf('PSD %s',ylbl))
ylim(ylims)
yticks(ytick)
box off
axis square
pbaspect([2 1 1])
hold on

grey  = [0.5 0.5 0.5];
[x, y] = clusts4plot(ppwr(:,1), 0.05, ylims);
if ~isempty(x)
    for i = 1:size(x,1)
        patch(f(x(i,:)), y, grey,'LineStyle','none','FaceAlpha',0.1);
        hold on
    end
end


% Plot Random Spectrum 
f2 = figure;
plot(f, mu(3,:),'color', [0.5 0.5 0.5], 'linewidth', 2)
hold on
fill(freqs, se(3,:), [0.5 0.5 0.5],'edgecolor','none','facealpha', 0.1)
xlim(xlims)
title('Pre-stimulus Power - FPZ (-700 ms)')
xlabel('Frequency (Hz)')
ylabel(sprintf('PSD %s',ylbl))
ylim(ylims)
yticks(ytick)
box off
axis square
pbaspect([2 1 1])
hold on



% Second set
f1 = figure;

% Get mean and standard error
mu = squeeze(mean(pyy(:,:,:,2)));
se = squeeze(std(pyy(:,:,:,2))./sqrt(20));

% Resize standard errors
se = [mu - se, fliplr(mu + se)];

% Plot Trough Spectrum
plot(f, mu(1,:),'color', color(1,:), 'linewidth', 2)
hold on
fill(freqs, se(1,:), color(1,:),'edgecolor','none','facealpha', 0.1)
xlim(xlims)
xlabel('Frequency (Hz)')
ylabel(sprintf('PSD %s',ylbl))
ylim(ylims)
yticks(ytick)
box off
axis square
pbaspect([2 1 1])
hold on

% Plot Peak Spectrum 
plot(f, mu(2,:),'color', color(2,:), 'linewidth', 2)
hold on
fill(freqs, se(2,:), color(2,:),'edgecolor','none','facealpha', 0.1)
xlim(xlims)
xlabel('Frequency (Hz)')
ylabel(sprintf('PSD %s',ylbl))
ylim(ylims)
yticks(ytick)
box off
axis square
pbaspect([2 1 1])
hold on

grey  = [0.5 0.5 0.5];
[x, y] = clusts4plot(ppwr(:,2), 0.05, ylims);
if ~isempty(x)
for i = 1:size(x,1)
    patch(f(x(i,:)), y, grey,'LineStyle','none','FaceAlpha',0.1);
    hold on
end
end

% Plot Random Spectrum 
f3 = figure;
plot(f, mu(3,:),'color', [0.5 0.5 0.5], 'linewidth', 2)
hold on
fill(freqs, se(3,:), [0.5 0.5 0.5],'edgecolor','none','facealpha', 0.1)
xlim(xlims)
yticks(ytick)
title('Pre-stimulus Power - FZ (-700 ms)')
xlabel('Frequency (Hz)')
ylabel(sprintf('PSD %s',ylbl))
ylim(ylims)
yticks(ytick)
box off
axis square
pbaspect([2 1 1])
hold on


% Save Figure
svpath = '/Applications/Toolbox/MATLAB/Figures';
print(f0, fullfile(svpath, 'non-Fooofed TwoPhase Pre-Stimulus Spectrogram FPZ  (-700 ms)'),'-dsvg');
print(f0, fullfile(svpath, 'non-Fooofed TwoPhase Pre-Stimulus Spectrogram FPZ (-700 ms)'),'-dpng');
print(f1, fullfile(svpath, 'non-Fooofed TwoPhase Pre-Stimulus Spectrogram Fz (-700 ms)'),'-dsvg');
print(f1, fullfile(svpath, 'non-Fooofed TwoPhase Pre-Stimulus Spectrogram Fz (-700 ms)'),'-dpng');

print(f2, fullfile(svpath, 'non-Fooofed Random Pre-Stimulus Spectrogram FPZ  (-700 ms)'),'-dsvg');
print(f2, fullfile(svpath, 'non-Fooofed Random Pre-Stimulus Spectrogram FPZ (-700 ms)'),'-dpng');
print(f3, fullfile(svpath, 'non-Fooofed Random Pre-Stimulus Spectrogram Fz (-700 ms)'),'-dsvg');
print(f3, fullfile(svpath, 'non-Fooofed Random Pre-Stimulus Spectrogram Fz (-700 ms)'),'-dpng');


% Timing is everything
sprintf('Time to run: %d', toc(h));