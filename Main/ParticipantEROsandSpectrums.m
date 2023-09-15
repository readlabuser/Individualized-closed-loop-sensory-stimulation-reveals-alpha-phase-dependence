%% Housekeeping
% Clear all figures, data, and set a timer.
close all
clear
clc

% Set timer
h = tic;

% Setup pyversion for fooof
py = pyenv;

%% Prep Work
% This is where we load all of our participant files directories', as well as 
% parameters for processing the output of the ECHT device and exclude participant's 
% with artifacts.

% Load current data table
datafilename = '/Applitcations/Toolbox/MATLAB/TwoPhaseExp_table.mat';    %% Changed for Windows
load(datafilename);
dataTable = sortrows(dataTable);

% Prep ECHT Parameters - speciify wavelet frequency range, resolution + 
% low cycles creates wide frequency ranges
params = echtparams();

% Load subject directories
pname = genpath('/Applications/Toolbox/SubjectData');  %% Changed for Windows
pname = strsplit(pname,':')'; %% Switched

% Select subjects to exclude from group analysis
exSub = {'151','152', '126'}; 
pname = pname(~contains(pname,exSub));

% Locate Random Phase ERP or TwoPhase datasets
pname0 =  pname(contains(pname,'IAFest') & ~contains(pname,'(Low'));
pname1 = pname(contains(pname,'ISI750ms85dB') & ~contains(pname,'(Low'));
pname2 = pname(contains(pname,'TwoPhaseERP') & ~contains(pname,'(Low'));

%% Gather data
% This loops through every participant, and saves int he proper directory figures 
% for each participant.  Each directory contains a figure of their detrended spectrogram, 
% their random ERP P1 aligned with individualized theoretical alpha, as well as 
% a set of observed evoked responses for for all conditions.


% Whats my sample?
n = length(pname1);

for k = 1:n

    % Get subject index
    tmp = strsplit(pname1{k},filesep);
    Sub_ID = tmp{contains(tmp,'Sub_')}
    subIdx = contains(dataTable.SubID,Sub_ID([5:end]));
    tbl(k) = dataTable(subIdx,:).IAF(1);
    
    % Process Random data
    [ERP,params] = batchERP(pname1{k},params);
    fs = ERP.fs(1);
    t = ERP.t0;

    % Get baseline period
    idxBL = find(t < -50 & t > -250);

    % Reformat ERP
    erp = squeeze(mean(permute(ERP.erp,[3 1 2])));

    % Baseline Normalize
    erp = erp - mean(erp(idxBL,:));

    % Get standard errors
    se = squeeze(std(permute(ERP.erp,[3 1 2]))./sqrt(20));

    % Process TwoPhase data
    [ERP,params] = batchERP2(pname2{k},params);

    % Reformat Morlet Transform
    tmp = permute(ERP.erp,[3 1 2]);

    % Check order of optimal/pessimal phase list in ERP structure against the
    % calculated phases in 'dataTable.' Fix wrap-around error if necessary
    if(dataTable(subIdx,:).OptPhase > dataTable(subIdx,:).PesPhase)
        ERP.setPhases = ERP.setPhases([2 1]);
    end

    % Sort trial types & good trials
    for i = 1:2
        pIdx{i} = intersect(ERP.goodTrials,...
            find(round(ERP.echtPhaseLabel)==ERP.setPhases(i)));

        % Get ERP
        tpmu(i,:,:) = squeeze(mean(tmp(pIdx{i},:,:)));
        tpmu(i,:,:) = tpmu(i,:,:) - mean(tpmu(i,idxBL,:));

        % Get standard errors
        tpse(i,:,:) = squeeze(std(permute(ERP.erp,[3 1 2]))./sqrt(20));

        % Get index of alpha center frequency
        [~,iaf] = min(abs(dataTable(subIdx,:).IAF(1) - ERP.f));

        % Set alpha
        aband(i,k,:,:) = squeeze(abs(mean(ERP.cwt(:,:,iaf, pIdx{i}), 4)));

        % Baseline
        aband(i,k,:,:) = aband(i,k,:,:) - mean(aband(i,k,idxBL,:), 3);


    end

    % Set CI time
    time = [ERP.t0 ; flipud(ERP.t0)]; 

    % Estimate intrinsic alpha activity from IAF alphaCF
    alphaSpline = normpdf(t,0,200);
    alphaSpline = alphaSpline/alphaSpline(t==0);

    % Get alpha center frequency
    alphaCF = dataTable(subIdx,:).IAF(1);
    P1lat = dataTable(subIdx,:).P1lat(1);

    % Get phase estimates
    optPhase = -(2*pi*(P1lat*1e-3)*alphaCF);
    pesPhase = optPhase-pi;

    % Get baseline period
    idxBL = find(t < 100 & t > -250);

    % Modeled alpha matrix
    modA(:,1) = max(abs(erp(idxBL,1)))*alphaSpline.*cos(2*pi*alphaCF*t*1e-3+optPhase);
    modA(:,2) = max(abs(erp(idxBL,1)))*alphaSpline.*cos(2*pi*alphaCF*t*1e-3+pesPhase);

    % Set trial type
    trial = {'Trough Alignment of P1', 'Peak Alignment of P1'};

    % Setup figure
    f0 = figure;
    for i = 1:2

        % Plot]
        subplot(2,1,i)
        plot(t, erp(:,1), 'color', [0.5 0.5 0.5],'linewidth',2)
        hold on
        plot(t,modA(:,i),'-','Color',0.7*ones(1,3),'LineWidth',1);
        box off
        grid on
        xlabel('Time (ms)')
        ylabel('Magnitude (uV)')
        ylim([-4 4])
        yticks([-4 0 4])
        xlim([min(t) max(t)])
        xticks([-200 -100 0 100 200 300 400 500])
        pbaspect([3 1 1])
        title(sprintf(trial{i}))


    end

    % Make Individual Subject Path
    svpath = fullfile('C:\Toolbox\MATLAB\Figures\All ERPs and EROs Not 126', Sub_ID, '\');
    if ~isdir(svpath)
        mkdir(svpath)
    end

    % Print to path
    % print(f0, fullfile(svpath, 'Random ERPs & EROs.svg'),'-dsvg')
    % print(f0, fullfile(svpath, 'Random ERPs & EROs.png'),'-dpng')

    % Get colors
    color = [0.0000 0.4470 0.7410
        0.8500 0.10 0.0980];

    % setup random
    time = [t ; flipud(t)];
    se = [erp + se ; flipud(erp - se)];

    f1(1) = figure;
    plot(t, erp(:,1), 'color', [0.3 0.3 0.3]*1,'linewidth',2)
    hold on
    fill(time,se(:,1), [0.3 0.3 0.3]*1,'edgecolor','none','facealpha', 0.1)
    box off
    xlabel('Time (ms)')
    ylabel('Magnitude (uV)')
    title('Random ERP FPZ')
    % ylim([-8 6])
    xlim([min(t) max(t)])
    xticks([-200 -100 0 100 200 300 400 500])
    pbaspect([2 1 1])

    f1(2) = figure;
    plot(t, erp(:,2), 'color', [0.3 0.3 0.3]*2,'linewidth',2)
    hold on
    fill(time,se(:,2), [0.3 0.3 0.3]*2,'edgecolor','none','facealpha', 0.1)
    box off
    xlabel('Time (ms)')
    ylabel('Magnitude (uV)')
    title('Random ERP FZ')
    % ylim([-8 6])
    xlim([min(t) max(t)])
    xticks([-200 -100 0 100 200 300 400 500])
    pbaspect([2 1 1])

    % Print to Path
    % print(f1(1), fullfile(svpath, 'Random FPZ ERPs'),'-dsvg')
    % print(f1(1), fullfile(svpath, 'Random FPZ ERPs'),'-dpng')
    % print(f1(2), fullfile(svpath, 'Random FZ ERPs'),'-dsvg')
    % print(f1(2), fullfile(svpath, 'Random FZ ERPs'),'-dpng')

    
    % Setup twophase
    mu = squeeze(tpmu(1,:,:));
    se2 = squeeze(tpse(1,:,:));
    se2 = [mu + se2 ; flipud(mu - se2)];

    % Plot
    f1(3) = figure;
    plot(t, mu(:,1), 'color', color(1,:),'linewidth',2)
    hold on
    fill(time,se2(:,1),color(1,:),'edgecolor','none','facealpha', 0.1)
    box off
    xlabel('Time (ms)')
    ylabel('Magnitude (uV)')
    title('TwoPhase Fpz')
    % ylim([-15 15])
    yticks([-15 -10 -5 0 5 10 15])
    xlim([min(t) max(t)])
    xticks([-200 -100 0 100 200 300 400 500])
    pbaspect([2 1 1])

    mu = squeeze(tpmu(2,:,:));
    se2 = squeeze(tpse(2,:,:));
    se2 = [mu + se2 ; flipud(mu - se2)];

    plot( t, mu(:,1), 'color', color(2,:),'linewidth',2)
    hold on
    fill(time,se2(:,1),color(2,:),'edgecolor','none','facealpha', 0.1)
    box off
    xlabel('Time (ms)')
    ylabel('Magnitude (uV)')
    title('TwoPhase Fpz')
    % ylim([-15 15])
    yticks([-15 -10 -5 0 5 10 15])
    xlim([min(t) max(t)])
    xticks([-200 -100 0 100 200 300 400 500])
    pbaspect([2 1 1])


    % Setup twophase
    mu = squeeze(tpmu(1,:,:));
    se2 = squeeze(tpse(1,:,:));
    se2 = [mu + se2 ; flipud(mu - se2)];

    
    f1(4) = figure;
    plot(t, mu(:,2), 'color', color(1,:),'linewidth',2)
    hold on
    fill(time,se2(:,2),color(1,:),'edgecolor','none','facealpha', 0.1)
    box off
    xlabel('Time (ms)')
    ylabel('Magnitude (uV)')
    title('TwoPhase Fz')
    % ylim([-15 15])
    yticks([-15 -10 -5 0 5 10 15])
    xlim([min(t) max(t)])
    xticks([-200 -100 0 100 200 300 400 500])
    pbaspect([2 1 1])

    mu = squeeze(tpmu(2,:,:));
    se2 = squeeze(tpse(2,:,:));
    se2 = [mu + se2 ; flipud(mu - se2)];

    plot( t, mu(:,1), 'color', color(2,:),'linewidth',2)
    hold on
    fill(time,se2(:,1),color(2,:),'edgecolor','none','facealpha', 0.1)
    box off
    xlabel('Time (ms)')
    ylabel('Magnitude (uV)')
    title('TwoPhase Fpz')
    % ylim([-15 15])
    yticks([-15 -10 -5 0 5 10 15])
    xlim([min(t) max(t)])
    xticks([-200 -100 0 100 200 300 400 500])
    pbaspect([2 1 1])

    % Print to Path
    % print(f1(3), fullfile(svpath, 'Observed FPZ ERPs'),'-dsvg')
    % print(f1(3), fullfile(svpath, 'Observed FPZ ERPs'),'-dpng')
    % print(f1(4), fullfile(svpath, 'Observed FZ ERPs'),'-dsvg')
    % print(f1(4), fullfile(svpath, 'Observed FZ ERPs'),'-dpng')


    % Get Iaf Estimate data
    [ERP,params] = batchERP(pname0{k},params);

    % Get MT Params
    N = 10; %12;
    df = 0.5; %0.5;
    dt = 0.15;
    TW = (N*df)/2; % time half-bandwidth
    L = floor(2*TW)-1; % number of tapers
    MTparams.tapers = [TW L]; %[TW L];
    MTparams.Fs = fs; % sampling frequency (Hz)
    MTparams.fpass = [1 30]; % frequencies represented in Hz [low high]

    % Reformat recording
    EEG = [];
    for i = 1:length(ERP.eeg)
        EEG = cat(1,EEG,ERP.eeg{i});
    end

    % MT spectrogram
    [Pxx,txx,fxx] = mtspecgramc(EEG,[N dt],MTparams);

    % Save multiparticipant power
    power(k,:,:,:) = 10*log10(Pxx(txx <= 120,:,:));

    % Get channels
    ch = {'FPZ', 'FZ'};

    for i = 1:2

        % Power law fit: fit to median of spectrogram
        p = polyfit(fxx,10*log10(median(Pxx(:,:,i))),3); % 3rd-order fit
        pfit = polyval(p,fxx);
        pyy = 10*log10(Pxx(:,:,i))-pfit;

        % Get alphaCF from 3rd order polynomial fit
        aRng = fxx>=7.5 & fxx<=15;
        [alphaPks,alphaFr,wA,pA] = findpeaks(median(pyy(:,aRng)),fxx(aRng),'MinPeakHeight',0.5,'MinPeakProminence',1);

        % Select strongest peak within the defined alpha range
        [~,isMaxAlpha] = max(alphaPks);
        acfpoly(k,i) = alphaFr(isMaxAlpha);

        % Save detrended power
        dpower(k,:,:,i) = pyy(txx <= 120,:);

        % Plot MT spectrogram
        f2(i) = figure;
        imagesc(txx,fxx,10*log10(Pxx(:,:,i))'); %pyy'
        xticks([0 20 40 60 80 100 120 140])
        colormap jet
        set(gca,'Ydir','normal')
        box off
        clim([-20 20])
        colorbar
        ylabel('Frequency (Hz)')
        xlabel('Time (s)')
        ylim([5 30])
        yticks([5 10 15 20 25 30])
        xlim([20 120])
        pbaspect([3 1 1])

        % Setup participant average spectrum
        polyspctrm(k,:,i) = squeeze(median(pyy));

        % Print to Path
        % print(f2(i), fullfile(svpath, sprintf('MT Spectrograms %s', ch{i})),'-dsvg')
        % print(f2(i), fullfile(svpath, sprintf('MT Spectrograms %s', ch{i})),'-dpng')

        % Plot polynomial detrended spectrogram
        f3(i) = figure;
        subplot(2,1,1)
        plot(fxx,median(pyy),'color',color(1,:),'linewidth',.5)
        hold on
        fill([fxx, fliplr(fxx)],...
            [mean(pyy) + abs(std(pyy)./sqrt(20)) , fliplr(mean(pyy) - abs(std(pyy)./sqrt(20)))],...
            color(1,:),'edgecolor','none','facealpha', 0.1)
        xlabel('Frequency (Hz)')
        ylabel('Power (dB)')
        xticks([0 10 20 30 40])
        xlim([0 30])
        ylim([-5 6])
        axis square
        grid on
        title('Polynomial Detrend')
        box off

        % Dot for IAF
        hold on
        scatter(alphaFr(isMaxAlpha), alphaPks(isMaxAlpha), 25, color(2,:),'filled')

        % setup
        N = size(EEG, 1); pnts = N;
        hz = linspace(0, params.fs/2, floor(N/2)+1);
        fc = fft(EEG(:,i)); 
        fc = fc(1:length(hz));
        pwr = abs(fc);

        % fooof spectrum
        frange = [1 30]; settings = struct();
        mdl = fooof(hz, pwr, frange, settings, 1);

        % Rename variables and subtract aperiodic fit
        freqs = mdl.freqs;
        rcrdd = mdl.power_spectrum;
        crrctd = mdl.fooofed_spectrum - mdl.ap_fit;

        % Setup participant average fooof spectrum
        fooofspctrm(k,:,i) = crrctd(1:4000);

        % Get alphaCF from 3rd order polynomial fit
        % aRng = freqs(freqs>=7.5 & freqs<=15);
        % [maxalpha,isMaxAlpha] = max(crrctd(freqs>=7.5 & freqs<=15));
        % acffooof(k,i) = aRng(isMaxAlpha);


        % Get alphaCF from 3rd order polynomial fit
        aRng = freqs>=7.5 & freqs<=15;
        [alphaPks,alphaFr,wA,pA] = findpeaks(crrctd(aRng),freqs(aRng));

        % Select strongest peak within the defined alpha range
        [~,isMaxAlpha] = max(alphaPks);
        acffooof(k,i) = alphaFr(isMaxAlpha);

        % Plot fooof detrend
        subplot(2,1,2)
        plot(freqs,crrctd,'color',color(1,:),'linewidth',2)
        xlabel('Frequency (Hz)')
        ylabel('Power (dB)')
        title('Fooof Detrend')
        xticks([0 10 20 30 40])
        grid on
        axis square
        xlim([0 30])
        box off
        hold on

        % Dot for IAF
        hold on
        scatter(alphaFr(isMaxAlpha), alphaPks(isMaxAlpha), 25, color(2,:),'filled')

        % Print to Path
        % print(f3(i), fullfile(svpath, sprintf('Detrended Spectrums %s', ch{i})),'-dsvg')
        % print(f3(i), fullfile(svpath, sprintf('Detrended Spectrums %s', ch{i})),'-dpng')


        % Welch Spectrogram
        N = size(EEG,1);
        [wpxx,wfxx] = pwelch(EEG,[],[],[1:fs/N:30],fs);
        p = polyfit(wfxx,10*log10(wpxx(:,i)),3); % 3rd-order fit
        pfit = polyval(p,wfxx);
        wpyy = 10*log10(wpxx) - pfit;
        wpyy = wpyy(1:4000); wfxx = wfxx(1:4000);

        % Get alphaCF from 3rd order polynomial fit
        % aRng = wfxx(wfxx>=7.5 & wfxx<=15);
        % [maxalpha,isMaxAlpha] = max(wpyy(wfxx>=7.5 & wfxx<=14));
        % acfwelch(k,i) = aRng(isMaxAlpha);

        % Get alphaCF from 3rd order polynomial fit
        aRng = wfxx>=7.5 & wfxx<=15;
        [alphaPks,alphaFr,wA,pA] = findpeaks(wpyy(aRng),wfxx(aRng));

        % Select strongest peak within the defined alpha range
        [~,isMaxAlpha] = max(alphaPks);
        acfwelch(k,i) = alphaFr(isMaxAlpha);

        % Save welch spectrum
        welchspctrm(k,:,i) = wpyy(1:4000);
        wfxx = wfxx(1:4000);

        % Plot
        wel(i) = figure;
        plot(wfxx, wpyy, 'color', color(1,:))
        xlim([0 30])
        title("Power Spectrum from Welch")
        xlabel('Frequency (Hz)')
        ylabel('Power (dB)')
        grid on
        axis square
        box off
        xticks([0 10 20 30 40])

        % Print to Path
        print(wel(i), fullfile(svpath, sprintf('Detrended Welch Spectrums %s', ch{i})),'-dsvg')
        print(wel(i), fullfile(svpath, sprintf('Detrended Welch Spectrums %s', ch{i})),'-dpng')

    end

    % Get alpha band-delimit for participant
    [tf, fr] = wavecwt(EEG, 501, [1 30], [4 8], 400);
    tf = real(tf);
    [~,idx] = min(abs(dataTable(subIdx,:).IAF(1) - fr));
    alpha = tf(:,:,idx);

    f4 = figure;
    subplot(2,1,1)
    plot(ERP.time{1}, alpha(1:length(ERP.time{1}),1),'color', [0.7 0.7 0.7], 'linewidth',2)
    pbaspect([3 1 1])
    xticks([20 20.5 21 21.5])
    xticklabels({'0', '500', '1000', '1500'})
    ylim([-10 10])
    yticks([-10 0 10])
    xlim([20 21.5])
    title('Alpha Fpz')
    ylabel('Amplitude (uV)')
    grid on
    box off

    subplot(2,1,2)
    plot(ERP.time{1}, alpha(1:length(ERP.time{1}),2),'color', [0.7 0.7 0.7], 'linewidth',2)
    xticks([20 20.5 21 21.5])
    xticklabels({'0', '500', '1000', '1500'})
    pbaspect([3 1 1])
    xlim([20 21.5])
    ylim([-10 10])
    yticks([-10 0 10])
    title('Alpha Fz')
    ylabel('Amplitude (uV)')
    grid on
    box off
    xlabel('Time (ms)')

    % Print
    print(f4, fullfile(svpath, 'Alpha IAF Estimate'),'-dpng')   
    print(f4, fullfile(svpath, 'Alpha IAF Estimate'),'-dsvg')

    % Close figures
    close all

end

%% Post Process

% Grab means
mufooof = squeeze(nanmean(fooofspctrm));
mupoly = squeeze(nanmean(polyspctrm));
muwelch = squeeze(nanmean(welchspctrm));

% Grab Standard Errors
sefooof = squeeze(nanstd(fooofspctrm))./sqrt(n);
sepoly =  abs(squeeze(nanstd(polyspctrm))./sqrt(n));
sewelch = squeeze(nanstd(welchspctrm))./sqrt(n);
sefooof = [mufooof + sefooof; flipud(mufooof - sefooof)];
sepoly =  [mupoly + sepoly; flipud(mupoly - sepoly)];
sewelch = [muwelch + sewelch; flipud(muwelch - sewelch)];

%% Figures
% Compare IAF Estimates
f = figure;
for i = 1:2
    
    % Get model
    lm = fitlm(acfpoly(:,i), acffooof(:,i));
    
    % Plot
    subplot(1,2,i)
    scatter(acfpoly(:,i), acffooof(:,i), 70, color(1,:), 'filled')
    xlabel('Polynomial IAF')
    ylabel('Fooof IAF')
    grid on
    axis square
    box off
    hold on
    xlim([7 15])
    ylim([7 15])

    % Plot Identity line
    plot([7 15], [7 15], 'color', color(1,:), 'linewidth',2)
    title(sprintf('Fooof vs. Polynomial Detrend %s', ch{i}));

    % Add Model text
    text(11,9,sprintf('R^2 = %.3f \np-value = %.4f',...
        lm.Rsquared.Adjusted, lm.ModelFitVsNullModel.Pvalue));
    
end

% Save
svpath = fullfile('C:\Toolbox\MATLAB\Figures\All ERPs and EROs\');
print(f, fullfile(svpath, 'Fooof vs poly'), '-dsvg')
print(f, fullfile(svpath, 'Fooof vs poly'), '-dpng')

% Save
save('Polynomial or Fooof.mat', 'acfpoly','acffooof','lm');

f0 = figure;
for i = 1:2

    % Plot
    subplot(1,2,i)
    plot(fxx, mupoly(:,i),'color', color(1,:), 'linewidth',.5)
    hold on
    fill([fxx , fliplr(fxx)], sepoly(:,i), color(1,:), ...
        'edgecolor','none','facealpha', 0.1);
    title(sprintf('Polynomial Detrended PSD %s', ch{i}))
    xticks([0 10 20 30 40])
    axis square
    box off
    grid on
    xlim([0 30])
    ylim([-4 6])


    % Get alphaCF from 3rd order polynomial fit
    aRng = fxx(fxx>=7.5 & fxx<=15);
    [maxalpha,isMaxAlpha] = max(mupoly(fxx>=7.5 & fxx<=15,i));

    % Dot for IAF
    hold on
    aRng(isMaxAlpha)
    scatter(aRng(isMaxAlpha), maxalpha, 25, color(2,:),'filled')

end
% 
print(f0, fullfile(svpath, 'Poly group average'), '-dsvg')
print(f0, fullfile(svpath, 'Poly group average'), '-dpng')

f1 = figure;
freqs = freqs(1:4000);
for i = 1:2

    % Plot
    subplot(1,2,i)
    plot(freqs, mufooof(:,i),'color', color(1,:), 'linewidth',2)
    hold on
    fill([freqs , fliplr(freqs)], sefooof(:,i), color(1,:), ...
        'edgecolor','none','facealpha', 0.1);
    title(sprintf('Foooof Detrended PSD %s', ch{i}));
    xticks([0 10 20 30 40])
    axis square
    box off
    grid on
    xlim([0 30])
    ylim([0 1])

    % Get alphaCF from 3rd order polynomial fit
    aRng = freqs(freqs>=7.5 & freqs<=15);
    [maxalpha,isMaxAlpha] = max(mufooof(freqs>=7.5 & freqs<=15,i));

    % Dot for IAF
    hold on
    scatter(aRng(isMaxAlpha), maxalpha, 25, color(2,:),'filled')

end

% Print
print(f1, fullfile(svpath, 'Fooof group average'), '-dsvg')
print(f1, fullfile(svpath, 'Fooof group average'), '-dpng')

f2 = figure;
for i = 1:2

    % Plot
    subplot(1,2,i)
    plot(fr, muwavelet(:,i),'color', color(1,:), 'linewidth',2)
    hold on
    fill([fr , fliplr(fr)], sewavelet(:,i), color(1,:), ...
        'edgecolor','none','facealpha', 0.1);
    title(sprintf('Wavelet Detrended PSD %s', ch{i}));
    xticks([0 10 20 30 40])
    axis square
    box off
    grid on
    xlim([0 30])

end

print(f2, fullfile(svpath, 'Wavelet group average'), '-dsvg')
print(f2, fullfile(svpath, 'Wavelet group average'), '-dpng')

% Get grand average
tmp = squeeze(mean(power));

f3 = figure;
for i = 1:2

    % Plot
    subplot(2, 1,i)
    imagesc(txx(txx <= 120), fxx, tmp(:,:,i)'); %pyy'
    title(sprintf('Grand Average MT Spectrums %s', ch{i}));
    xticks([0 20 40 60 80 100 120])
    colormap jet
    set(gca,'Ydir','normal')
    caxis([-10 10])
    ylim([5 30])
    pbaspect([3 1 1])
    yticks([5 10 15 20 25 30])
    colorbar
    xlim([20 120])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    box off

end
% % 
print(f3, fullfile(svpath, 'Grand Average MT Spectrums'), '-dsvg')
print(f3, fullfile(svpath, 'Grand Average MT Spectrums'), '-dpng')

% Get grand average
tmp = squeeze(mean(dpower));

f4 = figure;
for i = 1:2

    % Plot
    subplot(2, 1,i)
    imagesc(txx(txx <= 120), fxx, tmp(:,:,i)'); %pyy'
    title(sprintf('Grand Average detrended MT Spectrums %s', ch{i}));
    xticks([0 20 40 60 80 100 120])
    colormap jet
    set(gca,'Ydir','normal')
    caxis([-10 10])
    colorbar
    ylim([5 30])
    pbaspect([3 1 1])
    yticks([5 10 15 20 25 30])
    xlim([20 120])
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    box off

end

print(f4, fullfile(svpath, 'Grand Average detrended MT Spectrums'), '-dsvg')
print(f4, fullfile(svpath, 'Grand Average detrended MT Spectrums'), '-dpng')

fprintf('Mean poly: %.3f +/- %.3f\n', mean(acfpoly), std(acfpoly))
fprintf('Mean fooof: %.3f +/- %.3f\n', mean(acffooof), std(acffooof))

% Prep ero envelope comparisons
trough = squeeze(aband(1,:,:,:));
peak = squeeze(aband(2,:,:,:));

% Get means
mut = squeeze(mean(trough(:,:,2) - trough(:,:,1)));
mup = squeeze(mean(peak(:,:,2) - peak(:,:,1)));

% Get Standard errors
setrough = squeeze(std(trough(:,:,2) - trough(:,:,1)))./sqrt(20);
setrough = [mut + setrough , fliplr(mut - setrough)];
sepeak = squeeze(std(peak(:,:,2) - peak(:,:,1)))./sqrt(20);
sepeak = [mup + sepeak , fliplr(mup - sepeak)];

% Nonparametric permutation testing
circ = 0;
trial_level = 0;
iter = 1000;
correction = "cluster";
[pt,~,~] = npperm(trough(:,:,2), trough(:,:,1), iter, correction, circ, trial_level);
[pp,~,~] = npperm(peak(:,:,2), peak(:,:,1), iter, correction, circ, trial_level);


% Plot envelope comaprisons
f5 = figure;

% Plot signals
subplot(2,2,1)
plot(t, mut, 'color', color(1,:), 'linewidth',2)
hold on
fill(time, setrough, color(1,:), ...
    'edgecolor','none','facealpha', 0.1);
pbaspect([3 1 1])
title(sprintf('Envelope Difference Trough'))
ylim([-.2 .4] )
xlim([-250 500])
box off

% Plot permutation results
subplot(2,2,3)
fill(time, [pt,ones(1,length(pt))],[0 0 0], 'MarkerEdgeColor','none')
pbaspect([3 1 1])
xlim([-250 500])
set(gca,'ydir','reverse')
box off
ylim([0 0.05])

% Plot signals
subplot(2,2,2)
plot(t, mup, 'color', color(2,:), 'linewidth',2)
hold on
fill(time, sepeak, color(2,:), ...
    'edgecolor','none','facealpha', 0.1);
pbaspect([3 1 1])
title(sprintf('Envelope Difference Peak'))
ylim([-.2 .4] )
xlim([-250 500])
box off

% Plot permutation results
subplot(2,2,4)
fill(time, [pp,ones(1,length(pp))],[0 0 0], 'MarkerEdgeColor','none')
pbaspect([3 1 1])
xlim([-250 500])
set(gca,'ydir','reverse')
box off
ylim([0 0.05])

print(f5, fullfile(svpath, 'ERO envelope Comaprison'), '-dsvg')
print(f5, fullfile(svpath, 'ERO envelope Comaprison'), '-dpng')


% Compare IAF Estimates
f6 = figure;
for i = 1:2
    
    % Get model
    lm = fitlm(acfpoly(:,i), acfwelch(:,i));
    
    % Plot
    subplot(1,2,i)
    scatter(acfpoly(:,i), acfwelch(:,i), 70, color(1,:), 'filled')
    xlabel('Polynomial IAF')
    ylabel('Welch IAF')
    grid on
    axis square
    box off
    hold on
    xlim([7 15])
    ylim([7 15])

    % Plot Identity line
    plot([7 15], [7 15], 'color', color(1,:), 'linewidth',2)
    title(sprintf('Fooof vs. Polynomial Detrend %s', ch{i}));

    % Add Model text
    text(11,9,sprintf('R^2 = %.3f \np-value = %.4f',...
        lm.Rsquared.Adjusted, lm.ModelFitVsNullModel.Pvalue));
    
end

% Save
svpath = fullfile('C:\Toolbox\MATLAB\Figures\All ERPs and EROs\');
print(f6, fullfile(svpath, 'Welch vs poly'), '-dsvg')
print(f6, fullfile(svpath, 'Welch vs poly'), '-dpng')

f7 = figure;
for i = 1:2

    % Plot
    subplot(1,2,i)
    plot(wfxx, muwelch(:,i),'color', color(1,:), 'linewidth',.5)
    hold on
    fill([wfxx ; flipud(wfxx)], sewelch(:,i), color(1,:), ...
        'edgecolor','none','facealpha', 0.1);
    title(sprintf('Polynomial Detrended Welch PSD %s', ch{i}))
    xticks([0 10 20 30 40])
    axis square
    box off
    grid on
    xlim([0 30])
    ylim([-4 6])


    % Get alphaCF from 3rd order polynomial fit
    aRng = fxx(fxx>=7.5 & fxx<=15);
    [maxalpha,isMaxAlpha] = max(mupoly(fxx>=7.5 & fxx<=14,i));

    % Dot for IAF
    hold on
    scatter(aRng(isMaxAlpha), maxalpha, 25, color(2,:),'filled')
aRng(isMaxAlpha)
end
% 
print(f7, fullfile(svpath, 'Welch group average'), '-dsvg')
print(f7, fullfile(svpath, 'Welch group average'), '-dpng')


% Compare IAF Estimates
f8 = figure;
for i = 1:1
    
    % Get model
    lm = fitlm(tbl, acfpoly(:,1));
    
    % Plot
    scatter(tbl, acfpoly(:,1), 70, color(1,:), 'filled')
    xlabel('Original IAF')
    ylabel('Polynomial IAF')
    grid on
    axis square
    box off
    hold on
    xlim([7 15])
    ylim([7 15])

    % Plot Identity line
    plot([7 15], [7 15], 'color', color(1,:), 'linewidth',2)
    title(sprintf('Original vs. Polynomial Detrend %s', ch{i}));

    % Add Model text
    text(11,9,sprintf('R^2 = %.3f \np-value = %.4f',...
        lm.Rsquared.Adjusted, lm.ModelFitVsNullModel.Pvalue));
    
end

% Save
svpath = fullfile('C:\Toolbox\MATLAB\Figures\All ERPs and EROs\');
print(f8, fullfile(svpath, 'Poly vs Scott orignal'), '-dsvg')
print(f8, fullfile(svpath, 'Poly vs Scott orignal'), '-dpng')

