% Created 2021-05-25 Tylor J Harlow - tylor.harlow@uconn.edu

%% Housekeeping
tic;
clear;
close all;
clc;

svpath = '/Applications/Toolbox/MATLAB/Figures/Plots wout 126';
savePlots = 1;
updateTable = 0;

% Load current data table
datafilename = '/Applications/Toolbox/MATLAB/TwoPhaseExp_table.mat';
load(datafilename);

%% Analysis parameters

% Define pre-processing for raw EEG data
params = echtparams('filterImplementation',0,'zeroPhase',0,'tWin',[-0.250 0.900]);

CH = {'Fpz','Fz'}; % Electrode channel labels

% Anonymous function to calculate instantaneous frequency from phase
iPhaseToiFreq = @(fs,phi)(fs/(2*pi)*diff(unwrap(phi)));

syncToECHT = 1; % Use instantaneous phase values from ECHT device? [0:no, 1:yes]

%% Load data

% Get paths
pname = genpath('/Applications/Toolbox/SubjectData');
pname = strsplit(pname,':')';

% Locate Two-Phase ERP datasets
pname = pname(contains(pname,'TwoPhaseERP85dB') & ~contains(pname,'(Low'));

% Exclude Subjects
exSub = {'Sub_151','Sub_152', '126'};
if(~isempty(exSub))  
    pname = pname(~contains(pname,exSub)); % exclude subject(s) from analysis
end


%% Run analyses

% Whats my sample?
n = length(pname);

% Grab data
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

    % Store sampling rate
    fs = ERP.fs(1);

    % Separate by Phase
    for ph = 1:2
        pIdx{ph} = intersect(ERP.goodTrials,find(round(ERP.echtPhaseLabel)==ERP.setPhases(ph)));
    end

    % Loop through elctrodes
    for j = 1:2

        % Grab wavelet transform
        wt = squeeze(permute(ERP.cwt(:,j,:,:),[4 2 3 1]));

        % Separate by condition
        opt = (squeeze(wt(pIdx{1},:,:)));
        pes = (squeeze(wt(pIdx{2},:,:)));

        % Intertrial phase clustering
        optimalitc(k,:,:,j) = abs(mean(exp(1i*angle(opt))));
        pessimalitc(k,:,:,j) = abs(mean(exp(1i*angle(pes))));

        % Evoked Power calculation
        optimalevoked(k,:,:,j) = (mean(opt).*mean(conj(opt)));
        pessimalevoked(k,:,:,j) = (mean(pes).*mean(conj(pes)));

        % Total Power calculation
        optimal(k,:,:,j) = (mean(opt.*conj(opt)));
        pessimal(k,:,:,j) = (mean(pes.*conj(pes)));

        % ERD/ERS
        fs = ERP.fs(1);
        erp = mean(squeeze(permute(ERP.erp(:,j,:), [3 1 2])));
        noerp = (squeeze(permute(ERP.erp(:,j,:), [3 1 2])) - erp)';
        [ind f] = wavecwt(noerp, fs, [1 70], [4 12], 140);
        optighost(k,:,:,j) = squeeze(mean(abs(ind(:,pIdx{1},:)).^2,2))';
        pesighost(k,:,:,j) = squeeze(mean(abs(ind(:,pIdx{2},:)).^2,2))';
   

    end
end

%% Finalize Processing

% Baseline Period
idx = find(t< -50 & t>-250);

%================= CWT Spectrogram ======================%
optbc = squeeze(mean(10*log10(optimal./mean(optimal(:,:,idx,:),3)))); 
pesbc = squeeze(mean(10*log10(pessimal./mean(pessimal(:,:,idx,:),3))));

%================= Evoked ACtivity ======================%
% Prestimulus baseline correction
optevbc = squeeze(mean(10*log10(optimalevoked./mean(optimalevoked(:,:,idx,:),3)))); 
pesevbc = squeeze(mean(10*log10(pessimalevoked./mean(pessimalevoked(:,:,idx,:),3))));
% Whole-time average baseline correction
optevbcmu = squeeze(mean(10*log10(optimalevoked./mean(optimalevoked,3)))); 
pesevbcmu = squeeze(mean(10*log10(pessimalevoked./mean(pessimalevoked,3))));

%================= Induced ACtivity =====================%
optghost = squeeze(mean(10*log10(optighost./mean(optighost(:,:,idx,:),3)))); 
pesghost = squeeze(mean(10*log10(pesighost./mean(pesighost(:,:,idx,:),3))));

%================= ITPC =================================%
optimalitcav = squeeze(mean(optimalitc));
pessimalitcav = squeeze(mean(pessimalitc));


%% Figures

co2 = [0.0000 0.4470 0.7410
       0.8500 0.10 0.0980];

%================= ITPC Figures ===========================%
f1 = figure
contourf(t,f,(optimalitcav(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([0 0.5])
c = colorbar
c.Label.String = 'ITPC'
colormap jet
axis square
box off
ylabel("Frequency")
title("Trough ITCPC FPZ")

f2 = figure
contourf(t,f,(pessimalitcav(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([0 0.5])
c = colorbar
c.Label.String = 'ITPC'
colormap jet
axis square
box off
title("Peak ITPC FPZ")

f3 = figure
contourf(t,f,(optimalitcav(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([0 0.5])
c = colorbar
c.Label.String = 'ITPC'
colormap jet
axis square
box off
xlabel('Time (ms)')
ylabel("Frequency")
title("Trough ITCPC FZ")

f4 = figure
contourf(t,f,(pessimalitcav(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([0 0.5])
c = colorbar
c.Label.String = 'ITPC'
colormap jet
axis square
box off
xlabel('Time (ms)')
title("Peak ITPC FZ")

%%
%================= Spectrogram ========================%
f5 = figure
contourf(t,f,real(optbc(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-1 1])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
ylabel("Frequency")
title("Trough Total Power  FPZ")

f6 = figure
contourf(t,f,real(pesbc(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-1 1])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
title("Peak Total Power  FPZ")

f7 = figure
contourf(t,f,real(optbc(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-1 1])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
xlabel('Time (ms)')
ylabel("Frequency")
title("Trough Total Power  FZ")

f8 = figure
contourf(t,f,real(pesbc(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-1 1])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
xlabel('Time (ms)')
title("Peak Total Power  FZ")


%%
%================= Evoked Spectrogram ================%
f9 = figure
contourf(t,f,real(optevbc(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-20 20])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
ylabel("Frequency")
title("Trough Evoked Power  FPZ")

f10 = figure
contourf(t,f,real(pesevbc(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-20 20])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
title("Peak Evoked Power  FPZ")

f11 = figure
contourf(t,f,real(optevbc(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-20 20])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
xlabel('Time (ms)')
ylabel("Frequency")
title("Trough Evoked Power  FZ")

f12 = figure
contourf(t,f,real(pesevbc(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-20 20])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
xlabel('Time (ms)')
title("Peak Evoked Power  FZ")


%%
%================= Induced Spectrogram ================%
f13 = figure
contourf(t,f,real(optghost(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-1 1])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
ylabel("Frequency")
title("Trough Induced Power  FPZ")

f14 = figure
contourf(t,f,real(pesghost(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-1 1])
c = colorbar
c.Label.String = 'dB'
colormap jet
% caxis([0 1])
axis square
box off
title("Peak Induced Power  FPZ")

f15 = figure
contourf(t,f,real(optghost(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-1 1])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
xlabel('Time (ms)')
ylabel("Frequency")
title("Trough Induced Power  FZ")

f16 = figure
contourf(t,f,real(pesghost(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-1 1])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
xlabel('Time (ms)')
title("Peak Induced Power  FZ")


%% Evoked spectrogram - whole time mean baseline correction
%================= Evoked Spectrogram ================%
f17 = figure
contourf(t,f,real(optevbcmu(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-10 10])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
ylabel("Frequency")
title("Trough Evoked Power  FPZ (baseline = mean)")

f18 = figure
contourf(t,f,real(pesevbcmu(:,:,1)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-10 10])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
title("Peak Evoked Power  FPZ (baseline = mean)")

f19 = figure
contourf(t,f,real(optevbcmu(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-10 10])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
xlabel('Time (ms)')
ylabel("Frequency")
title("Trough Evoked Power  FZ (baseline = mean)")

f20 = figure
contourf(t,f,real(pesevbcmu(:,:,2)),30,'linecolor','none')
set(gca,'YDir','normal')
xlim([-200 900]); ylim([2 70])
hold on
plot(t, 8*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
plot(t, 12*ones(size(t)), '--','color', [0.5 0.5 0.5], 'linewidth',2)
caxis([-10 10])
c = colorbar
c.Label.String = 'dB'
colormap jet
axis square
box off
xlabel('Time (ms)')
title("Peak Evoked Power  FZ (baseline = mean)")


%% Save Figures
savePlots = 1;
if savePlots

    % Setup 

    fprintf('Saving plots...\n');
    svpath =  '/Applications/Toolbox/MATLAB/Figures/Spectrograms wout 126/';

    %================= Save Figures ==================%

    % SVGS
    print(f1,fullfile(svpath,'Trough ITCPC FPZ.svg'),'-dsvg','-painters');
    print(f2,fullfile(svpath,'Peak ITCPC FPZ.svg'),'-dsvg','-painters');
    print(f3,fullfile(svpath,'Trough ITCPC FZ.svg'),'-dsvg','-painters');
    print(f4,fullfile(svpath,'Peak ITCPC FZ.svg'),'-dsvg','-painters');
    
    print(f5,fullfile(svpath,'Trough Total Power  FPZ.svg'),'-dsvg','-painters');
    print(f6,fullfile(svpath,'Peak Total Power  FPZ.svg'),'-dsvg','-painters');
    print(f7,fullfile(svpath,'Trough Total Power  FZ.svg'),'-dsvg','-painters');
    print(f8,fullfile(svpath,'Peak Total Power  FZ.svg'),'-dsvg','-painters');
    
    print(f9,fullfile(svpath,'Trough Evoked Power  FPZ.svg'),'-dsvg','-painters');
    print(f10,fullfile(svpath,'Peak Evoked Power  FPZ.svg'),'-dsvg','-painters');
    print(f11,fullfile(svpath,'Trough Evoked Power  FZ.svg'),'-dsvg','-painters');
    print(f12,fullfile(svpath,'Peak Evoked Power  FZ.svg'),'-dsvg','-painters');
    
    print(f13,fullfile(svpath,'Trough Induced Power  FPZ.svg'),'-dsvg','-painters');
    print(f14,fullfile(svpath,'Peak Induced Power  FPZ.svg'),'-dsvg','-painters');
    print(f15,fullfile(svpath,'Trough Induced Power  FZ.svg'),'-dsvg','-painters');
    print(f16,fullfile(svpath,'Peak Induced Power  FZ.svg'),'-dsvg','-painters');

    print(f17,fullfile(svpath,'Trough Evoked Power  FPZ (no baseline).svg'),'-dsvg','-painters');
    print(f18,fullfile(svpath,'Peak Evoked Power  FPZ (no baseline).svg'),'-dsvg','-painters');
    print(f19,fullfile(svpath,'Trough Evoked Power  FZ (no baseline).svg'),'-dsvg','-painters');
    print(f20,fullfile(svpath,'Peak Evoked Power  FZ (no baseline).svg'),'-dsvg','-painters');

    % PNGs
    print(f1,fullfile(svpath,'Trough ITCPC FPZ.png'),'-dpng');
    print(f2,fullfile(svpath,'Peak ITCPC FPZ.png'),'-dpng');
    print(f3,fullfile(svpath,'Trough ITCPC FZ.png'),'-dpng');
    print(f4,fullfile(svpath,'Peak ITCPC FZ.png'),'-dpng');
    
    print(f5,fullfile(svpath,'Trough Total Power  FPZ.png'),'-dpng');
    print(f6,fullfile(svpath,'Peak Total Power  FPZ.png'),'-dpng');
    print(f7,fullfile(svpath,'Trough Total Power  FZ.png'),'-dpng');
    print(f8,fullfile(svpath,'Peak Total Power  FZ.png'),'-dpng');
    
    print(f9,fullfile(svpath,'Trough Evoked Power  FPZ.png'),'-dpng');
    print(f10,fullfile(svpath,'Peak Evoked Power  FPZ.png'),'-dpng');
    print(f11,fullfile(svpath,'Trough Evoked Power  FZ.png'),'-dpng');
    print(f12,fullfile(svpath,'Peak Evoked Power  FZ.png'),'-dpng');
    
    print(f13,fullfile(svpath,'Trough Induced Power  FPZ.png'),'-dpng');
    print(f14,fullfile(svpath,'Peak Induced Power  FPZ.png'),'-dpng');
    print(f15,fullfile(svpath,'Trough Induced Power  FZ.png'),'-dpng');
    print(f16,fullfile(svpath,'Peak Induced Power  FZ.png'),'-dpng');

    print(f17,fullfile(svpath,'Trough Evoked Power  FPZ (no baseline).svg'),'-dpng','-painters');
    print(f18,fullfile(svpath,'Peak Evoked Power  FPZ (no baseline).svg'),'-dpng','-painters');
    print(f19,fullfile(svpath,'Trough Evoked Power  FZ (no baseline).svg'),'-dpng','-painters');
    print(f20,fullfile(svpath,'Peak Evoked Power  FZ (no baseline).svg'),'-dpng','-painters');

end
