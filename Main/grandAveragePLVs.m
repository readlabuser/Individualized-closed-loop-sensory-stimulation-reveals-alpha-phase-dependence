% grandAveragePLVs.m
%
% Re-creation of "Mean Actual Phase histograms.png" from SBIR
%
% Created 2022-09-22 scott.bressler@elemindtech.com
% Adapted 2023-08-16 tylor.harlow@uconn.edu

%% Housekeeping
clear;
close all;
clc;

% Grab the default color order
co = [0.0000 0.4470 0.7410
       0.8500 0.10 0.0980];

%% Load Data

% Select subjects to exclude from group analysis
exSub = {'Sub_151','Sub_152','126'}; %,'Sub_102','Sub_104','Sub_126','Sub_138'};

%% Analysis parameters
% LabView streaming command
params.strcmd = 'str,m,dt,eeg1,eeg2,out,iphs,ps_pulse_count,ps_pulse_id,ps_onphs,ps_offphs,ps_onontime,ps_durbyphs,ps_durbytime,ps_isi';

% Define pre-processing for raw EEG data
params = echtparams('zeroPhase',0);

CH = {'Fpz','Fz'}; % Electrode channel labels

% Anonymous function to calculate instantaneous frequency from phase
iPhaseToiFreq = @(fs,phi)(fs/(2*pi)*diff(unwrap(phi)));

syncToECHT = 1; % Use instantaneous phase values from ECHT device? [0:no, 1:yes]

%% Load summarized data table
% Load current data table
datafilename = '/Applications/Toolbox/MATLAB/TwoPhaseExp_table.mat';      %% Changed for Windows
load(datafilename);
dataTable = sortrows(dataTable);

%% Load data
pname = genpath('/Applications/Toolbox/SubjectData');  %% Changed for Windows
pname = strsplit(pname,':')'; %% Switched

% Locate Random Phase ERP datasets
pname1 = pname(contains(pname,'ISI750ms85dB') & ~contains(pname,'(Low'));

% Locate Two-Phase ERP datasets
pname2 = pname(contains(pname,'TwoPhaseERP85dB') & ~contains(pname,'(Low'));

% Exclude subjects?
if(~isempty(exSub))  
    pname1 = pname1(~contains(pname1,exSub)); % exclude subject(s) from analysis
    pname2 = pname2(~contains(pname2,exSub)); % exclude subject(s) from analysis
end

% Start Grand error array
error = [];
echterror = [];

for k = 1:length(pname2)

    % Get subject info
    tmp = strsplit(pname2{k},filesep);
    Sub_ID = tmp{contains(tmp,'Sub_')};
    subIdx = contains(dataTable.SubID,Sub_ID([5:end]));
    sub{k} = Sub_ID;

    % Setup alpha center frequency
    params.alphaCF = dataTable(subIdx,:).IAF(1);

    % Run
    [ERP,params] = batchERP(pname2{k},params);

    % Check order of optimal/pessimal phase list in ERP structure against the
    % calculated phases in 'dataTable.' Fix wrap-around error if necessary
    if(dataTable(subIdx,:).OptPhase > dataTable(subIdx,:).PesPhase)
        ERP.setPhases = ERP.setPhases([2 1]);
        fprintf('Switcharoo!\n')
    end

    % Store IAF 
    iaf(k) = params.alphaCF;

    % Separate by Phase & Grab Alpha Stuff
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

        % Remove ERP
        fs = ERP.fs(1);
        for j = 1:2
            erp = mean(squeeze(permute(ERP.erp(:,j,:), [3 1 2])));
            noerp = (squeeze(permute(ERP.erp(:,j,:), [3 1 2])) - erp)';
            [indalpha(:,j,:,:), f] = wavecwt(noerp, fs, [1 70], [4 12], 140);
        end
        indalpha = permute(indalpha, [3 1 2 4]);

        % Get Induced
        [~,IAF] = min(abs(dataTable(subIdx,:).IAF(1) - f));
        indalpha = abs(squeeze(indalpha(pIdx{ph},:,:,IAF))).^2;
        ind(k,:,:,ph) = squeeze(mean(indalpha));

        % Cleanup
        clearvars indalpha
    end

    % Onset phase-locking by condition
    OptOnPLV(k) = abs(mean(exp(1i*(deg2rad(wrapTo180(rad2deg(ERP.aPhaseOn(pIdx{1}))) - ERP.setPhases(1))))));
    PesOnPLV(k) = abs(mean(exp(1i*(deg2rad(wrapTo180(rad2deg(ERP.aPhaseOn(pIdx{2}))) - ERP.setPhases(2))))));
    OptECHTPLV(k) = abs(mean(exp(1i*(deg2rad(wrapTo180(rad2deg(ERP.echtPhase(pIdx{1}))) - ERP.setPhases(1))))));
    PesECHTPLV(k) = abs(mean(exp(1i*(deg2rad(wrapTo180(rad2deg(ERP.echtPhase(pIdx{2}))) - ERP.setPhases(2))))));
    OptOn(k) = circ_mean(ERP.aPhaseOn(pIdx{1}));
    PesOn(k) = circ_mean(ERP.aPhaseOn(pIdx{2}));

    % Error by subject - both conditions for hilbert alpha
    accuracy(k) = abs(mean(exp(1i*(deg2rad(wrapTo360(rad2deg(ERP.aPhaseOn(ERP.goodTrials))) - ERP.echtPhaseLabel(ERP.goodTrials))))));

    % Error by subject - both conditions for echt
    echt(k) = abs(mean(exp(1i*(deg2rad(wrapTo360(rad2deg(ERP.echtPhase(ERP.goodTrials))) - ERP.echtPhaseLabel(ERP.goodTrials))))));

    % Phase lag for hilbert alpha
    phaselag(k) = circ_mean((deg2rad(rad2deg(ERP.aPhaseOn(ERP.goodTrials)) - ERP.echtPhaseLabel(ERP.goodTrials))));

    % Phase lag for ECHT alpha
    echtlag(k) = circ_mean((deg2rad(rad2deg(ERP.echtPhase(ERP.goodTrials)) - ERP.echtPhaseLabel(ERP.goodTrials))));

    % Get total error for hilbert alpha
    error = cat(1, error, (deg2rad(rad2deg(ERP.aPhaseOn(ERP.goodTrials)) - ERP.echtPhaseLabel(ERP.goodTrials))));

    % Get total error for echt alpha
    echterror = cat(1, error, (deg2rad(rad2deg(ERP.echtPhase(ERP.goodTrials)) - ERP.echtPhaseLabel(ERP.goodTrials))));

    % Error by subject - both conditions for wavelet alpha
    waccuracy(k) = abs(mean(exp(1i*deg2rad(wrapTo180(rad2deg(ERP.waPhaseOn(ERP.goodTrials))) - ERP.echtPhaseLabel(ERP.goodTrials)))));

    % Phase lag for wavelet alpha
    wphaselag(k) = circ_mean((deg2rad(rad2deg(ERP.waPhaseOn(ERP.goodTrials)) - ERP.echtPhaseLabel(ERP.goodTrials))));

    % Get total error for wavelet alpha
    werror = cat(1, error, (deg2rad(wrapTo180(rad2deg(ERP.waPhaseOn(ERP.goodTrials))) - ERP.echtPhaseLabel(ERP.goodTrials))));

    % Intertrial phase clusterings
    hitpc(k,:) = [abs(mean(1i*ERP.aPhaseOn(pIdx{1}))), abs(mean(1i*ERP.aPhaseOn(pIdx{2})))];
    witpc(k,:) = [abs(mean(1i*ERP.waPhaseOn(pIdx{1}))), abs(mean(1i*ERP.waPhaseOn(pIdx{2})))];

    % Process Twophase ISI
    twoISI(k) = mean(ERP.isi);

    % Grab Random Condition & ISI
    ERP = batchERP(pname1{k}, params);
    randISI(k) = mean(ERP.isi);

    % Count random trials
    trials(k) = size(ERP.erp,3);

    % Time length
    l(k) = max(ERP.time{1});

    % Subject ID
    ID(k) = dataTable(subIdx,:).SubID;

end


%% Circ-stat Details
nbins = deg2rad(0:20:360); % polar histogram bins

% % Summary stats -- Optimal (alpha trough)
% mAng_OptOn = rad2deg(circ_mean(OptOn'));
% sdAng_OptOn = rad2deg(circ_std(OptOn'));
% PLV_OptOn = mean(OptOnPLV');
% 
% % Summary stats -- Pessimal (alpha peak)
% mAng_PesOn = rad2deg(circ_mean(PesOn'));
% sdAng_PesOn = rad2deg(circ_std(PesOn'));
% PLV_PesOn = mean(PesOnPLV');
% PLV_total = cat(1,OptOnPLV',PesOnPLV');
% stdtotal = circ_std(phaselag');
% PLV_total = mean(PLV_total);
% mtotal = circ_mean(phaselag');
% PLV_ECHT = mean(cat(1,OptECHTPLV',PesECHTPLV'));
% mecht = circ_mean(echtlag');
% stdecht = circ_std(echtlag');

% Summary stats -- Optimal (alpha trough)
mAng_OptOn = rad2deg(circ_mean(OptOn'));
[~,sdAng_OptOn] = circ_std(OptOn');
sdAng_OptOn = rad2deg(sdAng_OptOn);
PLV_OptOn = mean(OptOnPLV');

% Summary stats -- Pessimal (alpha peak)
mAng_PesOn = rad2deg(circ_mean(PesOn'));
[~,sdAng_PesOn] = circ_std(PesOn');
sdAng_PesOn = rad2deg(sdAng_PesOn);
PLV_PesOn = mean(PesOnPLV');
PLV_total = cat(1,OptOnPLV',PesOnPLV');
[~,stdtotal] = circ_std(phaselag');
% stdtotal = rad2deg(stdtotal);
% PLV_total = mean(PLV_total);
mtotal = circ_mean(phaselag');
PLV_ECHT = mean(cat(1,OptECHTPLV,PesECHTPLV));
mecht = circ_mean(echtlag');
[~,stdecht] = circ_std(echtlag');
% stdecht = rad2deg(stdecht);

%% Power Details

% Get median split phase error
low = find(phaselag <= median(phaselag));
high = find(phaselag > median(phaselag));

% Baseline normalize power measures
t = ERP.t0;
n = length(phaselag);
time = [t ; flipud(t)];
idxBL = find(t< -50 & t>-250);
% evo = evo - mean(evo(:,idxBL,:,:),2);
tot = tot - mean(tot(:,idxBL,:,:),2);
ind = ind - mean(ind(:,idxBL,:,:),2);

% Mean and SEM for all measures - note: SEM are adapted for future use of 
% Loop through medians
medis = {low, high};close
for j = 1:2
    % fill functions
    evomu(:,:,:,j) = squeeze(mean(evo(medis{j},:,:,:))); % Evoked
    se = squeeze(std(evo(medis{j},:,:,:)))./sqrt(n);
    evose(:,:,:,j) = cat(1, evomu(:,:,:,j) + se,...
        flip(evomu(:,:,:,j) - se,1));

    totmu(:,:,:,j) = squeeze(mean(tot(medis{j},:,:,:))); % Total
    se = squeeze(std(tot(medis{j},:,:,:)))./sqrt(n);
    totse(:,:,:,j) = cat(1, totmu(:,:,:,j) + se,...
        flip(totmu(:,:,:,j) - se,1));

    indmu(:,:,:,j) = squeeze(mean(ind(medis{j},:,:,:))); % Induced
    se = squeeze(std(ind(medis{j},:,:,:)))./sqrt(n);
    indse(:,:,:,j) = cat(1, indmu(:,:,:,j) + se,...
        flip(indmu(:,:,:,j) - se,1));

    itcmu(:,:,:,j) = squeeze(mean(itc(medis{j},:,:,:))); % ITPC
    se = squeeze(std(itc(medis{j},:,:,:)))./sqrt(n);
    itcse(:,:,:,j) = cat(1, itcmu(:,:,:,j) + se,...
        flip(itcmu(:,:,:,j) - se,1));

    eromu(:,:,:,j) = squeeze(mean(ero(medis{j},:,:,:))); % ERO
    se = squeeze(std(ero(medis{j},:,:,:)))./sqrt(n);
    erose(:,:,:,j) = cat(1, eromu(:,:,:,j) + se,...
        flip(eromu(:,:,:,j) - se,1));

    envmu(:,:,:,j) = squeeze(mean(env(medis{j},:,:,:))); % Envelope
    se = squeeze(std(env(medis{j},:,:,:)))./sqrt(n);
    envse(:,:,:,j) = cat(1, envmu(:,:,:,j) + se,...
        flip(envmu(:,:,:,j) - se,1));
end

%% Figure 1: Optimal (alpha trough)
if(exist('f1','var') && ishandle(f1)), close(f1); end
f1 = figure('Position',[147 377 560 420]);

ax1 = polaraxes;
hold on;
polarhistogram(ax1,(OptOn),nbins,'Normalization','probability',...
    'FaceAlpha',0.5,'EdgeColor',co(1,:),'FaceColor',co(1,:));
% polarhistogram(ax1,(OptOnPLV),nbins,'Normalization','probability',...
%     'FaceAlpha',0.5,'EdgeColor',[0.5 0.5 0.5]);

% Mean angle and PLV
polarplot(ax1,[0 deg2rad(mAng_OptOn)],[0 PLV_OptOn],...
    'LineWidth',4,'Color',co(1,:));

% SD angle arc
polarplot(ax1,deg2rad(mAng_OptOn)+deg2rad(linspace(-sdAng_OptOn,sdAng_OptOn,100)),...
    PLV_OptOn*ones(size(linspace(-sdAng_OptOn,sdAng_OptOn,100))),...
    'LineWidth',2,'Color',co(1,:));

% Target Onset/Offsets
polarplot(ax1,[0 circ_mean(deg2rad(134))],[0 PLV_OptOn],...
    'LineWidth',2,'Color','k','LineStyle','--');

% Adjust axes properties
ax1.Title.String = sprintf('Alpha Trough Stimulation (n = %d)',numel(OptOn));
ax1.Title.FontSize = 16;
ax1.RAxisLocation = 190;
ax1.FontSize = 14;
rlim(ax1, [0 0.6]);



%% Figure 2: Pessimal (alpha peak)
if(exist('f2','var') && ishandle(f2)), close(f2); end
f2 = figure('Position',[708 377 560 420]);

ax2 = polaraxes;
hold on;
polarhistogram(ax2,PesOn,nbins,'Normalization','probability',...
    'FaceAlpha',0.5,'EdgeColor',co(2,:),'FaceColor',co(2,:));
% polarhistogram(ax1,(PesOnPLV),nbins,'Normalization','probability',...
%     'FaceAlpha',0.5,'EdgeColor',[0.5 0.5 0.5]);

% Mean angle and PLV
polarplot(ax2,[0 deg2rad(mAng_PesOn)],[0 PLV_PesOn],...
    'LineWidth',4,'Color',co(2,:));

% SD angle arc
polarplot(ax2,deg2rad(mAng_PesOn)+deg2rad(linspace(-sdAng_PesOn,sdAng_PesOn,100)),...
    PLV_PesOn*ones(size(linspace(-sdAng_PesOn,sdAng_PesOn,100))),...
    'LineWidth',2,'Color',co(2,:));

% Target Onset/Offsets
polarplot(ax2,[0 circ_mean(deg2rad(314))],[0 PLV_PesOn],...
    'LineWidth',2,'Color','k','LineStyle','--');

% Adjust axes properties
ax2.Title.String = sprintf('Alpha Peak Stimulation (n = %d)',numel(PesOn));
ax2.Title.FontSize = 16;
ax2.RAxisLocation = 15;
ax2.FontSize = 14;
rlim(ax2, [0 0.6]);


%% Figure 3: PLVs - Acausal Validity
if(exist('f3','var') && ishandle(f3)), close(f3); end
f3 = figure('Position',[708 377 560 420]);

ax3 = polaraxes;
hold on;
polarhistogram(ax3,phaselag,nbins,'Normalization','probability',...
    'FaceAlpha',0.5,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);

% Mean angle and PLV
polarplot(ax3,[0 mtotal],[0 PLV_total],...
    'LineWidth',4,'Color',[0.5 0.5 0.5]);

% SD angle arc
polarplot(ax3,mtotal+(linspace(-stdtotal,stdtotal,100)),...
    PLV_total*ones(size(linspace(-stdtotal,stdtotal,100))),...
    'LineWidth',2,'Color',[0.5 0.5 0.5]);

% Target Onset/Offsets
polarplot(ax3,[0 circ_mean(deg2rad(0))],[0 PLV_total],...
    'LineWidth',1,'Color','k','LineStyle','--');

text(-2,0.2,sprintf('PLV: %0.3f',PLV_total),'FontWeight','bold', 'FontSize', 14)
text(-2.2,0.3,sprintf('PE: %0.3f (\x00B1 %0.3f)',circ_mean(phaselag'), circ_std(phaselag')),'FontWeight','bold', 'FontSize', 14)

% Adjust axes properties
ax3.Title.String = sprintf('Acausal, Group Phase Error (n = %d)',numel(OptOn));
ax3.Title.FontSize = 16;
ax3.RAxisLocation = 15;
ax3.FontSize = 14;
rlim(ax3, [0 0.6]);

%% Figure 3b: PLVs - ECHT Reliability
f3b = figure('Position',[708 377 560 420]);

ax3b = polaraxes;
hold on;
polarhistogram(ax3b,echtlag,nbins,'Normalization','probability',...
    'FaceAlpha',0.5,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);

% Mean angle and PLV
polarplot(ax3b,[0 mecht],[0 PLV_ECHT],...
    'LineWidth',4,'Color',[0.5 0.5 0.5]);

% SD angle arc
polarplot(ax3b,mecht+(linspace(-stdecht,stdecht,100)),...
    PLV_ECHT*ones(size(linspace(-stdecht,stdecht,100))),...
    'LineWidth',2,'Color',[0.5 0.5 0.5]);

% Target Onset/Offsets
polarplot(ax3b,[0 circ_mean(deg2rad(0))],[0 PLV_ECHT],...
    'LineWidth',1,'Color','k','LineStyle','--');

text(-2,0.35,sprintf('PLV: %0.3f',PLV_ECHT),'FontWeight','bold', 'FontSize', 14)
text(-2.2,0.5,sprintf('PE: %0.3f (\x00B1 %0.3f)',circ_mean(echtlag'), circ_std(echtlag')),'FontWeight','bold', 'FontSize', 14)

% Adjust axes properties
ax3b.Title.String = sprintf('ECHT, Group Phase Error (n = %d)',numel(OptOn));
ax3b.Title.FontSize = 16;
ax3b.RAxisLocation = 15;
ax3b.FontSize = 14;
rlim(ax3b, [0 1]);

%% Figure 5: Peak Phase versus P50 latency
% Data for conceptualization of alpha center frequency/P50 latency relationship to Onset Phase
p50 = 1e-3*(min(dataTable.P1lat(:,1)):80)'; % P50 peak latency (seconds)
alphaCF = 8:12; % Alpha center frequency (Hz)

Y = wrapTo360(-360*p50*alphaCF); % matrix of latency:phase relationship
minY = min(Y);
maxY = max(Y);

if(exist('f5') & ishandle(f5)), close(f5); end

f4 = figure('Position',[350 1457 1000 500]);
colororder(copper(length(alphaCF)));

ax5a = axes('Position',[0.1100 0.1100 0.4000 0.8150]);
hold on;
for k = 1:size(Y,2)
    plot(ax5a,p50*1e3,Y(:,k),'.-','LineWidth',1,'color',[0 0 k*.2]);
end
plot(ax5a,dataTable.P1lat(:,1),dataTable.OptPhase,'ko','LineWidth',1,...
    'MarkerFaceColor','w','MarkerSize',8);

% Subject 136
plot(ax5a,dataTable(13,:).P1lat(:,1),dataTable(13,:).OptPhase,'o','Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',8);


ax5a.XLabel.String = 'P1 Peak Latency (ms)';
ax5a.YLabel.String = 'Alpha Trough Phase at Onset (degrees)';
% ax5a.XMinorGrid = 'on';
% ax5a.YGrid = 'on';
ax5a.Title.String = 'Alpha Trough Onset Phases Predicted by P1 Latency';

% Plot reference point
plot(ax5a, 50, 180, 'o', 'color', [1 0.5 0.8],'LineWidth',1,...
    'MarkerFaceColor',[1 0.5 0.8],'MarkerSize',8);

% Add prediction equation
text(ax5a,0.02,0.05,...
    '$\mathrm{onset}\angle^{\circ}=-360^{\circ}\times\mathrm{P1_{latency} (sec)}\times\alpha_{CF}(Hz)$',...
    'Interpreter','latex','FontSize',12,'Units','normalized');

% Add phase text labels
for k = 1:length(alphaCF)
   text(ax5a,1e3*p50(end)+0.5,Y(end,k),sprintf('%1.1f Hz',alphaCF(k))); 
end
ax5a.XLim = [min(p50*1e3) 75];

ax5b = polaraxes('Position',[0.5500 0.5500 0.6975/2 0.7335/2]);
hold on;
% for k = 1:length(alphaCF)
% %     polarplot(ax5b,[0 deg2rad(minY(k))],[0 1],'.--','LineWidth',1);
% %     polarplot(ax5b,[0 deg2rad(maxY(k))],[0 1],'o-','LineWidth',1);
%     polarplot(ax5b,[0 deg2rad(minY(k):maxY(k)) 0],alphaCF(k)*[0 ones(size(minY(k):maxY(k))),0],...
%         '-','LineWidth',2,'color',[0 0 k*.2]);
% end
polarplot(ax5b,cat(1,zeros(height(dataTable),1)',deg2rad(dataTable.OptPhase)'),...
    cat(1,zeros(1,height(dataTable),1),dataTable.ECHT_fac'),'o-','Color',co(1,:),'LineWidth',1,'MarkerFaceColor','w');

% Mean setphase 
polarplot(ax5b,[0 circ_mean(deg2rad(dataTable.OptPhase))],[0 max(dataTable.ECHT_fac)],...
    'LineWidth',4,'Color',co(1,:),'LineStyle','-');
ax5b.Title.String = 'P1 Estimated Trough Phase Ranges';

% Imaginary subject
polarplot(ax5b,[deg2rad(180)],[10],'o','Color',[1 0.5 0.8],'MarkerFaceColor',[1 0.5 0.8]);

% Subject 136
polarplot(ax5b,[deg2rad(dataTable(13,:).OptPhase)],[dataTable(13,:).ECHT_fac],'o','Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980]);

% STD
polarplot(ax5b,circ_mean(deg2rad(dataTable.OptPhase))+(linspace(-circ_std(deg2rad(dataTable.OptPhase)),circ_std(deg2rad(dataTable.OptPhase)),100)),...
    max(dataTable.ECHT_fac)*ones(size(linspace(-circ_std(deg2rad(dataTable.OptPhase)),circ_std(deg2rad(dataTable.OptPhase)),100))),...
    'LineWidth',2,'Color',co(1,:));


ax5c = polaraxes('Position',[0.5500 0.0970 0.6975/2 0.7335/2]);
hold on;
minY = min(Y) + 180;
maxY = max(Y) + 180;
% for k = 1:length(alphaCF)
% %     polarplot(ax5b,[0 deg2rad(minY(k))],[0 1],'.--','LineWidth',1);
% %     polarplot(ax5b,[0 deg2rad(maxY(k))],[0 1],'o-','LineWidth',1);
%     polarplot(ax5c,[0 deg2rad(minY(k):maxY(k)) 0],alphaCF(k)*[0 ones(size(minY(k'):maxY(k))),0],...
%         '-','LineWidth',2,'color',[k*.2 0 0]);
% end
polarplot(ax5c,cat(1,zeros(height(dataTable),1)',deg2rad(dataTable.PesPhase)'),...
    cat(1,zeros(1,height(dataTable),1),dataTable.ECHT_fac'),'o-','Color',co(2,:),'LineWidth',1,'MarkerFaceColor','w');

% Mean setphase 
polarplot(ax5c,[0 circ_mean(deg2rad(dataTable.PesPhase))],[0 max(dataTable.ECHT_fac)],...
    'LineWidth',4,'Color',co(2,:),'LineStyle','-');
% STD
polarplot(ax5c,circ_mean(deg2rad(dataTable.PesPhase))+(linspace(-circ_std(deg2rad(dataTable.PesPhase)),circ_std(deg2rad(dataTable.PesPhase)),100)),...
    max(dataTable.ECHT_fac)*ones(size(linspace(-circ_std(deg2rad(dataTable.PesPhase)),circ_std(deg2rad(dataTable.PesPhase)),100))),...
    'LineWidth',2,'Color',co(2,:));

% Imaginary subject
polarplot(ax5c,[deg2rad(0)],[10],'o','Color',[1 0.5 0.8],'MarkerFaceColor',[1 0.5 0.8]);


% Subject 136
polarplot(ax5c,[deg2rad(dataTable(13,:).PesPhase)],[dataTable(13,:).ECHT_fac],'o','Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980]);


ax5c.Title.String = 'P1 Estimated Peak Phase Ranges';

%% Figure 6: Optimal PLV by subject
f5 = figure('Position',[147 377 560 420]);

ax1 = polaraxes;
hold on;

% All
for i = 1:length(OptOn)
polarplot(ax1,[0 OptOn(i)],[0 OptOnPLV(i)],...
    'LineWidth',4,'Color',[co(1,:) 0.1]);
end
% Mean angle and PLV
polarplot(ax1,[0 deg2rad(mAng_OptOn)],[0 PLV_OptOn],...
    'LineWidth',4,'Color',co(1,:));

% SD angle arc
polarplot(ax1,deg2rad(mAng_OptOn)+deg2rad(linspace(-sdAng_OptOn,sdAng_OptOn,100)),...
    PLV_OptOn*ones(size(linspace(-sdAng_OptOn,sdAng_OptOn,100))),...
    'LineWidth',2,'Color',co(1,:));

% Target Onset/Offsets
polarplot(ax1,[0 circ_mean(deg2rad(134))],[0 PLV_OptOn],...
    'LineWidth',1,'Color','k','LineStyle','--');

% Adjust axes properties
ax1.Title.String = sprintf('Alpha Trough PLV (n = %d)',numel(OptOn));
ax1.Title.FontSize = 16;
ax1.RAxisLocation = 190;
ax1.FontSize = 14;



%% Figure 7: Pessimal PLV by subject
f6 = figure('Position',[708 377 560 420]);

ax2 = polaraxes;
hold on;

% All
for i = 1:length(PesOn)
    polarplot(ax2,[0, PesOn(i)],[0 PesOnPLV(i)],...
    'LineWidth',4,'Color',[co(2,:) 0.1]);
end
% Mean angle and PLV
polarplot(ax2,[0 deg2rad(mAng_PesOn)],[0 PLV_PesOn],...
    'LineWidth',4,'Color',co(2,:));

% SD angle arc
polarplot(ax2,deg2rad(mAng_PesOn)+deg2rad(linspace(-sdAng_PesOn,sdAng_PesOn,100)),...
    PLV_PesOn*ones(size(linspace(-sdAng_PesOn,sdAng_PesOn,100))),...
    'LineWidth',2,'Color',co(2,:));

% Target Onset/Offsets
polarplot(ax2,[0 circ_mean(deg2rad(314))],[0 PLV_PesOn],...
    'LineWidth',1,'Color','k','LineStyle','--');

% Adjust axes properties
ax2.Title.String = sprintf('Alpha Peak PLV (n = %d)',numel(PesOn));
ax2.Title.FontSize = 16;
ax2.RAxisLocation = 15;
ax2.FontSize = 14;

%% All phase error - figure 7

f7 = figure('Position',[708 377 560 420]);
nbins = deg2rad(0:20:360); % polar histogram bins
ax7 = polaraxes;
hold on;
polarhistogram(ax7, error,nbins,'Normalization','probability',...
    'FaceAlpha',0.5,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);

% Mean angle and PLV
polarplot(ax7,[0 circ_mean(error)],[0 0.1*PLV_total],...
    'LineWidth',4,'Color',[0.5 0.5 0.5]);

text(ax7,-0.3,0.06,sprintf('PLV: %0.3f',PLV_total), 'HorizontalAlignment', 'center','FontWeight','bold', 'FontSize', 14)
text(-1,0.03,sprintf('PE: %0.3f (\x00B1 %0.3f)',circ_mean(error), circ_std(error)),'FontWeight','bold', 'FontSize', 14)

% Adjust axes properties
ax7.Title.String = sprintf('Post-Hoc Acausal Phase Error');
ax7.Title.FontSize = 16;
ax7.RAxisLocation = 15;
ax7.FontSize = 14;
rlim(ax7, [0 0.1]);

%% All phase error - figure 7
% BIN SIZE PREVIOUSLY 5
f7b = figure('Position',[708 377 560 420]);
nbins = deg2rad(0:20:360); % polar histogram bins
ax7b = polaraxes;
hold on;
polarhistogram(ax7b, echterror ,nbins,'Normalization','probability',...
    'FaceAlpha',0.5,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);

% Mean angle and PLV
polarplot(ax7b,[0 mecht],[0 0.1*PLV_ECHT],...
    'LineWidth',4,'Color',[0.5 0.5 0.5]);

text(ax7b,-0.4,0.06,sprintf('PLV: %0.3f',PLV_ECHT), 'HorizontalAlignment', 'center','FontWeight','bold', 'FontSize', 14)
text(-1.1,0.035,sprintf('PE: %0.3f (\x00B1 %0.3f)',circ_mean(echterror), circ_std(echterror)),'FontWeight','bold', 'FontSize', 14)

% Adjust axes properties
ax7b.Title.String = sprintf('ECHT Phase Error');
ax7b.Title.FontSize = 16;
ax7b.RAxisLocation = 15;
ax7b.FontSize = 14;
rlim(ax7b, [0 0.1]);

%% IAF to Phase error

% Get colors
color = [0.0000 0.4470 0.7410
    0.8500 0.10 0.0980];

f8 = figure;
subplot(2,2,1)
mdl = fitlm(iaf, accuracy);
coefs = mdl.Coefficients.Estimate;
scatter(iaf, accuracy, 70, color(1,:), 'filled')
hold on
xlim([8 13])
r = refline(coefs(2),coefs(1)); 
r.Color = color(1,:);
r.LineWidth = 2;
xlabel('IAF')
x = linspace(8, 13, 1000)';
[~, ci] = predict(mdl, x);
fill([x; flipud(x)],[ci(:,1); flipud(ci(:,2))], color(1,:),'edgecolor','none','facealpha', 0.1)
ylabel('Acausal Phase-locking Accuracy')
title('IAF and Post-hoc Acausal PLV')
grid on
axis square
box off
text(11.5, 0.6,...
    sprintf('Adjusted R-squared: %0.2f', mdl.Rsquared.Adjusted),"HorizontalAlignment","center",...
    'FontSize', 12)
text(11.5, 0.56,...
    sprintf('P-value: %0.3f', mdl.ModelFitVsNullModel.Pvalue),"HorizontalAlignment","center",...
    'FontSize', 12)


subplot(2,2,2)
mdl = fitlm(iaf, echt);
coefs = mdl.Coefficients.Estimate;
scatter(iaf, echt, 70, color(1,:), 'filled')
hold on
xlim([8 13])
r = refline(coefs(2),coefs(1)); 
r.Color = color(1,:);
r.LineWidth = 2;
xlabel('IAF')
x = linspace(8, 13, 1000)';
[~, ci] = predict(mdl, x);
fill([x; flipud(x)],[ci(:,1); flipud(ci(:,2))], color(1,:),'edgecolor','none','facealpha', 0.1)
ylabel('ECHT Phase-locking Accuracy')
title('IAF and ECHT PLV')
grid on
axis square
box off
text(11.5, 0.97,...
    sprintf('Adjusted R-squared: %0.2f', mdl.Rsquared.Adjusted),"HorizontalAlignment","center",...
    'FontSize', 12)
text(11.5, 0.96,...
    sprintf('P-value: %0.3f', mdl.ModelFitVsNullModel.Pvalue),"HorizontalAlignment","center",...
    'FontSize', 12)

subplot(2,2,3)
mdl = fitlm(iaf, phaselag);
coefs = mdl.Coefficients.Estimate;
scatter(iaf, phaselag, 70, color(1,:), 'filled')
hold on
xlim([8 13])
r = refline(coefs(2),coefs(1)); 
r.Color = color(1,:);
r.LineWidth = 2;
xlabel('IAF')
x = linspace(8, 13, 1000)';
[~, ci] = predict(mdl, x);
fill([x; flipud(x)],[ci(:,1); flipud(ci(:,2))], color(1,:),'edgecolor','none','facealpha', 0.1)
ylabel('Post-hoc Phase Error (rad)')
title('IAF and Post-hoc Acausal Phase Error')
grid on
axis square
box off
text(11.5, 1.5,...
    sprintf('Adjusted R-squared: %0.2f', mdl.Rsquared.Adjusted),"HorizontalAlignment","center",...
    'FontSize', 12)
text(11.5, 1.25,...
    sprintf('P-value: %0.3f', mdl.ModelFitVsNullModel.Pvalue),"HorizontalAlignment","center",...
    'FontSize', 12)


subplot(2,2,4)
mdl = fitlm(iaf, echtlag);
coefs = mdl.Coefficients.Estimate;
scatter(iaf, echtlag, 70, color(1,:), 'filled')
hold on
xlim([8 13])
r = refline(coefs(2),coefs(1)); 
r.Color = color(1,:);
r.LineWidth = 2;
xlabel('IAF')
x = linspace(8, 13, 1000)';
[~, ci] = predict(mdl, x);
fill([x; flipud(x)],[ci(:,1); flipud(ci(:,2))], color(1,:),'edgecolor','none','facealpha', 0.1)
ylabel('ECHT Phase Error (rad)')
title('IAF and ECHT Phase Error')
grid on
axis square
box off
text(11.5, 0.05,...
    sprintf('Adjusted R-squared: %0.2f', mdl.Rsquared.Adjusted),"HorizontalAlignment","center",...
    'FontSize', 12)
text(11.5, 0,...
    sprintf('P-value: %0.3f', mdl.ModelFitVsNullModel.Pvalue),"HorizontalAlignment","center",...
    'FontSize', 12)

%% Plot Power

% Signifigance Testing
[pitc,~, ~] = npperm(squeeze(itc(medis{2},:,:,1)),...
              squeeze(itc(medis{2},:,:,2)), 2000,"cluster", 0, 0);
[pevo,~, ~] = npperm(squeeze(evo(medis{2},:,:,1)),...
              squeeze(evo(medis{2},:,:,2)), 2000,"cluster", 0, 0);
[ptot,~, ~] = npperm(squeeze(tot(medis{2},:,:,1)),...
              squeeze(tot(medis{2},:,:,2)), 2000,"cluster", 0, 0);
[pind,~, ~] = npperm(squeeze(ind(medis{2},:,:,1)),...
              squeeze(ind(medis{2},:,:,2)), 2000,"cluster", 0, 0);

% ITPC Figures
f9 = figure;
for j = 2:2
    plot(t, itcmu(:,1,1,j),'color',0.5*j*co(1,:),'linewidth',2)
    hold on
    fill(time,itcse(:,1,1,j),0.5*j*co(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    plot(t, itcmu(:,1,2,j),'color',0.5*j*co(2,:),'linewidth',2)
    fill(time,itcse(:,1,2,j),0.5*j*co(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    xlim([0 max(t)])
    ylim([0 0.7])
    xlabel('Time (ms)')
    ylabel('ITPC')
    axis square
    box off
    title('ITPC FPZ')

    if j == 2
        [x, y] = clusts4plot(pitc(:,1), 0.05, [-10 10]);
        for i = 1:size(x,1)
            patch(t(x(i,:)), y, [0.5 0.5 0.5],'LineStyle','none','FaceAlpha',0.05);
            hold on
        end
    end
end

f10 = figure;
for j = 2:2
    plot(t, itcmu(:,2,1,j),'color',0.5*j*co(1,:),'linewidth',2)
    hold on
    fill(time,itcse(:,2,1,j),0.5*j*co(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    plot(t, itcmu(:,2,2,j),'color',0.5*j*co(2,:),'linewidth',2)
    fill(time,itcse(:,2,2,j),0.5*j*co(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    xlim([0 max(t)])
    ylim([0 0.7])
    xlabel('Time (ms)')
    ylabel('ITPC')
    axis square
    box off
    title('ITPC FZ')

    if j == 2
        [x, y] = clusts4plot(pitc(:,2), 0.05, [-10 10]);
        for i = 1:size(x,1)
            patch(t(x(i,:)), y, [0.5 0.5 0.5],'LineStyle','none','FaceAlpha',0.05);
            hold on
        end
    end
end


% Evoked Power Figures
f11 = figure;
for j = 2:2
    plot(t, evomu(:,1,1,j),'color',0.5*j*co(1,:),'linewidth',2)
    hold on
    fill(time,evose(:,1,1,j),0.5*j*co(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    plot(t, evomu(:,1,2,j),'color',0.5*j*co(2,:),'linewidth',2)
    fill(time,evose(:,1,2,j),0.5*j*co(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    xlim([0 max(t)])
    xlabel('Time (ms)')
    ylabel('Evoked Power')
    axis square
    box off
    title('Evoked FPZ')

    if j == 2
        [x, y] = clusts4plot(pevo(:,1), 0.05, [-10 10]);
        for i = 1:size(x,1)
            patch(t(x(i,:)), y, [0.5 0.5 0.5],'LineStyle','none','FaceAlpha',0.05);
            hold on
        end
    end
end

f12 = figure;
for j = 2:2
    plot(t, evomu(:,2,1,j),'color',0.5*j*co(1,:),'linewidth',2)
    hold on
    fill(time,evose(:,2,1,j),0.5*j*co(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    plot(t, evomu(:,2,2,j),'color',0.5*j*co(2,:),'linewidth',2)
    fill(time,evose(:,2,2,j),0.5*j*co(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    xlim([0 max(t)])
    xlabel('Time (ms)')
    ylabel('Evoked Power')
    axis square
    box off
    title('Evoked FZ')

    if j == 2
        [x, y] = clusts4plot(pevo(:,2), 0.05, [-10 10]);
        for i = 1:size(x,1)
            patch(t(x(i,:)), y, [0.5 0.5 0.5],'LineStyle','none','FaceAlpha',0.05);
            hold on
        end
    end
end


% Total Power Figures
f13 = figure;
for j = 2:2
    plot(t, totmu(:,1,1,j),'color',0.5*j*co(1,:),'linewidth',2)
    hold on
    fill(time,totse(:,1,1,j),0.5*j*co(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    plot(t, totmu(:,1,2,j),'color',0.5*j*co(2,:),'linewidth',2)
    fill(time,totse(:,1,2,j),0.5*j*co(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    xlim([0 max(t)])
    xlabel('Time (ms)')
    ylabel('Total Power')
    axis square
    box off
    title('Total FPZ')

    if j == 2
        [x, y] = clusts4plot(ptot(:,1), 0.05, [-10 10]);
        for i = 1:size(x,1)
            patch(t(x(i,:)), y, [0.5 0.5 0.5],'LineStyle','none','FaceAlpha',0.05);
            hold on
        end
    end
end

f14 = figure;
for j = 2:2
    plot(t, totmu(:,2,1,j),'color',0.5*j*co(1,:),'linewidth',2)
    hold on
    fill(time,totse(:,2,1,j),0.5*j*co(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    plot(t, totmu(:,2,2,j),'color',0.5*j*co(2,:),'linewidth',2)
    fill(time,totse(:,2,2,j),0.5*j*co(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    xlim([0 max(t)])
    xlabel('Time (ms)')
    ylabel('Total Power')
    axis square
    box off
    title('Total FZ')

    if j == 2
        [x, y] = clusts4plot(ptot(:,2), 0.05, [-10 10]);
        for i = 1:size(x,1)
            patch(t(x(i,:)), y, [0.5 0.5 0.5],'LineStyle','none','FaceAlpha',0.05);
            hold on
        end
    end
end



% Induced Power Figures
f15 = figure;
for j = 2:2
    plot(t, indmu(:,1,1,j),'color',0.5*j*co(1,:),'linewidth',2)
    hold on
    fill(time,indse(:,1,1,j),0.5*j*co(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    plot(t, indmu(:,1,2,j),'color',0.5*j*co(2,:),'linewidth',2)
    fill(time,indse(:,1,2,j),0.5*j*co(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    xlim([0 max(t)])
    xlabel('Time (ms)')
    ylabel('Induced Power')
    axis square
    box off
    title('Induced FPZ')

    if j == 2
        [x, y] = clusts4plot(pind(:,1), 0.05, [-10 10]);
        for i = 1:size(x,1)
            patch(t(x(i,:)), y, [0.5 0.5 0.5],'LineStyle','none','FaceAlpha',0.05);
            hold on
        end
    end
end


f16 = figure;
for j = 2:2
    plot(t, indmu(:,2,1,j),'color',0.5*j*co(1,:),'linewidth',2)
    hold on
    fill(time,indse(:,2,1,j),0.5*j*co(1,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    plot(t, indmu(:,2,2,j),'color',0.5*j*co(2,:),'linewidth',2)
    fill(time,indse(:,2,2,j),0.5*j*co(2,:), ...
        'edgecolor','none', ...
        'facealpha', 0.1)
    xlim([0 max(t)])
    xlabel('Time (ms)')
    ylabel('Induced Power')
    axis square
    box off
    title('Induced FZ')

    if j == 2
        [x, y] = clusts4plot(pind(:,2), 0.05, [-10 10]);
        for i = 1:size(x,1)
            patch(t(x(i,:)), y, [0.5 0.5 0.5],'LineStyle','none','FaceAlpha',0.05);
            hold on
        end
    end
end


%% Phase Error Final Plot

figure;
set(0,'DefaultLineLineWidth',2)
h = cdfplot(abs(wrapTo180(rad2deg(phaselag)))); % swap for error
h.Color = co(1,:);
hold on
h1 = cdfplot(abs(wrapTo180(rad2deg(echtlag)))); % swap for echterror
h1.Color = [0.5 0.5 0.5];
plot([0 180],[0 1],'color',co(2,:));
title('Cumualtive Distribution of Average Phase Error');
ylabel('Cumulative Probability')
xlabel('Phase Error in degrees')
legend({'Acausal Post-Hoc','ecHT Causal PE',''})
box off
axis square

figure;
set(0,'DefaultLineLineWidth',2)
h = cdfplot(abs(wrapTo180(rad2deg(error)))); % swap for error
h.Color = co(1,:);
hold on
h1 = cdfplot(abs(wrapTo180(rad2deg(echterror)))); % swap for echterror
h1.Color = [0.5 0.5 0.5];
plot([0 180],[0 1],'color',co(2,:));
title('Cumualtive Distribution of Average Phase Error');
ylabel('Cumulative Probability')
xlabel('Phase Error in degrees')
legend({'Acausal Post-Hoc','ecHT Causal PE',''})
box off
axis square

%% Save
% % save("PLV.mat")
svpath = '/Applications/Toolbox/MATLAB/Figures/Polar Plots';
print(f1,fullfile(svpath,'Trough Phases.svg'),'-dsvg');
print(f2,fullfile(svpath,'Peak Phases.svg'),'-dsvg');
print(f3,fullfile(svpath,'Group Post-hoc Phase Error.svg'),'-dsvg');
print(f3b,fullfile(svpath,'Group ECHT Phase Error.svg'),'-dsvg');
print(f4,fullfile(svpath,'IAF, Latency, Peak & Trough Phase.svg'),'-dsvg');
print(f5, fullfile(svpath, 'Trough PLVs.svg'),'-dsvg');
print(f6, fullfile(svpath, 'Peak PLVs.svg'),'-dsvg');
print(f7, fullfile(svpath, 'All Post-hoc Phase Error.svg'),'-dsvg');
print(f7b, fullfile(svpath, 'All ECHT Phase Error.svg'),'-dsvg');
print(f8, fullfile(svpath, 'IAF and PLV.svg'),'-dsvg');
print(f9, fullfile(svpath,'ITPC FPZ.svg'),'-dsvg');
print(f10, fullfile(svpath,'ITPC FZ.svg'),'-dsvg');
print(f11, fullfile(svpath,'Evoked Power FPZ.svg'),'-dsvg');
print(f12, fullfile(svpath,'Evoked Power FZ.svg'),'-dsvg');
print(f13, fullfile(svpath,'Total Power FPZ.svg'),'-dsvg');
print(f14, fullfile(svpath,'Total Power FZ.svg'),'-dsvg');
print(f15, fullfile(svpath,'Induced Power FPZ.svg'),'-dsvg');
print(f16, fullfile(svpath,'Induced Power FZ.svg'),'-dsvg');
print(f1,fullfile(svpath,'Trough Phases.png'),'-dpng');
print(f2,fullfile(svpath,'Peak Phases.png'),'-dpng');
print(f3,fullfile(svpath,'Group Post-hoc Phase Error.png'),'-dpng');
print(f3b,fullfile(svpath,'Group ECHT Phase Error.png'),'-dpng');
print(f4,fullfile(svpath,'IAF, Latency, Peak & Trough Phase.png'),'-dpng');
print(f5, fullfile(svpath, 'Trough PLVs.png'),'-dpng');
print(f6, fullfile(svpath, 'Peak PLVs.png'),'-dpng');
print(f7, fullfile(svpath, 'All Post-hoc Phase Error.png'),'-dpng');
print(f7b, fullfile(svpath, 'All ECHT Phase Error.png'),'-dpng');
print(f8, fullfile(svpath, 'IAF and PLV.png'),'-dpng');
print(f9, fullfile(svpath,'ITPC FPZv'),'-dpng');
print(f10, fullfile(svpath,'ITPC FZ.png'),'-dpng');
print(f11, fullfile(svpath,'Evoked Power FPZ.png'),'-dpng');
print(f12, fullfile(svpath,'Evoked Power FZ.png'),'-dpng');
print(f13, fullfile(svpath,'Total Power FPZ.png'),'-dpng');
print(f14, fullfile(svpath,'Total Power FZ.png'),'-dpng');
print(f15, fullfile(svpath,'Induced Power FPZ.png'),'-dpng');
print(f16, fullfile(svpath,'Induced Power FZ.png'),'-dpng');

