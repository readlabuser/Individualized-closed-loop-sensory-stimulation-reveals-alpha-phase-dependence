function [ERP,params] = batchERP(pname,params)
useZeroPhaseFiltering = params.zeroPhase;

%% Initialize variables
ERP.eeg = [];
ERP.erp = [];
ERP.stim = [];
ERP.echtPhase = [];
ERP.echtPhaseLabel = [];
ERP.aPhaseOn = [];
ERP.aPhaseOf = [];
ERP.aPhase = [];
ERP.tPhase = [];
ERP.erpStim = [];
ERP.halpha = [];
ERP.walpha = [];
ERP.isi = [];
ERP.waPhaseOn = [];

%% Read .csv ECHT log file
c = dir(fullfile(pname,'*.csv'));
fnames = {c.name}';

%% Iterate through .csv log files and combine
for k = 1:length(fnames)
    [rawData,Echo,opts] = readECHTLog(fullfile(pname,fnames{k}),params.strcmd);
    
    data(k).t = double(rawData.Time_us-rawData.Time_us(1))*1e-6; % time vector
    data(k).fs = round(1/mean(rawData.dt*1e-6)); % average sampling frequency
    ERP.fs(k) = data(k).fs;
    data(k).rawEEG = [rawData.eeg1,rawData.eeg2]*1e6; % scalp EEG voltages
    %     data(k).rawEEG = randn(size(data(k).rawEEG)); % random noise control (delete later)
    data(k).stim = rawData.out/max(rawData.out); % pulse train gating stimulus [0,1]
    data(k).pID = rawData.ps_pulse_count;

    ERP.eeg{k} = data(k).rawEEG;
    
    % Grab ECHT-set alpha center frequency setting !fac,##
    % params.alphaCF = str2double(Echo.ExtraVar2(contains(Echo.cmd,'!fac')))
    params.fpassA = params.alphaCF*[0.75 1.25];
    % params.fpassA = params.alphaCF*[0.80 1.20];
    params.useECHTAlpha = 1;
    params.fs = ERP.fs(k);
    
    % ECHT set phase(s)
    idx = find(contains(Echo.cmd,'!ps_onphs'));
    onPh = str2double(char(Echo{idx,2})):str2double(char(Echo{idx,3})):str2double(char(Echo{idx,4}));
    
    % Bandpass filter for pre-processing of raw EEG data
    fprintf('Bandpass filtering raw EEG signal to %d-%d Hz...',params.fpass1);
    switch params.filterImplementation
        case 0
            data(k).eeg = data(k).rawEEG;
        case 1
            Hbp = kaiserBPF(params);
            data(k).eeg = filtfilt(Hbp.Numerator,1,data(k).rawEEG);
        case 2
            Hbp = getBPFilter;
            set(Hbp,'PersistentMemory',true); 
%             data(k).eeg = filter(Hbp,data(k).rawEEG);
            data(k).eeg = circshift(filter(Hbp,data(k).rawEEG),-62,1); % adjust for group delay
            data(k).stim = circshift(data(k).stim,-62);
    end
    fprintf('COMPLETE\n');

    % Identify pulse onsets and divide continuous BPF EEG data into epochs
    fprintf('Epoch continuous EEG data...');
    kOn = cat(1,false,diff(data(k).stim)==1);
    kOf = cat(1,false,diff(data(k).stim)==-1);
    kOn1 = kOn & rawData.ps_pulse_count==1;
    kOf1 = kOf & rawData.ps_pulse_count==1;
    data(k).kOn = find(kOn);
    data(k).kOn1 = find(kOn1);
    data(k).kOf1 = find(kOf1);
    
    data(k).time{k} = zeros(size(data(k).eeg));
    data(k).time{k}(kOn1) = 1;
    
    % Define ERP analysis window
    data(k).params.tWin = params.tWin; % analysis window (seconds)
    kWin = round(data(k).params.tWin*data(k).fs);
    t0 = ([kWin(1):kWin(2)]'/data(k).fs)*1e3;
    
    % Check for epoching time overun and adjust
    if (data(k).kOn1(end)+kWin(2)>length(data(k).eeg))
        warning('Epoching window overun: removing last onset event');
        kOn1(data(k).kOn1(end)) = false;
        data(k).kOn1 = find(kOn1);
    end

    % Bandpass filter for alpha band (7.5-12.5 Hz)
    fprintf('Bandpass filtering raw EEG signal to alpha %1.2f-%1.2f Hz ',params.fpassA);
    [bA,aA] = butter(params.ord,params.fpassA/(data(k).fs/2),'bandpass');
    [gdA,fr] = grpdelay(bA,aA,1024,data(k).fs);
    if useZeroPhaseFiltering
        fprintf('(zero-phase)...');
        data(k).alphaEEG = filtfilt(bA,aA,data(k).rawEEG);
    else
        fprintf('(causal IIR)...');
%         mgdA = round(mean(gdA(fr>=params.fpassA(1) & fr<=params.fpassA(2))));
%         mgdA = 102; % fudge factor for alpha group delay
%         tmp = filter(bA,aA,data(k).rawEEG);
%         data(k).alphaEEG = cat(1,tmp(mgdA+1:end,:),zeros(mgdA,size(tmp,2)));
        data(k).alphaEEG = filter(bA,aA,data(k).rawEEG);
    end
    hAlpha = hilbert(data(k).alphaEEG);
    data(k).alphaPhaseOn = angle(hAlpha(data(k).kOn1(1:end),:));
    data(k).alphaPhaseOf = angle(hAlpha(data(k).kOf1(1:end),:));
    fprintf('COMPLETE\n');
    
    % Bandpass fitler for theta band (4-7 Hz)
    fprintf('Bandpass filtering raw EEG signal to theta %1.2f-%1.2f Hz ',params.fpassT);
    [bT,aT] = butter(params.ord,params.fpassT/(data(k).fs/2),'bandpass');
    [gdT,fr] = grpdelay(bT,aT,1024,data(k).fs);
    if useZeroPhaseFiltering
        fprintf('(zero-phase)...');
        data(k).thetaEEG = filtfilt(bT,aT,data(k).rawEEG);
    else
        fprintf('(causal IIR)...');
%         mgdT = round(mean(gdT(fr>=params.fpassT(1) & fr<=params.fpassT(2))));
%         mgdT = 205; % fudge factor for theta group delay
%         tmp = filter(bT,aT,data(k).rawEEG);
%         data(k).thetaEEG = cat(1,tmp(mgdT+1:end,:),zeros(mgdT,size(tmp,2)));
        data(k).thetaEEG = filter(bT,aT,data(k).rawEEG);
    end
    hTheta = hilbert(data(k).thetaEEG);
    data(k).thetaPhase = angle(hTheta(data(k).kOn1(1:end),:));
    fprintf('COMPLETE\n');

    % Wavelet alpha
    c = linspace(4, 12, 140);
    f = linspace(1, 70, 140);
    [~,idx] = min(abs(params.alphaCF - f)); 
    [waves,~] = wavecwt(data(k).eeg, data(k).fs, params.alphaCF, c(idx), 1);

    % Check wavelet phases
    data(k).walphaPhaseOn = angle(waves(data(k).kOn1(1:end),:));

    m = 0;
    %     for n = 1:length(data(k).kOn)
    for n = 1:length(data(k).kOn1)
        m = m+1;
        try
            data(k).erp(:,:,m) = data(k).eeg(data(k).kOn1(n)+[kWin(1):kWin(2)],:);
            data(k).halpha(:,:,m) = hAlpha(data(k).kOn1(n)+[kWin(1):kWin(2)],:);
            data(k).walpha(:,:,m) = waves(data(k).kOn1(n)+[kWin(1):kWin(2)],:);
            data(k).erpStim(:,m) = data(k).stim(data(k).kOn1(n)+[kWin(1):kWin(2)]);
        catch
            fprintf('Whoa\n');
            data(k).kOn = data(k).kOn1(1:end-1);
            break;
        end
    end
    fprintf('COMPLETE\n');
  
    % Get interstimulus interval
    n = diff(data(k).kOn1) * (1/data(k).fs);
    
    % Grab instantaneous phase information from ECHT device
%     ERP.echtPhase
    
    % Combine epochs across session
    ERP.isi = cat(1, ERP.isi, n);
    ERP.time{k} = data(k).time;
    ERP.erp = cat(3,ERP.erp,data(k).erp);
    ERP.halpha = cat(3, ERP.halpha, data(k).halpha);
    ERP.walpha = cat(3, ERP.walpha, data(k).walpha);
    ERP.stim{k} = data(k).stim;
    ERP.tPhase = cat(1,ERP.tPhase,data(k).thetaPhase);
    ERP.echtPhase = cat(1,ERP.echtPhase,wrapToPi(rawData.iphs(kOn1)));
    ERP.echtPhaseLabel = cat(1,ERP.echtPhaseLabel,rawData.ps_onphs(kOn1));
    ERP.erpStim = cat(2,ERP.erpStim,data(k).erpStim);
    ERP.alphaCF = params.alphaCF;
    ERP.setPhases = onPh;
    ERP.time{k} = data(k).t;
    ERP.aPhaseOn = cat(1,ERP.aPhaseOn,data(k).alphaPhaseOn);
    ERP.waPhaseOn = cat(1,ERP.waPhaseOn,data(k).walphaPhaseOn);
    ERP.aPhaseOf = cat(1,ERP.aPhaseOf,data(k).alphaPhaseOf);
    fprintf('COMPLETE\n');
    
end
ERP.t0 = t0;
ERP.nSessions = length(fnames);

% Find epoch-free (good) and artifact-contaminated (rej) trials
ERP.goodTrials = find(max(max(abs(ERP.erp)),[],2)<params.thresh);
ERP.rejTrials = find(max(max(abs(ERP.erp)),[],2)>=params.thresh);

%% Phase sort epochs to alpha and theta
PHIs = (-180:30:180)-15;

tIDK{1} = find(ERP.tPhase(:,1) >=deg2rad(165) & ERP.tPhase(:,1)<=deg2rad(180) | ERP.tPhase(:,1)>=(-180) & ERP.tPhase(:,1)<=deg2rad(-165));
aIDK{1} = find(ERP.aPhaseOn(:,1) >=deg2rad(165) & ERP.aPhaseOn(:,1)<=deg2rad(180) | ERP.aPhaseOn(:,1)>=(-180) & ERP.aPhaseOn(:,1)<=deg2rad(-165));
eIDK{1} = find(ERP.echtPhase(:,1) >=deg2rad(165) & ERP.echtPhase(:,1)<=deg2rad(180) | ERP.echtPhase(:,1)>=(-180) & ERP.echtPhase(:,1)<=deg2rad(-165));

for k = 2:length(PHIs)-1
    tIDK{k} = find(ERP.tPhase(:,1)>=deg2rad(PHIs(k)) & ERP.tPhase(:,1)<=deg2rad(PHIs(k+1)));
    aIDK{k} = find(ERP.aPhaseOn(:,1)>=deg2rad(PHIs(k)) & ERP.aPhaseOn(:,1)<=deg2rad(PHIs(k+1)));
    eIDK{k} = find(ERP.echtPhase(:,1)>=deg2rad(PHIs(k)) & ERP.echtPhase(:,1)<=deg2rad(PHIs(k+1)));
end

ERP.thetaByPhase = tIDK;
ERP.alphaByPhase = aIDK;
ERP.echtByPhase = eIDK;
ERP.phi = PHIs+15;

end % EOF