function params = echtparams(varargin)

% This function loads basic parameters for processing the ECHT output,
% unless otherwise specified.  Requires the finputcheck function from
% EEGLAB
%
% Written Tylor J Harlow 9/12/2023 - tylor.harlow@uconn.edu


params = finputcheck(varargin,...
    {
        % LabView streaming command
        'strcmd' 'string' 'str,m,dt,eeg1,eeg2,out,iphs,ps_pulse_count,ps_pulse_id,ps_onphs,ps_offphs,ps_onontime,ps_durbyphs,ps_durbytime,ps_isi'...
                           'str,m,dt,eeg1,eeg2,out,iphs,ps_pulse_count,ps_pulse_id,ps_onphs,ps_offphs,ps_onontime,ps_durbyphs,ps_durbytime,ps_isi';
        
        % Define pre-processing bandpass filter for raw EEG data
        'ord'                      'real'         []                                     4; % Filter Order
        'taps'                     'real'         [0 inf]                             1024; % number of taps for pre-processing FIR filter
        'gd1'                      'real'         [0 inf]                           1024/2;
        
        % Wavelets
        'cycles'                   'real'         []                                [4 12]; % Wavelet cycles
        'freqs'                    'real'         []                                [1 70]; % Wavelet frequency range
        'steps'                    'real'         []                                   140; % Numer of wavelets

        % Frequency sliding
        'detrend'                  'real'         []                                     1; % Whether or not to detrend FFT before

        % Define passbands for filters
        'fpass1'                   'real'         []                                [2 34]; % pre-processing passband
        'fpassA'                   'real'         []                        10*[0.75 1.25]; % alpha passband
        'useECHTAlpha'             'real'         [0 inf]                                0; % Alpha center frequency defined by ECHT device? [0:no, 1:yes]
        'thetaCF'                  'real'         [0 inf]                                5;
        'fpassT'                   'real'         []                         5*[0.75 1.25]; % theta passband
        
        % Define Analysis and Quality parameters
        'tWin'                     'real'         []                        [-0.250 0.500]; % default;
        'analysisWin'              'real'         []                               [0 200]; % Post-stimulus analysis comparison window (ms)
        'thresh'                   'real'         [0 inf]                              100;
        'filterImplementation'     'real'         [0 inf]                                1; %[0:unfiltered, 1:zero-phase FIR, 2: causal IIR]
        'zeroPhase'                'real'         [0 inf]                                1; % Use zero-phase bandpass filtering for alpha/theta? [0:no 1:yes]
        'BaseLineCorrection'       'real'         [0 inf]                                0; % Apply pre-stimulus basline correction to ERPs? [0:no, 1:yes]

    }, 'echtparams');

end