function [tf, frex] = wavecwt(signal, srate, frange, ncycles, nsteps)

% This function performs a highly customizable complex morlet wavelet transform.
% Written Tylor J Harlow 9/12/2023 - tylor.harlow@uconn.edu
%
% adapted from code by Michael X Cohen
% "A better way to define and describe Morlet wavelets for time-frequency analysis,
% NeuroImage,
% Volume 199,
% 2019,
% Pages 81-86,
% ISSN 1053-8119,
% https://doi.org/10.1016/j.neuroimage.2019.05.048.
% (https://www.sciencedirect.com/science/article/pii/S1053811919304409)"


% Inputs
%
%     signal  = recording in the form [time x trial/channel]
%     srate   = sampling rate
%     nfrex   = number of center frequencies 
%     ncycles = number of cycles specified as one or two digit vector
%     frange  = frequency range specied as a vector of two digits
% Outputs
%
%     tf      = time-frequency response dimensions [trial x f x time]


if length(ncycles) <2
    ncycles = [ncycles ncycles];
end

if isempty(nsteps)
    nfrex = 40;
end

if length(frange) <2
    frange = [frange frange];
end

if length(ncycles) <2
    ncycles = [ncycles ncycles];
end

%%% time-frequency analysis parameters
frex = linspace(frange(1),frange(2),nsteps);
cycles = linspace(ncycles(1),ncycles(2),nsteps);


% setup wavelet and convolution parameters
pnts  = length(signal(:,1));
trials  = length(signal(1,:));
wavet = [-1 : 1/srate : 1];
halfw = (length(wavet) - 1)/2;
nWave = length(wavet);
nData = pnts*trials;
nConv = nWave + nData - 1;

% initialize time-frequency matrix
tf = zeros(pnts,length(signal(1,:)),length(frex));


% loop over frequencies
for fi=1:length(frex)
    % Cycles & Wavelet
    s = cycles(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*wavet) .* exp(-wavet.^2./(2*s^2));

    % Convolvs
    waveX = fft( wavelet, nConv);
    waveX = waveX./max(real(waveX));
    
    % Reshape signal
    timeS  = 0:1/srate:(1/frex(fi));
    sig = reshape(signal,1,[]);

    % spectrum of data
    dataX = fft( sig, nConv);
    as = ifft( waveX.*dataX, nConv);
    
    % trim and reshape
    as = as(1:nConv);
    as = as(halfw + 1 : end - halfw);
    as = reshape(as,pnts,1,trials);
    tf(:, :, fi) = as;

end
end
