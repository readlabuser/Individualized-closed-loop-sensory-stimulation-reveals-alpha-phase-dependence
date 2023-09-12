function [dataTable,echoTable,opts] = readECHTLog(filename,strcmd)
% [dataTable,echoTable] = readECHTLog(filename,strcmd)
%
% Parses the ECHT .csv log file into two different MATLAB tables,
% "dataTable" and "echoTable".
%
% readECHTLog properly assigns variable names and types to the correct
% columns based on the !str command line.
%
% This function is firmware-specific (2020-12-25)
% Use for "TwoPulseSleep" and "TwoPulseAwake" experiments.
%
% INPUT VARIABLES
%    filename : full path and file name of the .csv log file
%      strcmd : !str command line from LabView script (optional input)
%
% OUTPUT VARIABLES
%   dataTable : ECHT data table
%   echoTable : command echo table
%
% Created 2020-12-29 scott.bressler@elemindtech.com

%% Define table variables
if nargin<2
    strcmd = 'str,m,dt,eeg1,eeg2,out,iphs,ps_pulse_count,ps_pulse_id,ps_onphs,ps_offphs,ps_onontime,ps_durbyphs,ps_durbytime,ps_isi';
end

str = strsplit(strcmd,','); %%% Switched????

varNames = {'DateTime','Time_us','Type','Mode'};
varNames = cat(2,varNames,str(3:end));

opts = delimitedTextImportOptions(...
    'NumVariables',32,...
    'VariableNames',varNames,...
    'Delimiter',',');

% Define variable types for log data
opts.VariableTypes([1,3]) = {'char'};
opts.VariableTypes([2,5,6,7,9,10:13]) = {'double'};
opts.VariableTypes([4,8,14:17]) = {'int16'};

%% Load log data and sub-divide table into data and echo
fprintf('Loading .csv log data...');
rawData = readtable(filename,opts);

dataTable = rawData(strcmp(rawData.Type,'data0'),1:length(varNames));
echoTable = rawData(strcmp(rawData.Type,'echo'),18:end);
echoTable.Properties.VariableNames(1) = {'cmd'};

fprintf('COMPLETE\n');



