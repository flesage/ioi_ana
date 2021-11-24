function OutData = NormalisationFiltering(FolderData, FileData, lowFreq, highFreq, bDivide, varargin)
%%%%% Data Normalisation by low-pass filtering %%%%%
% 
% General Infos:
%
% This function can be used to normalise channels (delta F/F or delta R/R),
% or to do a low-pass filtering.
%
% Inputs:
%
% Option A: data to be normalised will be opened within this function
%
%   1- FolderData:  Folder containing the data to be oppened
%   2- FileData:    Channel to open('red', 'green', 'fluo_475', etc.)
%   3- lowFreq:     low frequency cut-off, set to 0 to ignore
%   4- highFreq:    high frequency cut-off, set to 0 to ignore
%   5- bDivide:     if 1, the data returned (below highFreq) is normalised
%                   by the low freq signal (below lowFreq)
%                   if 0, the low freq signal (below lowFreq) is
%                   substracted from the data returned (below highFreq)
%
%   Ex: Dat = NormalisationFiltering(pwd, 'red', 1/120, 1, 1);
%   
%   This call would return the red channel with a low-pass at 1 Hz,
%   normalized (through a division) by a low-pass at 1/120 Hz.
%
% Option B: data to be normalised is given by one of the argument
%
%   1- FolderData:  Folder containing the data to be oppened
%   2- FileData:    Data, as a 3D matrix (Y, X, Time).
%   3- lowFreq:     low frequency cut-off, set to 0 to ignore
%   4- highFreq:    high frequency cut-off, set to 0 to ignore 
%   5- bDivide:     if 1, the data returned (below highFreq) is normalised
%                   by the low freq signal (below lowFreq)
%                   if 0, the low freq signal (below lowFreq) is
%                   substracted from the data returned (below highFreq)
%   6- Freq (Optional): Data sample rate.
%
%   Ex: Dat = NormalisationFiltering(pwd, dat, 0, 1, 1);
%   
%   This call would return the data contained in the variable dat with a
%   low-pass at 1 Hz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'FolderData',@isfolder)
addRequired(p,'FileData', @(x) ischar(x) || isnumeric(x))
addRequired(p,'lowFreq',@(x) isscalar(x) & isnumeric(x))
addRequired(p,'highFreq', @(x) isscalar(x) & isnumeric(x))
addRequired(p,'bDivide', @(x) isscalar(x) & (isnumeric(x) | islogical(x)))
addOptional(p,'Freq', [], @(x) (isscalar(x) & isnumeric(x)) | isempty(x))
% Parse inputs:
parse(p, FolderData, FileData, lowFreq, highFreq, bDivide, varargin{:});
% Select function based on "FileData" variable type:
if( isa(p.Results.FileData, 'char') )
    OutData = NormFiltFromFile(p.Results.FolderData, p.Results.FileData,...
        p.Results.lowFreq, p.Results.highFreq, p.Results.bDivide);
else
    OutData = NormFiltDirect(p.Results.FolderData, p.Results.FileData,...
        p.Results.lowFreq, p.Results.highFreq, p.Results.bDivide, p.Results.Freq);
end

end


function OutData = NormFiltFromFile(FolderData, FileName, lowFreq, highFreq, bDivide)

if( ~strcmp(FolderData(end),filesep) )
    FolderData = strcat(FolderData, filesep);
end

fprintf('Opening: %s \n', [FolderData FileName '.dat']);
Infos = matfile([FolderData FileName '.mat']);
fid = fopen([FolderData FileName '.dat']);
OutData = fread(fid, inf, '*single');
OutData = reshape(OutData, Infos.datSize(1,1), Infos.datSize(1,2), []);


% Temporal filtering butterworth
if( lowFreq > 0 )
    UseLPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, lowFreq, Infos.Freq); %Fluo lower Freq
    lpass = design(f,'butter');
else
    UseLPFilt = 0;
end
if( highFreq > 0 )
    UseHPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, highFreq, Infos.Freq);   %Fluo Higher Freq
    hpass = design(f,'butter');
else
    UseHPFilt = 0;
end


dims = size(OutData);
% Hd = zeros(dims,'single');
PrcLims = round(linspace(1, dims(1), 11));
fprintf('Progress: ');
for ind = 1:dims(1)
    Signal = double(squeeze(OutData(ind,:,:)));
    if( UseLPFilt )
        LP_lowCutOff = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
    else
        LP_lowCutOff = ones(size(Signal));
    end
    if( UseHPFilt )
        LP_highCutOff = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
    else
        LP_highCutOff = Signal;
    end
    
    if( bDivide )
        OutData(ind,:,:) = single(LP_highCutOff./LP_lowCutOff);
    else
        OutData(ind,:,:) = single(LP_highCutOff-LP_lowCutOff);
    end
    
    if( any(ind == PrcLims) )
        idx = find(ind == PrcLims);
        fprintf('%d%%..', 10*(idx-1));
    end
    
end
fprintf('\n');    

end

function OutData = NormFiltDirect(FolderData, OutData, lowFreq, highFreq, bDivide, Freq)

if isempty(Freq)
    FileData = dir(fullfile(FolderData, '*.mat'));
    Tags = {'red','green','yellow','fluo_'};
    idx = arrayfun(@(x) contains(FileData(x).name, Tags), 1:size(FileData,1));
    idx = find(idx,1,'first');
    Infos = matfile(fullfile(FolderData, FileData(idx).name));
    Freq = Infos.Freq;
end
% Temporal filtering butterworth
if( lowFreq > 0 )
    UseLPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, lowFreq, Freq); %Fluo lower Freq
    lpass = design(f,'butter');
else
    UseLPFilt = 0;
end
if( highFreq > 0 )
    UseHPFilt = 1;
    f = fdesign.lowpass('N,F3dB', 4, highFreq, Freq);   %Fluo Higher Freq
    hpass = design(f,'butter');
else
    UseHPFilt = 0;
end


dims = size(OutData);
% Hd = zeros(dims,'single');
PrcLims = round(linspace(1, dims(1), 11));
fprintf('Progress: ');
for ind = 1:dims(1)
    Signal = double(squeeze(OutData(ind,:,:)));
    if( UseLPFilt )
        LP_lowCutOff = filtfilt(lpass.sosMatrix, lpass.ScaleValues, Signal')';
    else
        LP_lowCutOff = ones(size(Signal));
    end
    if( UseHPFilt )
        LP_highCutOff = filtfilt(hpass.sosMatrix, hpass.ScaleValues, Signal')';
    else
        LP_highCutOff = Signal;
    end
    
    if( bDivide )
        OutData(ind,:,:) = single(LP_highCutOff./LP_lowCutOff);
    else
        OutData(ind,:,:) = single(LP_highCutOff-LP_lowCutOff);
    end
    
    if( any(ind == PrcLims) )
        idx = find(ind == PrcLims);
        fprintf('%d%%..', 10*(idx-1));
    end
    
end
fprintf('\n');   
end