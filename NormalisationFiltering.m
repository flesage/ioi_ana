function OutData = NormalisationFiltering(varargin)
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
%   1- Data:        Data, as a 3D matrix (Y, X, Time)
%   2- Freq:        Sample rate of Data
%   3- lowFreq:     low frequency cut-off, set to 0 to ignore
%   4- highFreq:    high frequency cut-off, set to 0 to ignore 
%   5- bDivide:     if 1, the data returned (below highFreq) is normalised
%                   by the low freq signal (below lowFreq)
%                   if 0, the low freq signal (below lowFreq) is
%                   substracted from the data returned (below highFreq)
%
%   Ex: Dat = NormalisationFiltering(dat, 10, 0, 1, 1);
%   
%   This call would return the data contained in the variable dat with a
%   low-pass at 1 Hz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( isa(varargin{2}, 'char') )
    FolderData = varargin{1};
    FileData = varargin{2};
    lowFreq = varargin{3};
    highFreq = varargin{4};
    bDivide = varargin{5};
    
    OutData = NormFiltFromFile(FolderData, FileData, lowFreq, highFreq, bDivide);
else
    FolderData = varargin{1};
    Data = varargin{2};
    lowFreq = varargin{3};
    highFreq = varargin{4};
    bDivide = varargin{5};
    
    OutData = NormFiltDirect(FolderData, Data, lowFreq, highFreq, bDivide);    
end

end


function OutData = NormFiltFromFile(FolderData, FileData, lowFreq, highFreq, bDivide)

if( ~strcmp(FolderData(end),filesep) )
    FolderData = strcat(FolderData, filesep);
end

fprintf('Opening: %s \n', [FolderData FileData '.dat']);
Infos = matfile([FolderData FileData '.mat']);
fid = fopen([FolderData FileData '.dat']);
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

function OutData = NormFiltDirect(OutData, Freq, lowFreq, highFreq, bDivide)

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