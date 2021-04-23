function varargout = HemoCorrection(Folder, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Infos:
%
% This function is used to remove the hemodynamic fluctuations from any
% fluorescence signal. 
%
% Inputs:
% 1. Folder: Folder contaning the dataset to work with.
% 2. Varargin -> if empty: a dialog box will be prompt to ask user which
%                   channels to use to do the correction.
%             -> cell array of string: to specify which channels to use.
%                   Ex: HemoCorrection(pwd, {'Red', 'Green'});
% Ouput:
% - If an output is set, the result of the correction will be given back
% through this output. All the data in the folder will remain unchanged.
% - If no output is specified, the fChan.dat files in Folder will be overwritten
% with the corrected data.
% 
% Exemples:
%
% 1- HemoCorrection(pwd); 
% The fluorescence dat files in the folder will be overwriten and a dialog
% box will be used in order to select which channels must be used to
% compute the correction.
% 2- NewFluo = HemoCorrection(pwd, {'Green'});
% The dat files in the folder won't be overwriten. The new corrected
% fluorescence data will be in NewFluo. Only the green channel will be used
% to compute the correction.
% 3- HemoCorrection(pwd, {'Red, 'Green', 'Amber'});
% fChan_475.dat will be overwriten with the corrected fluorescence data.
% All three hemodynamic channels will be used to compute the correction.

if( ~strcmp(Folder(end),filesep) )
    Folder = strcat(Folder, filesep);
end

if( nargin <= 1 )
    cList = dir([Folder '*Chan.dat']);
    fn = {};
    for ind = 1:size(cList,1)
        if( ~strcmp(cList(ind).name(1),'f') )
            fn{end+1} = cList(ind).name;
        end
    end
    
    [idx, tf] = listdlg('PromptString',{'Select channels to be used to',...
    'compute hemodynamic correction.',''},...
    'ListString',fn);
    
    if( tf == 0 )
        return;
    end
    
    fn = fn(idx);
    clear cList idx ind tf;
else
    fn = {};
    tmp = varargin{1};
    for ind = 1:size(tmp,2)
        tag = lower(tmp{ind});
        switch tag
            case 'red'
                fn{end+1} = 'rChan.dat';
            case {'amber', 'yellow'}
                fn{end+1} = 'yChan.dat';
            case 'green'
                fn{end+1} = 'gChan.dat';
        end
    end
end

%Load data:
DataH = [];
Infos = matfile('Data_Fluo.mat');
for ind = 1:size(fn,2)
   eval(['fid = fopen(''' fn{ind} ''');']);
   tmp = fread(fid, inf, 'single=>single');
   fclose(fid);
   tmp = reshape(tmp, Infos.datSize(1,1), Infos.datSize(1,2), []);
   tmp = imgaussfilt(tmp,1, 'Padding', 'symmetric');
   tmp = reshape(tmp, [], size(tmp,3));
   tmp = tmp./mean(tmp,2) - 1;
   DataH = cat(3, DataH, tmp);
end
DataH = movmean(DataH,Infos.Freq,2); 
clear fid fn ind tmp;
fList = dir([Folder 'fChan*.dat']);

%Correction:
for ind = 1:size(fList,1)
    eval(['fid = fopen(''' fList(ind).name ''');']);
    fData = fread(fid, inf, 'single=>single');
    fclose(fid);
    fData = reshape(fData,Infos.datSize(1,1)*Infos.datSize(1,2), []);
    tmp = fData./mean(fData,2);
    
    A = zeros(size(fData));
    %Bf = zeros(size(fData,1),(2 + size(DataH,3)));
    warning('off', 'MATLAB:rankDeficientMatrix');
    for indF = 1:size(fData,1)
        if( length(size(DataH)) == 2 )
            X = [ones(1, size(fData,2)); linspace(0,1,size(fData,2)); squeeze(DataH(indF,:))];
        else
            X = [ones(1, size(fData,2)); linspace(0,1,size(fData,2)); squeeze(DataH(indF,:,:))'];
        end
        B = X'\tmp(indF,:)';
        A(indF,:) = (X'*B)';
        %Bf(indF,:) = B;
    end
    warning('on', 'MATLAB:rankDeficientMatrix');
    clear B X tmp;
    
    fData = fData./A;
    
    fData = reshape(fData, Infos.datSize(1,1), Infos.datSize(1,2), []);
    if( nargout == 0 )
        eval(['fid = fopen(''' fList(ind).name ''');']);
        fwrite(fid, fData, 'single');
        fclose(fid);
    else
        varargout{ind} = fData;
    end 
end

end