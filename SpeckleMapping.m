function DatOut = SpeckleMapping(folderPath, sType, channel)
%%%%%%%%%%%%%%%%%%%% Speckle Mapping function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the standard deviation (spatialy or temporaly) of speckle
% acquisition. This measure is proportional to the strength of blood flow
% in vessels.
%
% INPUTS:
%
% 1- folderPath: Folder containing the speckle data (called speckle.dat)
%
% 2- sType: how the stdev should be computed. Two options:
%       - Spatial: stdev will be computed on a 5x5 area in the XY plane
%       - Temporal: stdev will be computed on a 5x1 vector in the time dimension 
%
% 3- channel (optional): Channel to analyse, for example 'green', 'red',
% etc. (speckle by default)
%
% OUTPUT:
%
% 1- DatOut: StDev variation over time.

if(nargin < 3)
    channel = 'speckle';
end

if( ~strcmp(folderPath(end), filesep) )
    folderPath = strcat(folderPath, filesep);
end

if(~exist([folderPath channel '.dat'],'file') )
    disp([channel '.dat file is missing. Did you run ImagesClassificiation?']);
    return;
end

disp(['Opening ' channel '.dat']);

Infos = matfile([folderPath channel '.mat']);
fid = fopen([folderPath channel '.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat, Infos.datSize(1,1), Infos.datSize(1,2),[]);
dat = dat./mean(dat,3);

disp('Mapping Computation');
switch sType
    case 'Spatial'
        Kernel = zeros(5,5,1);
        Kernel(:,:,1) = fspecial('disk',2)>0;
        DatOut = stdfilt(dat,Kernel);
    case 'Temporal'
        Kernel = ones(1,1,5);
        
        DatOut = stdfilt(dat,Kernel);
     
end

DatOut = -log10(DatOut);

%Remove outliers
sorted = sort(DatOut(:));
outlier = sorted(round(0.99*size(sorted(:),1)));
DatOut(DatOut>outlier) = outlier;

disp('Done');

%Generate output
% copyfile([folderPath channel '.mat'], [folderPath flow '.mat'])
mFileOut = matfile([folderPath 'flow.mat'], 'Writable', true);
mFileOut.FirstDim = Infos.FirstDim;
mFileOut.Freq = Infos.Freq;
mFileOut.Stim = Infos.Stim;
mFileOut.datLength = Infos.datLength;
mFileOut.datSize = Infos.datSize;
mFileOut.datFile = 'flow.dat';

fid = fopen([folderPath 'flow.dat'],'w'); 
fwrite(fid, single(DatOut), 'single');
fclose(fid);

end
