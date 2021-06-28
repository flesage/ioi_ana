function ConvertToTiff(FolderPath)

if( ~strcmp(FolderPath(end),filesep) )
    FolderPath = strcat(FolderPath, filesep);
end
disp(['Starting Tiff conversion for: ' FolderPath]);
ChanList = dir([FolderPath '*.dat']);
if( isempty(ChanList) )
    disp(['No dat file found. End of convertion.']);
    return;
end

if( exist([FolderPath 'red.mat'],'file') )
    InfList = dir([FolderPath 'red.mat']);
end
if( exist([FolderPath 'green.mat'],'file') )
    if( ~exist('InfList') )
        InfList = dir([FolderPath 'green.mat']);
    end
    InfList(end+1) = dir([FolderPath 'green.mat']);
end
if( exist([FolderPath 'yellow.mat'],'file') )
    if( ~exist('InfList') )
        InfList = dir([FolderPath 'yellow.mat']);
    end
    InfList(end+1) = dir([FolderPath 'yellow.mat']);
end
if( exist([FolderPath 'fluo*.mat'],'file') )
    fluoList = dir([FolderPath 'fluo*.mat']);
    for ind = 1:size(fluoList,2)
        if( ~exist('InfList') )
            InfList = dir([FolderPath fluoList(ind).name]);
        end
        InfList(end+1) = dir([FolderPath fluoList(ind).name]);
    end
end
if( isempty(InfList) )
    InfList = dir([FolderPath 'Data_*.mat']);
end
if( isempty( InfList ) )
    disp(['No mat file found. End of convertion.']);
    return;
end

for indC = 1:size(ChanList,1)
    ChanName = ChanList(indC).name;
    fid = fopen([FolderPath ChanName]);
    dat = fread(fid,inf,'single');
    tag = ChanName(1:(strfind(ChanName,'.') - 1));
    idx = arrayfun(@(X) contains(InfList(X).name, tag,'IgnoreCase',true), 1:size(InfList,2)); 
    Infos = matfile([FolderPath InfList(find(idx)).name]);
    dat = reshape(dat, Infos.datSize(1,1), Infos.datSize(1,2), []);
    
    outFName = [FolderPath 'img_'];
    switch(ChanName(1))
        case 'g'
            disp('Saving green channel.');
            outFName = strcat(outFName, 'Green.tif');
        case 'r'
            disp('Saving red channel.');
            outFName = strcat(outFName, 'Red.tif');
        case 'y'
            disp('Saving amber channel.');
            outFName = strcat(outFName, 'Yellow.tif');
        case 'f'
            disp('Saving fluo channel.');
            outFName = strcat(outFName, 'Fluo.tif');
        case 's'
            disp('Saving speckle channel.');
            outFName = strcat(outFName, 'Speckle.tif');
    end
    
    dat = uint16(dat);
    imwrite(dat(:, :, 1), outFName, 'Compression','none');        
    for indF = 2:size(dat,3)
        imwrite(dat(:, :, indF), outFName, 'WriteMode', 'append',  'Compression','none');        
    end
end
disp('Tiff conversion is done.')
end