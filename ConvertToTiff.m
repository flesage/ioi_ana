function ConvertToTiff(FolderPath)

if( ~strcmp(FolderPath(end),filesep) )
    FolderPath = strcat(FolderPath, filesep);
end
disp(['Starting Tiff conversion for: ' FolderPath]);
ChanList = dir([FolderPath '*Chan.dat']);
InfList = dir([FolderPath 'Data_*.mat']);
for indC = 1:size(ChanList,1)
    ChanName = ChanList(indC).name;
    fid = fopen([FolderPath ChanName]);
    dat = fread(fid,inf,'single');
    idx = arrayfun(@(X) contains(InfList(X).name,['_' ChanName(1)]), 1:size(InfList,1)); 
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
    end
    
    dat = uint16(dat);
    imwrite(dat(:, :, 1), outFName, 'Compression','none');        
    for indF = 2:size(dat,3)
        imwrite(dat(:, :, indF), outFName, 'WriteMode', 'append',  'Compression','none');        
    end
end
disp('Tiff conversion is done.')
end