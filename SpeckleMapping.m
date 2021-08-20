function DatOut = SpeckleMapping(folderPath, sType)

if( ~strcmp(folderPath(end), filesep) )
    folderPath = strcat(folderPath, filesep);
end

disp(['Starting .mat conversion for: ' folderPath]);

if(~exist([folderPath 'speckle.dat'],'file') )
    disp('Speckle.dat is missing. Have you run ImagesClassificiation?');
    return;
end

disp('Opening speckle.dat');

Infos = matfile([folderPath 'speckle.mat']);
fid = fopen([folderPath 'speckle.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat, Infos.datSize(1,1), Infos.datSize(1,2),[]);
dat = -log10(dat./mean(dat,3));

switch sType
    case 'Spatial'
        Kernel = zeros(5,5,1);
        Kernel(:,:,1) = fspecial('disk',2)>0;
        DatOut = stdfilt(dat,Kernel);
    case 'Temporal'
        Kernel = zeros(1,1,5);
        Kernel(3,3,:) = 1;
        DatOut = stdfilt(dat,Kernel);
end

end



