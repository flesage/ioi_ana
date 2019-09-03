function out = Ana_Fluo(FolderName, b_HbCorr, b_HilbertF)
%%%%%%%%%%%%%%%%%%%%%%
% Validation & opening
%%%%%%%%%%%%%%%%%%%%%%

if( ~strcmp(FolderName(end), filesep) )
    FolderName = strcat(FolderName, filesep);
end

fprintf('Opening files.\n');

FileList = dir([FolderName 'Data_Fluo.mat']);
if( isempty(FileList) )
    disp(['No Fluorescence data files found in ' FolderName ' Folder.']);
    disp('Fluorescence Analysis will not run');
    return;
end

FInfo =  matfile([FolderName 'Data_Fluo.mat'], 'Writable', true);
fid = fopen([FolderName 'fChan.dat']);
Fluo = fread(fid, inf, 'single=>single');
fclose(fid);

if( exist([FolderName 'Data_Hbs.mat'], 'file') && b_HbCorr )
    FluoCorrGen(FolderName);
    
    HBInfo = matfile([FolderName 'Data_Hbs.mat']);
    fid = fopen([FolderName 'hCorr.dat']);
    Corr = fread(fid, inf, 'single');
    Corr = single(Corr);
    fclose(fid);

    fprintf('Applying HB Correction.\n');
    
    Fluo = Fluo(1:length(Corr));
    Fluo = Fluo.*Corr;
end

Fluo = reshape(Fluo, FInfo.datSize(1,1), FInfo.datSize(1,2),[]);
FInfo.datLength = size(Fluo,3);

fprintf('Filtering.\n');
if( ~b_HilbertF )
    dims = size(Fluo);
    f = fdesign.lowpass('N,F3dB', 4, FInfo.Freq/2, FInfo.Freq);
    hpass = design(f,'butter');
    f = fdesign.lowpass('N,F3dB', 4, 1/10, FInfo.Freq);
    lpass = design(f,'butter');
    
    dFh = single(filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(reshape(permute(Fluo, [3 1 2]), dims(3),[]))));
    dFh = reshape(dFh', FInfo.datSize(1,1), FInfo.datSize(1,2), []);
    dFl = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(reshape(permute(Fluo, [3 1 2]), dims(3),[]))));
    dFl = reshape(dFl', dims);
    Fluo = (dFh - dFl)./dFl;
    clear dFh dFl f;
else
    dims = size(Fluo);
    [~, ylower] = envelope(reshape(Fluo,[],dims(3))',4*FInfo.Freq,'peak');
    ylower = reshape(ylower', dims);
    Fluo = (Fluo - ylower)./ylower;
    for ind = 1:dims(3)
        Fluo(:,:,ind) = medfilt2(squeeze(Fluo(:,:,ind)),[5 5],'symmetric');
    end
end

fprintf('Saving.\n');
fid = fopen([FolderName 'fChan.dat'],'w');
fwrite(fid,Fluo,'single');
fclose(fid);

end

