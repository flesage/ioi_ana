function out = Ana_Fluo(FolderName, b_HbCorr)
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
Fluo = fread(fid, inf, 'single');
Fluo = single(Fluo);
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
fbase = ceil(60*FInfo.Freq);
Fbase = medfilt1(Fluo, fbase, [], 3, 'truncate');
Fluo = (Fluo - Fbase)./Fbase;

fprintf('Saving.\n');
fid = fopen([FolderName 'fChan.dat'],'w');
fwrite(fid,Fluo,'single');
fclose(fid);

end

