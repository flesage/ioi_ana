function Ana_Redirection(FolderName)

FileList = dir([FolderName filesep 'Data_speckle.mat']);
if( ~isempty(FileList) )
    disp(['Speckle data files found in ' FolderName ' Folder.']);
    disp('Speckle Analysis starting...');
    Ana_Speckle(FolderName);
else
    disp(['Fluorescence data files found in ' FolderName ' Folder.']);
    disp('Fluorescence Analysis starting...');
    Ana_Fluo(FolderName);
end

end