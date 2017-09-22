function out = BandingNoiseFilter(FolderName)

%%%%%%%%%%% 
% Opening %
%%%%%%%%%%%
%List of files to be included in the analysis
FileList = dir([FolderName filesep 'Data_*.mat']);
if( isempty(FileList) )
    %No compatible files were found
    disp(['No data files found in ' FolderName ' Folder.']);
    disp('Analysis will not run');
    return;
end

AcqInfoStream = readtable([FolderName filesep 'info.txt'],...
    'Delimiter',':','ReadVariableNames',false, 'ReadRowNames',true);

%Green channel detected
if( ~isempty(strfind([FileList.name],'green')) )
    fprintf('Green channel banding noise filtering...')
    Dat_ptr = matfile([FolderName filesep 'Data_green.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    Fs = Dat_ptr.Freq;
        
    DatPtr = memmapfile(Dat_ptr.datFile,...
        'Format', 'single');
    
    fGreen = BN_Filter(DatPtr, nrows, ncols, Fs);
    fprintf('done\n')
    clear DatPtr;
    fprintf('Saving...\n')
    fid = fopen(Dat_ptr.datFile);
    fwrite(fid, fGreen, 'single');
    fclose(fid);
    
    clear Dat_ptr nrows ncols Fs fid fGreen;
end
%Yellow channel detected
if( ~isempty(strfind([FileList.name],'yellow')) )
    fprintf('Yellow channel banding noise filtering...')
    Dat_ptr = matfile([FolderName filesep 'Data_yellow.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    Fs = Dat_ptr.Freq;
    DatPtr = memmapfile(Dat_ptr.datFile,...
        'Format', 'single', 'Writable', true);
    fYellow = BN_Filter(DatPtr, nrows, ncols, Fs);
    fprintf('done\n')
    fprintf('Saving...\n')
    DatPtr.Data(:) = fYellow(:);        
    clear Dat_ptr nrows ncols Fs DatPtr fYellow;
end
%Red channel detected
if( ~isempty(strfind([FileList.name],'red')) )
    fprintf('Red channel banding noise filtering...')
    Dat_ptr = matfile([FolderName filesep 'Data_red.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    Fs = Dat_ptr.Freq;
    DatPtr = memmapfile(Dat_ptr.datFile,...
        'Format', 'single', 'Writable', true);
    fRed = BN_Filter(DatPtr, nrows, ncols, Fs);
    fprintf('done\n')
    fprintf('Saving...\n')
    DatPtr.Data(:) = fRed(:);    
        
    clear Dat_ptr nrows ncols Fs DatPtr fRed;
end
%Speckle channel detected
if( ~isempty(strfind([FileList.name],'speckle')) )
    fprintf('Speckle channel banding noise filtering...')
    Dat_ptr = matfile([FolderName filesep 'Data_speckle.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    Fs = Dat_ptr.Freq;
    DatPtr = memmapfile(Dat_ptr.datFile,...
        'Format', 'single','Writable', true);
    fSpeckle = BN_Filter(DatPtr, nrows, ncols, Fs);
    fprintf('done\n')
    fprintf('Saving...\n')
    DatPtr.Data(:) = fSpeckle(:);    
        
    clear Dat_ptr nrows ncols Fs DatPtr fSpeckle;
end
%Fluo channel detected
if( ~isempty(strfind([FileList.name],'Fluo')) )
    fprintf('Fluo channel banding noise filtering...')
    Dat_ptr = matfile([FolderName filesep 'Data_Fluo.mat']);
    nrows = Dat_ptr.datSize(1,1);
    ncols = Dat_ptr.datSize(1,2);
    Fs = Dat_ptr.Freq;
    DatPtr = memmapfile(Dat_ptr.datFile,...
        'Format', 'single','Writable', true);
    fFluo = BN_Filter(DatPtr, nrows, ncols, Fs);
    fprintf('done\n')
    fprintf('Saving...\n')
    DatPtr.Data(:) = fFluo(:);   
        
    clear Dat_ptr nrows ncols Fs DatPtr fFluo;
end


function DeNoisedDat = BN_Filter(DatPtr, nrows, ncols, Fs)
f = fdesign.highpass('N,F3dB', 4, 0.1, Fs);
hpass = design(f,'butter');
clear f;
    
Frame = DatPtr.Data(:);
Frame = reshape(Frame, nrows, ncols, []);
msk = imdilate(max(Frame,[],3) >= 4095, strel('disk',5));

Frame = reshape(Frame,[],size(Frame,3));
Frame(msk(:),:) = mean(reshape(Frame(~msk(:),:),[],1));
Frame = reshape(Frame,nrows,ncols,[]);

Frame = permute(Frame,[3 1 2]);
Lim = 1024^3; 
Size = length(Frame(:))*4;
NIter = Size/Lim;
Iter = round(linspace(1, double(nrows*ncols+1), ceil(NIter)));
hf_Frame = zeros(size(Frame),'single');
for ind = 2:length(Iter)
    Tmp = Frame(:,Iter(ind-1):(Iter(ind)-1));
    hf_Frame(:,Iter(ind-1):(Iter(ind)-1)) = filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(Tmp));
end
clear Tmp;
Frame = permute(Frame,[2 3 1]);
hf_Frame = permute(hf_Frame,[2 3 1]);

VertData = Frame;
minD = min(size(Frame,1),size(Frame,2));
for indF = 1:size(Frame,3)
    Tmp = squeeze(hf_Frame(:,:,indF));
    VertData(:,:,indF) = Tmp - single(xRemoveStripesVertical(Tmp, nextpow2(minD)-4, 'db4', 2));
end

DeNoisedDat = Frame - VertData;

end
end