function out = Ana_IOI_FullFrame(FolderName, verbose)

%%%%%%%%%%% 
% Opening %
%%%%%%%%%%%
%List of files to be included in the analysis
FileList = dir([FolderName filesep 'Data_*.mat']);
Ws = []; Hs = []; Ts = []; Fs = [];
if( isempty(FileList) )
    %No compatible files were found
    disp(['No data files found in ' FolderName ' Folder.']);
    disp('Analysis will not run');
    return;
end

%Green channel detected
IsThereGreen = false;
if( ~isempty(strfind([FileList.name],'green')) )
    IsThereGreen = true;
    Dat_Gptr = matfile([FolderName filesep 'Data_green.mat']);
    nrows = Dat_Gptr.datSize(1,1);
    ncols = Dat_Gptr.datSize(1,2);
    nframes = Dat_Gptr.datLength;
    Ws = ncols;
    Hs = nrows;
    Ts = nframes - 1;
    Fs = Dat_Gptr.Freq;
    gDatPtr = memmapfile(Dat_Gptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes 
end
%Yellow channel detected
IsThereYellow = false;
if( ~isempty(strfind([FileList.name],'yellow')) )
    IsThereYellow = true;
    Dat_Yptr = matfile([FolderName filesep 'Data_yellow.mat']);
    nrows = Dat_Yptr.datSize(1,1);
    ncols = Dat_Yptr.datSize(1,2);
    nframes = Dat_Yptr.datLength;
    Ws = [Ws, ncols];
    Hs = [Hs, nrows];
    Ts = [Ts, nframes-1];
    Fs = [Fs, Dat_Yptr.Freq];
    yDatPtr = memmapfile(Dat_Yptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes 
end
%Red channel detected
IsThereRed = false;
if( ~isempty(strfind([FileList.name],'red')) )
    IsThereRed = true;
    Dat_Rptr = matfile([FolderName filesep 'Data_red.mat']);
    nrows = Dat_Rptr.datSize(1,1);
    ncols = Dat_Rptr.datSize(1,2);
    nframes = Dat_Rptr.datLength;
    Ws = [Ws, ncols];
    Hs = [Hs, nrows];
    Ts = [Ts, nframes-1];
    Fs = [Fs, Dat_Rptr.Freq];
    rDatPtr = memmapfile(Dat_Rptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes
end

%Confirm data dimensions
if( length(unique(Ws)) > 1 )
    disp('Channels have unmatching dimensions');
    disp('Analysis will stop here.');
    return;
end
iWidth = Ws(1);
if( length(unique(Hs)) > 1 )
    disp('Channels have unmatching dimensions');
    disp('Analysis will stop here.');
    return;
end
iHeight = Hs(1);
NbFrames = min(Ts); 
FreqHb = min(Fs);
clear Ws Hs Ts Fs;

%%%%%%%%%%%%%%%%% 
% Output Config %
%%%%%%%%%%%%%%%%%
%Creation of output file
if( exist([FolderName filesep 'Data_Hbs.mat'], 'file') )
    delete([FolderName filesep 'Data_Hbs.mat']);
    delete([FolderName filesep 'HbO.dat']);
    delete([FolderName filesep 'HbR.dat']);
end
OutputFile = matfile([FolderName filesep 'Data_Hbs.mat']);
OutputFile.datFileHbO = [FolderName filesep 'HbO.dat'];
OutputFile.datFileHbR = [FolderName filesep 'HbR.dat'];
OutputFile.datLength = NbFrames;
OutputFile.datSize = [iWidth, iHeight];
OutputFile.Freq = FreqHb;
fHbO = fopen([FolderName filesep 'HbO.dat'], 'w');
fHbR = fopen([FolderName filesep 'HbR.dat'], 'w');

%Filtering Parameters
f = fdesign.lowpass('N,F3dB', 4, 0.5, 5);
hpass = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, 1/60, 5);
lpass = design(f,'butter');
clear f;

clear FileList nframes

%%%%%%%%%%%%%%%%%%%%%%%%%  
% Hb Concentration Calc %
%%%%%%%%%%%%%%%%%%%%%%%%%
%HbO and HbR computation parameters
whichSystem = 0;
whichCurve = 'Dunn';
rescaling_factor = 1e6;
lambda1=450;
lambda2=700;
npoints=1000;
baseline_hbt = 100;
baseline_hbo = 60;
baseline_hbr = 40;

eps_pathlength = ioi_epsilon_pathlength(lambda1,lambda2,npoints,whichSystem,whichCurve,baseline_hbt,baseline_hbo,baseline_hbr);
if( IsThereGreen )
    A = eps_pathlength;
else
    A = [eps_pathlength(1,:); eps_pathlength(3,:)];
end
Ainv=single(rescaling_factor*pinv(A)); % A Inv devrait etre 3x2 ou 2x3 (3=couleurs, 2=hbo,hbr)
clear which* rescaling_factor lambda* npoints baseline_* eps_pathlength A

%For each row, filter each channels and compute HbO, HbR...
HbO = zeros(iHeight, iWidth, NbFrames,'single');  
HbR = zeros(iHeight, iWidth, NbFrames,'single');  
PrcTags = linspace(1, double(iHeight), 11); indPr = 2;
for ind = 1:iHeight
    if( verbose )
        disp(ind);
    end
    %Read, filter and detrend data
    if( IsThereRed )
        pR = rDatPtr.Data(ind:iHeight:(iHeight*iWidth*NbFrames));
        pR = reshape(pR, iWidth, [])';
        Rioi= filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(pR));
        Rbase = filtfilt(lpass.sosMatrix,lpass.ScaleValues, double(pR));
        Rnorm = Rioi./Rbase;
        Rnorm = single(Rnorm);
    end
    if( IsThereYellow )
        pY = yDatPtr.Data(ind:iHeight:(iHeight*iWidth*NbFrames));
        pY = reshape(pY, iWidth, [])';
        Yioi= filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(pY));
        Ybase = filtfilt(lpass.sosMatrix,lpass.ScaleValues, double(pY));
        Ynorm = Yioi./Ybase;        
        Ynorm = single(Ynorm);
    end
    if( IsThereGreen )
        pG = gDatPtr.Data(ind:iHeight:(iHeight*iWidth*NbFrames));
        pG = reshape(pG, iWidth, [])';
        Gioi= filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(pG));
        Gbase = filtfilt(lpass.sosMatrix,lpass.ScaleValues, double(pG));
        Gnorm = Gioi./Gbase;        
        Gnorm = single(Gnorm);    
    end
    
    %Compute HbO & HbR
    if( IsThereGreen )
        Cchan = cat(2, Rnorm(:), Gnorm(:), Ynorm(:));
    else
        Cchan = cat(2, Rnorm(:), Ynorm(:));
    end
    LogCchan = -log10(Cchan);
    Hbs = Ainv*LogCchan';

    %Save
    HbO(ind,:,:) = reshape(Hbs(1,:), [], iWidth)';
    HbR(ind,:,:) = reshape(Hbs(2,:), [], iWidth)';
         
    if( ind >= PrcTags(indPr) )
        fprintf('%d%%...', 10*(indPr-1));
        indPr = indPr + 1;
    end
end
fprintf('\n');
disp('Saving Hb values');
fwrite(fHbO, HbO, 'single');
fwrite(fHbR, HbR, 'single');
OutputFile.Stim = Dat_Rptr.Stim;
fclose(fHbO);
fclose(fHbR);
end
