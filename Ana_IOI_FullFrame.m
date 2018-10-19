function out = Ana_IOI_FullFrame(FolderName, verbose, OStream)

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
if( contains([FileList.name],'green') )
    IsThereGreen = true;
    Dat_Gptr = matfile([FolderName filesep 'Data_green.mat'],...
        'Writable', true);
    nrows = Dat_Gptr.datSize(1,1);
    ncols = Dat_Gptr.datSize(1,2);
    nframes = Dat_Gptr.datLength;
    Ws = ncols;
    Hs = nrows;
    Ts = nframes;
    Fs = Dat_Gptr.Freq;
   %Flip Data (Time x Xaxis x Yaxis)
    if( ~strcmp( Dat_Gptr.FirstDim, 't') )
        fid = fopen(Dat_Gptr.datFile, 'r');
        dat = zeros(nframes, nrows, ncols, 'single');
        for indT = 1:nframes
            dat(indT,:,:) = reshape(fread(fid, nrows*ncols, 'single'),nrows,[]);
        end
        fclose(fid);
        fid = fopen(Dat_Gptr.datFile, 'w');
        Dat_Gptr.FirstDim = 't';
        fwrite(fid, dat, 'single');
        fclose(fid);
        clear dat
    end
    gDatPtr = memmapfile(Dat_Gptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes 
end
%Yellow channel detected
IsThereYellow = false;
if( contains([FileList.name],'yellow') )
    IsThereYellow = true;
    Dat_Yptr = matfile([FolderName filesep 'Data_yellow.mat'],...
        'Writable', true);
    nrows = double(Dat_Yptr.datSize(1,1));
    ncols = double(Dat_Yptr.datSize(1,2));
    nframes = double(Dat_Yptr.datLength);
    Ws = [Ws, ncols];
    Hs = [Hs, nrows];
    Ts = [Ts, nframes];
    Fs = [Fs, Dat_Yptr.Freq];
    %Flip Data (Time x Xaxis x Yaxis)
    if( ~strcmp( Dat_Yptr.FirstDim, 't') )
        fid = fopen(Dat_Yptr.datFile, 'r');
        dat = zeros(nframes, nrows, ncols, 'single');
        for indT = 1:nframes
            dat(indT,:,:) = reshape(fread(fid, nrows*ncols, 'single'),nrows,[]);
        end
        fclose(fid);
        fid = fopen(Dat_Yptr.datFile, 'w');
        Dat_Yptr.FirstDim = 't';
        fwrite(fid, dat, 'single');
        fclose(fid);
        clear dat
    end
    yDatPtr = memmapfile(Dat_Yptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes 
end
%Red channel detected
IsThereRed = false;
if( contains([FileList.name],'red') )
    IsThereRed = true;
    Dat_Rptr = matfile([FolderName filesep 'Data_red.mat'],...
        'Writable', true);
    nrows = Dat_Rptr.datSize(1,1);
    ncols = Dat_Rptr.datSize(1,2);
    nframes = Dat_Rptr.datLength;
    Ws = [Ws, ncols];
    Hs = [Hs, nrows];
    Ts = [Ts, nframes];
    Fs = [Fs, Dat_Rptr.Freq];
    if( ~strcmp( Dat_Rptr.FirstDim, 't') )
        fid = fopen(Dat_Rptr.datFile, 'r');
        dat = zeros(nframes, nrows, ncols, 'single');
        for indT = 1:nframes
            dat(indT,:,:) = reshape(fread(fid, nrows*ncols, 'single'),nrows,[]);
        end
        fclose(fid);
        fid = fopen(Dat_Rptr.datFile, 'w');
        Dat_Rptr.FirstDim = 't';
        fwrite(fid, dat, 'single');
        fclose(fid);
        clear dat
    end
    rDatPtr = memmapfile(Dat_Rptr.datFile,...
        'Format', 'single');
    clear nrows ncols cframes
end

%Is all required colors available for HB calculation?
if( IsThereRed + IsThereYellow + IsThereGreen < 2 )
    disp('*** Impossible to compute Hb concentrations. More color channels needed.');
    fprintf('\n');
    return;
end

%Confirm data dimensions
if( length(unique(Ws)) > 1 )
    disp('Channels have unmatching dimensions');
    disp('Analysis will stop here.');
    return;
end
iWidth = double(Ws(1));
if( length(unique(Hs)) > 1 )
    disp('Channels have unmatching dimensions');
    disp('Analysis will stop here.');
    return;
end
iHeight = double(Hs(1));
NbFrames = double(min(Ts)); 
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
f = fdesign.lowpass('N,F3dB', 4, 0.5, FreqHb);
hpass = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, 1/60, FreqHb);
lpass = design(f,'butter');
clear f FileList nframes

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
if( IsThereGreen && IsThereYellow && IsThereRed )
    A = eps_pathlength;
elseif( IsThereYellow && IsThereRed )
    A = [eps_pathlength(1,:); eps_pathlength(3,:)];
elseif( IsThereGreen && IsThereRed )
    A = [eps_pathlength(1,:); eps_pathlength(2,:)];
elseif( IsThereGreen && IsThereYellow )
    A = [eps_pathlength(2,:); eps_pathlength(3,:)];
end
Ainv=single(rescaling_factor*pinv(A)); % A Inv devrait etre 3x2 ou 2x3 (3=couleurs, 2=hbo,hbr)
clear which* rescaling_factor lambda* npoints baseline_* eps_pathlength A

%For each row, filter each channels and compute HbO, HbR...
PrcTags = linspace(1, double(iWidth), 11); indPr = 2;
if( ~isempty(OStream) )
    StaticStr = OStream.String;
end
for ind = 1:iWidth
    %Read, filter and detrend data
    if( IsThereRed )
        rDatPtr = memmapfile(Dat_Rptr.datFile,...
            'Offset', 4*double(ind - 1)*iHeight*Dat_Rptr.datLength,...
            'Format', 'single');
        pR = rDatPtr.Data(1:(iHeight*Dat_Rptr.datLength));
        pR = reshape(pR, [], iHeight);
        pR = pR(1:NbFrames,:);
        Rioi= single(filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(pR)));
        Rbase = single(filtfilt(lpass.sosMatrix,lpass.ScaleValues, double(pR)));
        Rnorm = single(Rioi./Rbase);
        clear Rioi Rbase pR;
    end
    if( IsThereYellow )
        yDatPtr = memmapfile(Dat_Yptr.datFile,...
            'Offset', 4*double(ind - 1)*iHeight*Dat_Yptr.datLength,...
            'Format', 'single');
        pY = yDatPtr.Data(1:(iHeight*Dat_Yptr.datLength));
        pY = reshape(pY, [], iHeight);
        pY = pY(1:NbFrames,:);
        Yioi= single(filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(pY)));
        Ybase = single(filtfilt(lpass.sosMatrix,lpass.ScaleValues, double(pY)));
        Ynorm = single(Yioi./Ybase);      
        clear Yioi Ybase pY;
    end
    if( IsThereGreen )
        gDatPtr = memmapfile(Dat_Gptr.datFile,...
            'Offset', 4*double(ind - 1)*iHeight*Dat_Gptr.datLength,...
            'Format', 'single');
        pG = gDatPtr.Data(1:(iHeight*Dat_Gptr.datLength));
        pG = reshape(pG, [], iHeight);
        pG = pG(1:NbFrames,:);
        Gioi= single(filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(pG)));
        Gbase = single(filtfilt(lpass.sosMatrix,lpass.ScaleValues, double(pG)));
        Gnorm = single(Gioi./Gbase);    
        clear Gioi Gbase pG;
    end
    
    %Compute HbO & HbR
    if( IsThereGreen && IsThereYellow && IsThereRed )
        Cchan = cat(2, Rnorm(:), Gnorm(:), Ynorm(:));
    elseif( IsThereYellow && IsThereRed )
        Cchan = cat(2, Rnorm(:), Ynorm(:));
    elseif( IsThereGreen && IsThereRed )
        Cchan = cat(2, Rnorm(:), Gnorm(:));
    elseif( IsThereGreen && IsThereYellow )
        Cchan = cat(2, Gnorm(:), Ynorm(:));
    end
  
    LogCchan = -log10(Cchan);
    Hbs = Ainv*LogCchan';

    %Save
    if( ind > 1)
        fHbO = fopen([FolderName filesep 'HbO.dat'], 'a');
        fHbR = fopen([FolderName filesep 'HbR.dat'], 'a');
    end
    fwrite(fHbO,reshape(Hbs(1,:), [], iHeight), 'single');
    %HbO(ind,:,:) = reshape(Hbs(1,:), [], iWidth)';
    fwrite(fHbR,reshape(Hbs(2,:), [], iHeight), 'single');
    %HbR(ind,:,:) = reshape(Hbs(2,:), [], iWidth)';
    fclose(fHbO);
    fclose(fHbR);
    
    if( ind >= PrcTags(indPr) )
        if( isempty(OStream) )
            fprintf('%d%% .. ', 10*(indPr-1));            
        else
            OStream.String = sprintf('%s\r%s',...
                ['Completion: ' int2str(10*(indPr-1)) '%'],...
                StaticStr);
            drawnow;
        end
        indPr = indPr + 1;
    end
end
if( isempty(OStream) )
    fprintf('\n');
else
    OStream.String = StaticStr;
    OStream.String = sprintf('%s\r%s',...
        'Done.',...
        OStream.String);
    drawnow;
end
OutputFile.Stim = Dat_Rptr.Stim;

fHbO = fopen([FolderName filesep 'HbO.dat'], 'r+');
dat = fread(fHbO,inf,'single');
frewind(fHbO);
dat = reshape(dat, NbFrames, iWidth, iHeight);
dat = permute(dat, [1 3 2]);
fwrite(fHbO, dat,'single'); 
fclose(fHbO);
fHbR = fopen([FolderName filesep 'HbR.dat'], 'r+');
dat = fread(fHbR,inf,'single');
frewind(fHbR);
dat = reshape(dat, NbFrames, iWidth, iHeight);
dat = permute(dat, [1 3 2]);
fwrite(fHbR, dat,'single'); 
fclose(fHbR);

%Re-Flip Data:
if( IsThereGreen && strcmp(Dat_Gptr.FirstDim, 't') )
    clear gDatPtr;
    fid = fopen(Dat_Gptr.datFile, 'r');
    dat = fread(fid, inf, 'single');
    dat = reshape(dat, [], iWidth, iHeight);
    dat = permute(dat, [2 3 1]);
    fclose(fid);
    fid = fopen(Dat_Gptr.datFile, 'w+');
    Dat_Gptr.FirstDim = 'y';
    fwrite(fid, dat, 'single');
    fclose(fid);
    clear dat
end
if( IsThereYellow && strcmp(Dat_Yptr.FirstDim, 't') )
    clear yDatPtr;
    fid = fopen(Dat_Yptr.datFile, 'r');
    dat = fread(fid, inf, 'single');
    dat = reshape(dat, [], iWidth, iHeight);
    dat = permute(dat, [2 3 1]);
    fclose(fid);
    fid = fopen(Dat_Yptr.datFile, 'w');
    Dat_Yptr.FirstDim = 'y';
    fwrite(fid, dat, 'single');
    fclose(fid);
    clear dat
end
if( IsThereRed && strcmp(Dat_Rptr.FirstDim, 't') )
    clear rDatPtr;
    fid = fopen(Dat_Rptr.datFile, 'r');
    dat = fread(fid, inf, 'single');
    dat = reshape(dat, [], iWidth, iHeight);
    dat = permute(dat, [2 3 1]);
    fclose(fid);
    fid = fopen(Dat_Rptr.datFile, 'w');
    Dat_Rptr.FirstDim = 'y';
    fwrite(fid, dat, 'single');
    fclose(fid);
    clear dat
end
end
