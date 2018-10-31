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
fHbO = fopen([FolderName filesep 'HbO.dat'], 'W');
fHbR = fopen([FolderName filesep 'HbR.dat'], 'W');

%Filtering Parameters
fbase = ceil(60*FreqHb);
fioi = ceil(2*FreqHb);
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
PrcTags = linspace(1, double(iHeight), 11); indPr = 2;
if( ~isempty(OStream) )
    StaticStr = OStream.String;
end
for ind = 1:iHeight
    %Read, filter and detrend data
    if( IsThereRed )
        pR = rDatPtr.Data(double(ind):iHeight:(iHeight*iWidth*NbFrames));
        pR = reshape(pR, iWidth, [])';
        Rioi= medfilt1(pR,fioi,[],1,'truncate');
        Rbase = medfilt1(pR,fbase,[],1,'truncate');
        Rnorm = Rioi./Rbase;
        clear Rioi Rbase pR;
    end
    if( IsThereYellow )
        pY = yDatPtr.Data(double(ind):iHeight:(iHeight*iWidth*NbFrames));
        pY = reshape(pY, iWidth, [])';
        Yioi = medfilt1(pY,fioi,[],1,'truncate');
        Ybase = medfilt1(pY,fbase,[],1,'truncate');
        Ynorm = Yioi./Ybase;
        clear Yioi Ybase pY;
    end
    if( IsThereGreen )
        pG = gDatPtr.Data(ind:iHeight:(iHeight*iWidth*NbFrames));
        pG = reshape(pG, iWidth, [])';
        Gioi = medfilt1(pG,fioi,[],1,'truncate');
        Gbase = medfilt1(pG,fbase,[],1,'truncate');
        Gnorm = Gioi./Gbase;
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
    fwrite(fHbO,reshape(Hbs(1,:), [], iWidth), 'single');
    %HbO(ind,:,:) = reshape(Hbs(1,:), [], iWidth)';
    fwrite(fHbR,reshape(Hbs(2,:), [], iWidth), 'single');
    %HbR(ind,:,:) = reshape(Hbs(2,:), [], iWidth)';
    
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
if( IsThereRed )
   OutputFile.Stim = Dat_Rptr.Stim;
elseif( IsThereGreen )
   OutputFile.Stim = Dat_Gptr.Stim;
else
   OutputFile.Stim = Dat_Yptr.Stim;
end

fclose(fHbO);
fclose(fHbR);
fHbO = fopen([FolderName filesep 'HbO.dat'], 'r+');
dat = fread(fHbO,inf,'single');
frewind(fHbO);
dat = reshape(dat, NbFrames, iWidth, iHeight);
dat = permute(dat, [3 2 1]);
fwrite(fHbO, dat,'single');
fclose(fHbO);
fHbR = fopen([FolderName filesep 'HbR.dat'], 'r+');
dat = fread(fHbR,inf,'single');
frewind(fHbR);
dat = reshape(dat, NbFrames, iWidth, iHeight);
dat = permute(dat, [3 2 1]);
fwrite(fHbR, dat,'single');
fclose(fHbR);
end
