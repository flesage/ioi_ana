function out = OpenIOI_NewSyst(FolderName, Binning, Version)

%%%%DEFINE -> THIS VALUE IS HARDCODED!!!!
NOFPF = 256;

disp('Computing stimulation parameters')
disp('**************************');
IOIReadStimFile_NS(FolderName);

imgFilesList = dir([FolderName filesep 'img_*.bin']);
aiFilesList = dir([FolderName filesep 'ai_*.bin']);

AcqInfoStream = readtable([FolderName filesep 'info.txt'],...
    'Delimiter',':','ReadVariableNames',false, 'ReadRowNames',true);

%Version Management Here:
if( Version == 2)
    hWima = 5;
    hWai = 5;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint64', 1, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
elseif( Version == 1 )
    %TODO: To be validated
    hWima = 4;
    hWai = 5;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint16', [double(nx), double(ny)], 'imgj'};
else
    disp(['Error! System version ' int2str(Version) ' is not suported by this software']);
end

ImRes_XY = [nx, ny];
SizeImage = nx*ny*2 + 8;
NombreImage = 0;
for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    NombreImage = NombreImage+size(data.Data,1);
end
idImg = zeros(NombreImage, 1);
for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    idImg((NOFPF*(ind-1)+1):(NOFPF*(ind-1)+size(data.Data,1))) = arrayfun(@(x) data.Data(x).framej, 1:size(data.Data,1));
end

clear nx ny data header ind;

% Verbose
disp(['Opening of: ' FolderName]);
disp(['Number of Frames acquired: ' int2str(NombreImage)]);
disp(['Frames'' resolution: ' int2str(ImRes_XY(1)) ' pix X ' int2str(ImRes_XY(2)) ' pix']);
% end of Verbose

%If binning...
if( Binning )
    %Verbose
    disp('Binning option is ON');
    %end of Verbose
end

%%%%
%Stimulation Params
%%%%
if( exist([FolderName filesep 'StimParameters.mat'], 'file') )
    load([FolderName filesep 'StimParameters.mat']);
else
    disp('Something went wrong!');
    return;
end

AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', hWai*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, 1e4, 11, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],11);
    AnalogIN = [AnalogIN; tmp];
end
clear tmp ind data;

CamTrig = find((AnalogIN(1:(end-1),1) < 2.5) & (AnalogIN(2:end,1) >= 2.5))+1;
StartDelay = round(CamTrig(1)/10);
EndDelay = round((length(AnalogIN(:,1)) - CamTrig(end))/10);

% Verbose
disp(['Camera Trigs detected: ' int2str(length(CamTrig))]);
disp(['Recording of analog inputs starts ' int2str(StartDelay) ' ms before the first trigger.']) 
disp(['Recording of analog inputs ends ' int2str(EndDelay) ' ms after the last trigger.']) 
% end of Verbose
clear StartDelay EndDelay
%Missing frames?
%Not taken care of for now...

%Less trig than images... something's wrong!
if( length(CamTrig) < NombreImage  ) 
    disp('IOI Error: Analog recordings and Image files don''t match. Impossible to continue further.');
    out = 'Error';
    return
end

%%%%
% Color Sequence
%%%%
%Red                   = 0001;
%Yellow/Amber = 0010;
%Green               = 0100;
%Fluo/Laser       = 1000;

tColor = AcqInfoStream{'Illumination',1};
if( iscell(tColor) )
    tColor = str2double(cell2mat(tColor));
end
bFluo = (tColor > 7); tColor = mod(tColor,8);
bGreen = (tColor > 3); tColor = mod(tColor,4);
bYellow = (tColor > 1); tColor = mod(tColor,2);
bRed = (tColor > 0); 
clear fInfo tColor;

if( Binning )
    Rx = round(ImRes_XY(1)/2);
    Ry = round(ImRes_XY(2)/2);
else
    Rx = ImRes_XY(1);
    Ry = ImRes_XY(2);
end
Freq = AcqInfoStream{'FrameRateHz',1};
if( iscell(Freq) )
    Freq = str2double(cell2mat(Freq));
end

nbColors = (bFluo + bGreen + bYellow + bRed);
idx = 1;
if( bFluo )
    disp('Speckle illumination detected');
    if( exist([FolderName filesep 'Data_speckle.mat'],'file') )
        delete([FolderName filesep 'Data_speckle.mat']);
    end
    fSpeckle = matfile([FolderName filesep 'Data_speckle.mat'],'Writable',true);
    fSpeckle.datFile = [FolderName filesep 'sChan.dat'];
    fSpeckle.datSize = [Rx, Ry];
    fSpeckle.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
    fSpeckle.Freq = Freq/nbColors;
    cSpeckle = 1;
    fidS = fopen([FolderName filesep 'sChan.dat'],'w');
    sExpectedID = idx:nbColors:NombreImage;
    idx = idx + 1;
end
if( bRed )
    disp('Red illumination detected');
    if( exist([FolderName filesep 'Data_red.mat'],'file') )
        delete([FolderName filesep 'Data_red.mat']);
    end
    fRed = matfile([FolderName filesep 'Data_red.mat'],'Writable',true);
    fRed.datFile = [FolderName filesep 'rChan.dat'];
    fRed.datSize = [Rx, Ry];
    fRed.Stim = zeros(floor(NombreImage/nbColors), 1, 'single');
    fRed.Freq = Freq/nbColors;
    cRed = 1;
    fidR = fopen([FolderName filesep 'rChan.dat'],'w');
    rExpectedID = idx:nbColors:NombreImage;
    idx = idx + 1;
end
if( bYellow )
    disp('Yellow illumination detected');
    if( exist([FolderName filesep 'Data_yellow.mat'],'file') )
        delete([FolderName filesep 'Data_yellow.mat']);
    end
    fYellow = matfile([FolderName filesep 'Data_yellow.mat'],'Writable',true);
    fYellow.datFile = [FolderName filesep 'yChan.dat'];
    fYellow.datSize = [Rx, Ry];
    fYellow.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
    fYellow.Freq = Freq/nbColors;
    cYellow = 1;
    fidY = fopen([FolderName filesep 'yChan.dat'],'w');
    yExpectedID = idx:nbColors:NombreImage;
    idx = idx + 1;
end
if( bGreen )
    disp('Green illumination detected');
    if( exist([FolderName filesep 'Data_green.mat'],'file') )
        delete([FolderName filesep 'Data_green.mat']);
    end
    fGreen = matfile([FolderName filesep 'Data_green.mat'],'Writable',true);
    fGreen.datFile = [FolderName filesep 'gChan.dat'];
    fGreen.datSize = [Rx, Ry];
    fGreen.Stim = zeros(floor(NombreImage/nbColors), 1, 'single');
    fGreen.Freq = Freq/nbColors;
    cGreen = 1;
    fidG = fopen([FolderName filesep 'gChan.dat'],'w');
    gExpectedID = idx:nbColors:NombreImage;
    idx = idx + 1;
end

%Interpolation for bad or missing frames
disp(['Number of bad frames: ' int2str(sum(diff(idImg,1,1) == 0)) ' (' num2str(100*sum(diff(idImg,1,1) == 0)/NombreImage) '%)']);
uniqueFramesID = unique(idImg);
badFrames = find(~ismember(1:NombreImage, uniqueFramesID));
fprintf('Missing IDs: ');
arrayfun(@(x) fprintf('%d, ', x), badFrames); fprintf('\n');
fFrames = 1:nbColors;
lFrames = (NombreImage-nbColors+1):NombreImage;
if( any(ismember(fFrames, badFrames)) || any(ismember(lFrames, badFrames)) )
    fprintf('Extrapolating: ');
    extList = fFrames(ismember(fFrames, badFrames));
    if( ~isempty(extList) )
        for ind = 1:length(extList)
            fprintf(' #%d', extList(ind));
            iAfter = extList(ind);
            while( ismember(iAfter, badFrames) )
                iAfter = iAfter + nbColors;
            end
            iVeryAfter = iAfter + nbColors;
            while( ismember(iVeryAfter, badFrames) )
                iVeryAfter = iVeryAfter + nbColors;
            end
            
            fileNumber = floor(iAfter/NOFPF);
            fName = sprintf('img_%05d.bin',fileNumber);
            fAfter = memmapfile([FolderName filesep fName],...
                'Offset', 5*4 + (mod(iAfter,NOFPF)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
            
            fileNumber = floor(iVeryAfter/NOFPF);
            fName = sprintf('img_%05d.bin',fileNumber);
            fVeryAfter = memmapfile([FolderName filesep fName],...
                'Offset', 5*4 + (mod(iVeryAfter,NOFPF)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
            
            ratio = abs(extList(ind) - iVeryAfter)/abs(iAfter - iVeryAfter);
            eFrame = fVeryAfter.Data.imgj + ratio*(fAfter.Data.imgj - fVeryAfter.Data.imgj);
            
            fileNumber = floor(extList(ind)/NOFPF);
            fName = sprintf('img_%05d.bin',fileNumber);
            fExtra = memmapfile([FolderName filesep fName],...
                'Offset', 5*4 + (mod(extList(ind),NOFPF)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1, 'Writable', true);
            fExtra.Data.framej = uint64(extList(ind));
            fExtra.Data.imgj = eFrame;
        end
        clear iAfter iVeryAfter fileNumber fName fAfter fVeryAfter
        clear fExtra eFrame ratio extList ind lFrames
    end
    extList = lFrames(ismember(lFrames, badFrames));
    if( ~isempty(extList) )
        for ind = 1:length(extList)
            fprintf(' #%d', extList(ind));
            iBefore = extList(ind);
            while( ismember(iBefore, badFrames) )
                iBefore = iBefore - nbColors;
            end
            iVeryBefore = iBefore - nbColors;
            while( ismember(iVeryBefore, badFrames) )
                iVeryBefore = iVeryBefore - nbColors;
            end
            
            fileNumber = floor(iBefore/NOFPF);
            fName = sprintf('img_%05d.bin',fileNumber);
            fBefore = memmapfile([FolderName filesep fName],...
                'Offset', 5*4 + (mod(iBefore,NOFPF)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
            
            fileNumber = floor(iVeryBefore/NOFPF);
            fName = sprintf('img_%05d.bin',fileNumber);
            fVeryBefore = memmapfile([FolderName filesep fName],...
                'Offset', 5*4 + (mod(iVeryBefore,NOFPF)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
            
            ratio = (extList(ind) - iVeryBefore)/(iBefore - iVeryBefore);
            eFrame = fVeryBefore.Data.imgj + ratio*(fBefore.Data.imgj - fVeryBefore.Data.imgj);
            
            fileNumber = floor(extList(ind)/NOFPF);
            fName = sprintf('img_%05d.bin',fileNumber);
            fExtra = memmapfile([FolderName filesep fName],...
                'Offset', 5*4 + (mod(extList(ind),NOFPF)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1, 'Writable', true);
            fExtra.Data.framej = uint64(extList(ind));
            fExtra.Data.imgj = eFrame;
        end
        clear iBefore iVeryBefore fileNumber fName fBefore fVeryBefore
        clear fExtra eFrame ratio extList ind
    end
end
clear fFrames lFrames

mFrames = (nbColors+1):(NombreImage-nbColors);
if( any(ismember(mFrames, badFrames)) )
    fprintf('Interpolating: ');
    intList = mFrames(ismember(mFrames, badFrames));
    for ind = 1:length(intList)
        fprintf(' #%d', intList(ind));
        iAfter = intList(ind);
        while( ismember(iAfter, badFrames) )
            iAfter = iAfter + nbColors;
        end
        iBefore = intList(ind);
        while( ismember(iBefore, badFrames) )
            iBefore = iBefore - nbColors;
        end
        
        fileNumber = floor(iAfter/NOFPF);
        fName = sprintf('img_%05d.bin',fileNumber);
        fAfter = memmapfile([FolderName filesep fName],...
            'Offset', 5*4 + (mod(iAfter,NOFPF)-1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1);
        
        fileNumber = floor(iBefore/NOFPF);
        fName = sprintf('img_%05d.bin',fileNumber);
        fBefore = memmapfile([FolderName filesep fName],...
            'Offset', 5*4 + (mod(iBefore,NOFPF)-1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1);
        
        ratio = (intList(ind) - iBefore)/(iAfter - iBefore);
        eFrame = fBefore.Data.imgj + ratio*(fAfter.Data.imgj - fBefore.Data.imgj);
        
        fileNumber = floor(intList(ind)/NOFPF);
        fName = sprintf('img_%05d.bin',fileNumber);
        fIntra = memmapfile([FolderName filesep fName],...
            'Offset', 5*4 + (mod(intList(ind),NOFPF)-1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1, 'Writable', true);
        fIntra.Data.framej = uint64(intList(ind));
        fIntra.Data.imgj = eFrame;
    end
    clear iAfter iBefore fileNumber fName fAfter fBefore
    clear fIntra eFrame ratio intList ind mFrames
end
  
% Verbose
if( sum(Stim(:)) > 0 )
    disp('Stim detected: yes');
    disp(['Number of events: ' int2str(NbStim)]);
    disp(['Length of each event: ' int2str(StimLength) 'sec']);
    bStim = 1;
else
    bStim = 0;
end
% end of Verbose

%%%%
% Images Classification and filtering
%%%%
Marks = size(imgFilesList,1);
indPrc = 2;
flagPrc = round(linspace(1,100,Marks+1));
fprintf('Progress: ');
for ind = 1:Marks
    data = memmapfile([FolderName filesep imgFilesList(ind).name],...
        'Offset',5*4, 'Format', frameFormat, 'repeat', inf);

    fID = arrayfun(@(x) data.Data(x).framej, 1:size(data.Data,1));
    Frame = data.Data;
    Frame = reshape([Frame(:).imgj],ImRes_XY(1),ImRes_XY(2),[]);
        
    if( Binning )
        Frame = imresize(single(Frame), 0.5);
    else
        Frame = single(Frame);
    end

    for indF = 1:size(Frame,3)
        Frame(:,:,indF) = single(xRemoveStripesVertical(single(squeeze(Frame(:,:,indF))), 8, 'db4', 2));
    end
       
    iFrame = 1;
    if( bFluo )
        idx = iFrame:nbColors:size(data.Data,1);
        iFrame = iFrame + 1;
        
        fwrite(fidS, Frame(:, :, idx), 'single');
        if( bStim )
            fSpeckle.Stim(cSpeckle:(cSpeckle + length(idx) - 1),1) = Stim(idx + (ind-1)*NOFPF);
        else
            fSpeckle.Stim(cSpeckle:(cSpeckle + length(idx) - 1),1) = zeros(length(idx),1);
        end
        cSpeckle = (cSpeckle + length(idx));
    end
    if( bRed )
        idx = iFrame:nbColors:size(data.Data,1);
        iFrame = iFrame + 1;
              
        fwrite(fidR, Frame(:, :, idx), 'single');
        if( bStim )
            fRed.Stim(cRed:(cRed + length(idx) - 1),1) = Stim(idx + (ind-1)*NOFPF);
        else
            fRed.Stim(cRed:(cRed + length(idx) - 1),1) = zeros(length(idx),1);
        end
        cRed = (cRed + length(idx));
    end
    if( bYellow )
        idx = iFrame:nbColors:size(data.Data,1);
        iFrame = iFrame + 1;
               
        fwrite(fidY, Frame(:, :, idx), 'single');
        if( bStim )
            fYellow.Stim(cYellow:(cYellow + length(idx) - 1),1) = Stim(idx + (ind-1)*NOFPF);
        else
            fYellow.Stim(cYellow:(cYellow + length(idx) - 1),1) =  zeros(length(idx),1);
        end
        cYellow = (cYellow + length(idx));
    end
    if( bGreen )
        idx = iFrame:nbColors:size(data.Data,1);
        iFrame = iFrame + 1;
               
        fwrite(fidG, Frame(:, :, idx), 'single');
        if( bStim )
            fGreen.Stim(cGreen:(cGreen + length(idx) - 1),1) = Stim(idx + (ind-1)*NOFPF);
        else
            fGreen.Stim(cGreen:(cGreen + length(idx) - 1),1) =  zeros(length(idx),1);
        end
        cGreen = (cGreen + length(idx));
    end
  
    fprintf('%d%%...', flagPrc(indPrc));
    indPrc = indPrc + 1;
    if( mod(indPrc,15) == 0 )
        fprintf('\n');
    end
end
if( bFluo )
    fSpeckle.datLength = cSpeckle;
    fclose(fidS);
end
if( bRed )
    fRed.datLength = cRed;
    fclose(fidR);
end
if( bYellow )
    fYellow.datLength = cYellow;
    fclose(fidY);
end
if( bGreen )
    fGreen.datLength = cGreen;
    fclose(fidG);
end
  
fprintf('\n');
%Verbose
disp(['Done with file ' FolderName]);
str = ['************* ' sprintf('\r')];
disp(str);
%end of Verbose
end
