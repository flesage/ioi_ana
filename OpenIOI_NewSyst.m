function out = OpenIOI_NewSyst(FolderName, Binning)

imgFilesList = dir([FolderName filesep 'img_*.bin']);
aiFilesList = dir([FolderName filesep 'ai_*.bin']);

header = memmapfile([FolderName filesep imgFilesList(1).name], ...
    'Offset', 0, 'Format', {'int32', 5, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);

nx=header.Data.header(2);
ny=header.Data.header(3);
ImRes_XY = [nx, ny];
SizeImage = nx*ny*2 + 8;
NombreImage = 0;
for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',5*4,'Format',{'uint64', 1, 'framej';'uint16', [1024,1024], 'imgj'},'repeat',inf);
    NombreImage = NombreImage+size(data.Data,1);
end
idImg = zeros(NombreImage, 1);
for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',5*4,'Format',{'uint64', 1, 'framej';'uint16', [1024,1024], 'imgj'},'repeat',inf);
    idImg((256*(ind-1)+1):(256*(ind-1)+size(data.Data,1))) = arrayfun(@(x) data.Data(x).framej, 1:size(data.Data,1));
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
    %TODO....
end

AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, 1e4, 11, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],11);
    AnalogIN = [AnalogIN; tmp];
end
clear tmp ind data;

CamTrig = find((AnalogIN(1:(end-1),1) < 2.5) & (AnalogIN(2:end,1) >= 2.5));
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
fInfo = fopen([FolderName filesep 'info.txt']);
formatSpec = '%*s FrameRateHz: %d Width:%d Height:%d ExposureMsec:%f AISampleRate:%d AINChannels:%d Stimulation:%d Illumination:%d';
AcqInfos =textscan(fInfo,formatSpec, 'Delimiter', '\n');
fclose(fInfo);
tColor = AcqInfos{8};
bFluo = (tColor > 8); tColor = mod(tColor,8);
bGreen = (tColor > 4); tColor = mod(tColor,4);
bYellow = (tColor > 2); tColor = mod(tColor,2);
bRed = (tColor > 0); 
clear fInfo tColor;

if( Binning )
    Rx = round(ImRes_XY(1)/2);
    Ry = round(ImRes_XY(2)/2);
else
    Rx = ImRes_XY(1);
    Ry = ImRes_XY(2);
end

nbColors = (bFluo + bGreen + bYellow + bRed);
if( bRed )
    disp('Red illumination detected');
    if( exist([FolderName filesep 'Data_red.mat'],'file') )
        delete([FolderName filesep 'Data_red.mat']);
    end
    fRed = matfile([FolderName filesep 'Data_red.mat'],'Writable',true);
    fRed.datFile = [FolderName filesep 'rChan.dat'];
    fRed.datSize = [Rx, Ry];
    fRed.Stim = zeros(1, floor(NombreImage/nbColors), 'single');
    cRed = 1;
    fidR = fopen([FolderName filesep 'rChan.dat'],'w');
end
if( bGreen )
    disp('Green illumination detected');
    if( exist([FolderName filesep 'Data_green.mat'],'file') )
        delete([FolderName filesep 'Data_green.mat']);
    end
    fGreen = matfile([FolderName filesep 'Data_green.mat'],'Writable',true);
    fGreen.datFile = [FolderName filesep 'gChan.dat'];
    fGreen.datSize = [Rx, Ry];
    fGreen.Stim = zeros(1, floor(NombreImage/nbColors), 'single');
    cGreen = 1;
    fidG = fopen([FolderName filesep 'gChan.dat'],'w');
end
if( bYellow )
    disp('Yellow illumination detected');
    if( exist([FolderName filesep 'Data_yellow.mat'],'file') )
        delete([FolderName filesep 'Data_yellow.mat']);
    end
    fYellow = matfile([FolderName filesep 'Data_yellow.mat'],'Writable',true);
    fYellow.datFile = [FolderName filesep 'yChan.dat'];
    fYellow.datSize = [Rx, Ry];
    fYellow.Stim = zeros(1, floor(NombreImage/nbColors), 'single');
    cYellow = 1;
    fidY = fopen([FolderName filesep 'yChan.dat'],'w');
end
if( bFluo )
    disp('Speckle illumination detected');
    if( exist([FolderName filesep 'Data_speckle.mat'],'file') )
        delete([FolderName filesep 'Data_speckle.mat']);
    end
    fSpeckle = matfile([FolderName filesep 'Data_speckle.mat'],'Writable',true);
    fSpeckle.datFile = [FolderName filesep 'sChan.dat'];
    fSpeckle.datSize = [Rx, Ry];
    fSpeckle.Stim = zeros(1, floor(NombreImage/nbColors), 'single');
    cSpeckle = 1;
    fidS = fopen([FolderName filesep 'sChan.dat'],'w');
end

% Verbose
if( sum(Stim(:)) > 0 )
    disp('Stim detected: yes');
    disp(['Number of events: ' int2str(NbStim)]);
    disp(['Length of each event: ' int2str(StimLength) 'sec']);
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
        'Offset',5*4, 'Format', {'uint64', 1, 'framej';'uint16', [1024,1024], 'imgj'},...
        'repeat',inf);

    Frame = data.Data;
    Frame = reshape([Frame(:).imgj],ImRes_XY(1),ImRes_XY(2),[]);
        
    if( Binning )
        Frame = imresize(single(Frame), 0.5);
    else
        Frame = single(Frame);
    end

    for indF = 1:size(Frame,3)
        Frame(:,:,indF) = single(xRemoveStripesVertical(squeeze(Frame(:,:,indF)), 8, 'db4', 2));
    end
    
    iFrame = 1;
    if( bFluo )
        idx = iFrame:nbColors:size(data.Data,1);
        iFrame = iFrame + 1;
        fwrite(fidS, Frame(:, :, idx), 'single');
        fSpeckle.Stim(1, cSpeckle:(cSpeckle + length(idx) - 1)) = Stim(idx + (ind-1)*256);
    end
    if( bRed )
        idx = iFrame:nbColors:size(data.Data,1);
        iFrame = iFrame + 1;
        fwrite(fidR, Frame(:, :, idx), 'single');
        fRed.Stim(1, cRed:(cRed + length(idx) - 1)) = Stim(idx + (ind-1)*256);
    end
    if( bYellow )
        idx = iFrame:nbColors:size(data.Data,1);
        iFrame = iFrame + 1;
        fwrite(fidY, Frame(:, :, idx), 'single');
        fYellow.Stim(1, cYellow:(cYellow + length(idx) - 1)) = Stim(idx + (ind-1)*256);
    end
    if( bGreen )
        idx = iFrame:nbColors:size(data.Data,1);
        iFrame = iFrame + 1;
        fwrite(fidG, Frame(:, :, idx), 'single');
        fGreen.Stim(1, cGreen:(cGreen + length(idx) - 1)) = Stim(idx + (ind-1)*256);
    end
  
    fprintf('%d%%...', flagPrc(indPrc));
    indPrc = indPrc + 1;
end
if( bFluo )
    fSpeckle.datLength = cSpeckle;
end
if( bRed )
    fRed.datLength = cRed;
end
if( bYellow )
    fYellow.datLength = cYellow;
end
if( bGreen )
    fGreen.datLength = cGreen;
end
  
fprintf('\n');
%Verbose
disp(['Done with file ' FolderName]);
str = ['************* ' sprintf('\r')];
disp(str);
%end of Verbose
end
