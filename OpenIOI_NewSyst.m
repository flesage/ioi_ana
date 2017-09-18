function out = OpenIOI_NewSyst(FolderName, Binning, Version)
%%%%DEFINE -> THIS VALUE IS HARDCODED!!!!
NOFPF = 256;
DEF_VISUEL = 0;
if( exist('Config.m','file') )
    Config;
end

AcqInfoStream = readtable([FolderName filesep 'info.txt'],...
    'Delimiter',':','ReadVariableNames',false, 'ReadRowNames',true);

disp('Computing stimulation parameters')
disp('**************************');
if( DEF_VISUEL )
    tStim = AcqInfoStream{'Stimulation',1};
    tAIChan = AcqInfoStream{'AINChannels',1};
else
    if( Version > 1)
        tStim = AcqInfoStream{'Stimulation1',1};
    else
        tStim = AcqInfoStream{'Stimulation',1};
    end
    tAIChan = AcqInfoStream{'AINChannels',1};
end
if( iscell(tStim) )
    tStim = str2double(cell2mat(tStim));
    tAIChan =  str2double(cell2mat(tAIChan));
end
if( tStim )
    IOIReadStimFile_NS(FolderName);
else
    fprintf('No stimulation detected. \n');
end

imgFilesList = dir([FolderName filesep 'img_0*.bin']);
aiFilesList = dir([FolderName filesep 'ai_*.bin']);

%Version Management Here:
if( Version == 3)
    hWima = 5;
    hWai = 5;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint64', 3, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
    ImRes_XY = [nx, ny];
    SizeImage = nx*ny*2 + 3*8;
elseif( Version == 2)
    hWima = 5;
    hWai = 5;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint64', 1, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
    ImRes_XY = [nx, ny];
    SizeImage = nx*ny*2 + 8;
elseif( Version == 1 )
    %TODO: To be validated
    hWima = 4;
    hWai = 5;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint16', [double(nx), double(ny)], 'imgj'};
    ImRes_XY = [nx, ny];
    SizeImage = nx*ny*2 + 8;
else
    disp(['Error! System version ' int2str(Version) ' is not suported by this software']);
end


NombreImage = 0;
for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    NombreImage = NombreImage+size(data.Data,1);
end

if( Version == 3 )
   idImg = zeros(NombreImage, 3);
else
   idImg = zeros(NombreImage, 1);
end

for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    idImg((NOFPF*(ind-1)+1):(NOFPF*(ind-1)+size(data.Data,1)),:) = cell2mat(arrayfun(@(x) data.Data(x).framej, 1:size(data.Data,1),'UniformOutput',false))';
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

AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', hWai*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, 1e4, tAIChan, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],tAIChan);
    AnalogIN = [AnalogIN; tmp];
end
clear tmp ind data;

%%%%
%Stimulation Params
%%%%
if( DEF_VISUEL )
    FrameTiming = 0:0.2:NombreImage*0.2;
    StimTiming = 0:1e-4:length(AnalogIN(:,1))*1e-4;
    Stim = zeros(size(FrameTiming));
    for ind = 1:length(FrameTiming)
        idx = find(StimTiming <= FrameTiming(ind), 1, 'last');
        Stim(ind) = AnalogIN(idx,1);
    end
    bStim = 1;
    Stim = [0, Stim];
    Stim = Stim';
else
    if( exist([FolderName filesep 'StimParameters.mat'], 'file') )
        load([FolderName filesep 'StimParameters.mat']);
        disp('Stim detected: yes');
        disp(['Number of events: ' int2str(NbStim)]);
        disp(['Length of each event: ' int2str(StimLength) 'sec']);
        bStim = 1;
    else
        bStim = 0;
    end
end


if( ~DEF_VISUEL )
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
if( DEF_VISUEL )
    bFluo = (tColor > 3); tColor = mod(tColor,4);
    bYellow = (tColor > 1); tColor = mod(tColor,2);
    bRed = (tColor > 0);
    bGreen = 0;
    clear fInfo tColor;
else
    bFluo = (tColor > 7); tColor = mod(tColor,8);
    bGreen = (tColor > 3); tColor = mod(tColor,4);
    bYellow = (tColor > 1); tColor = mod(tColor,2);
    bRed = (tColor > 0);
    clear fInfo tColor;
end

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
    if( DEF_VISUEL )
        disp('Fluo illumination detected');
        if( exist([FolderName filesep 'Data_Fluo.mat'],'file') )
            delete([FolderName filesep 'Data_Fluo.mat']);
        end
        fSpeckle = matfile([FolderName filesep 'Data_Fluo.mat'],'Writable',true);
        fSpeckle.datFile = [FolderName filesep 'fChan.dat'];
        fSpeckle.datSize = [Rx, Ry];
        fSpeckle.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
        fSpeckle.Freq = Freq/nbColors;
        cSpeckle = 1;
        fidS = fopen([FolderName filesep 'fChan.dat'],'w');
    else
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
    end
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
[~, idxOri] = unique(idImg(:,1));
%badFrames = find(~ismember(1:NombreImage, uniqueFramesID));
Conseq = conv(diff(idImg(:,1)), [1 1 1]) == 3;
idxS = find(diff(Conseq,1,1)==-1) + 1; idxS(end) = [];
idxE = find(diff(Conseq,1,1)==1); idxE(1) = [];
badFrames = [];
for ind = 1:length(idxS)
    badFrames = [badFrames, (idImg(idxS(ind)) + 1):(idImg(idxE(ind)) - 1)];
end
idxOri(ismember(idImg(idxOri,1),badFrames)) = [];

%%% Lookup Table For missing frames
if( ~isempty(badFrames) )
    InterpLUT = zeros(8,size(badFrames,2));
    % 1: Frame before
    % 2: File Number where to find this frame
    % 3: Image ID in this file
    % 4: Frame After
    % 5: File Number where to find this frame
    % 6: Image ID in this file
    % 7: Ratio between these frame and the one missing
    % 8: Frame tag id
    for ind = 1:size(badFrames,2)
        tmpID = badFrames(ind);
        tmpBefore = tmpID - (nbColors:nbColors:(tmpID-1));
        idx = find(ismember(tmpBefore,idImg(:,1))&~ismember(tmpBefore,badFrames),1,'first');
        tmpBefore = tmpBefore(idx);
        InterpLUT(1,ind) = tmpBefore;
        idx = find(tmpBefore == idImg);
        InterpLUT(2,ind) = floor((idx-1)/256) + 1;
        InterpLUT(3,ind) = rem((idx-1),256) + 1;
        
        tmpAfter = tmpID + (nbColors:nbColors:(NombreImage));
        idx = find(ismember(tmpAfter,idImg(:,1))&~ismember(tmpAfter,badFrames),1,'first');
        tmpAfter = tmpAfter(idx);
        InterpLUT(4,ind) = tmpAfter;
        idx = find(tmpAfter == idImg);
        InterpLUT(5,ind) = floor((idx-1)/256) + 1;
        InterpLUT(6,ind) = rem((idx-1),256) + 1;
        
        tmpRatio =  (tmpID - idImg(InterpLUT(1,ind),1))./...
            (idImg(InterpLUT(4,ind),1) - idImg(InterpLUT(1,ind),1));
        InterpLUT(7,ind) = tmpRatio;
        InterpLUT(8,ind) = badFrames(ind);
    end
    clear tmpRatio tmpAfter tmpBefore;
    %%% Interpolation of missing frames
    TmpFrames = struct('framej',[], 'imgj',[]);
    for ind = 1:size(InterpLUT,2)
        dBefore = memmapfile([FolderName filesep...
            imgFilesList(InterpLUT(2,ind)).name],...
            'Offset', hWima*4 + (InterpLUT(3,ind) - 1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1);
        dAfter = memmapfile([FolderName filesep...
            imgFilesList(InterpLUT(5,ind)).name],...
            'Offset', hWima*4 + (InterpLUT(6,ind) - 1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1);
        TmpFrames(ind).imgj = uint16(round(InterpLUT(7,ind)*double(dAfter.Data.imgj - dBefore.Data.imgj))) + dBefore.Data.imgj;
        TmpFrames(ind).framej = uint64([InterpLUT(8,ind), 1, 1]);
    end
    
    fid = fopen([FolderName filesep 'img_interp.bin'],'w');
    for ind = 1:size(InterpLUT,2)
        fwrite(fid, TmpFrames(ind).framej, 'uint64');
        fwrite(fid, TmpFrames(ind).imgj, 'uint16');
    end
    fclose(fid);
end

%Rebuilding addresses for each frames...
ImAddressBook = zeros(NombreImage,2);
for ind = 1:NombreImage
    if( ismember(ind, InterpLUT(8,:)) )
        fidx = find( ind == InterpLUT(8,:), 1, 'first');
        ImAddressBook(ind,1) = size(imgFilesList,1) + 1;
        ImAddressBook(ind,2) = fidx;
    elseif( ismember(ind, idImg(idxOri)) )
        fidx = find( ind == idImg, 1, 'first');
        ImAddressBook(ind,1) = floor((fidx-1)/256) + 1;
        ImAddressBook(ind,2) = rem(fidx-1, 256) + 1;
    end
end

%Saving infos...
save([FolderName 'ImagesLUT.mat'], 'ImAddressBook');


%%%%
% Images Classification and filtering
%%%%
ind = 1;
if( bFluo )
    tags = ind:nbColors:NombreImage;
    Images = zeros(ImRes_XY(1), ImRes_XY(2), length(tags), 'single');
    for indI = 1:length(tags)
        indF = tags(indI);
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
               
        Images(:,:,indI) = dat.Data.imgj;
        if( bStim )
            fSpeckle.Stim(cSpeckle,1) = Stim(indF);
        else
            fSpeckle.Stim(cSpeckle,1) = 0;
        end
        cSpeckle = cSpeckle + 1;
    end
    fwrite(fidS, Images, 'single');
    ind = ind + 1;
    fSpeckle.datLength = cSpeckle - 1;
    fclose(fidS);
end
if( bRed )
    tags = ind:nbColors:NombreImage;
    Images = zeros(ImRes_XY(1), ImRes_XY(2), length(tags), 'single');
  
    for indI = 1:length(tags)
        indF = tags(indI);
        
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
         
        Images(:,:,indI) = dat.Data.imgj;
        if( bStim )
            fRed.Stim(cRed,1) = Stim(indF);
        else
            fRed.Stim(cRed,1) = 0;
        end
        cRed = cRed + 1;
    end
    fwrite(fidR, Images, 'single');
    ind = ind + 1;
    fRed.datLength = cRed - 1;
    fclose(fidR);
end
if( bYellow )
    tags = ind:nbColors:NombreImage;
    Images = zeros(ImRes_XY(1), ImRes_XY(2), length(tags), 'single');
   
    for indI = 1:length(tags)
        indF = tags(indI);
        
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
       
        
        Images(:,:,indI) = dat.Data.imgj;
        if( bStim )
            fYellow.Stim(cYellow,1) = Stim(indF);
        else
            fYellow.Stim(cYellow,1) = 0;
        end
        cYellow = cYellow + 1;
    end
    fwrite(fidY, Images, 'single');
    ind = ind + 1;
    fYellow.datLength = cYellow - 1;
    fclose(fidY);
end
if( bGreen )
    tags = ind:nbColors:NombreImage;
    Images = zeros(ImRes_XY(1), ImRes_XY(2), length(tags), 'single');
  
    for indI = 1:length(tags)
        indF = tags(indI);
   
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
       
        Images(:,:,indI) = dat.Data.imgj;
        if( bStim )
            fGreen.Stim(cGreen,1) = Stim(indF);
        else
            fGreen.Stim(cGreen,1) = 0;
        end
        cGreen = cGreen + 1;
    end
    fwrite(fidG, Images, 'single');
    ind = ind + 1;
    fGreen.datLength = cGreen - 1;
    fclose(fidG);
end

fprintf('\n');
%Verbose
disp(['Done with file ' FolderName]);
str = ['************* ' sprintf('\r')];
disp(str);
%end of Verbose
end
