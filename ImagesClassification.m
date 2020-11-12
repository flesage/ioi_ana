function out = ImagesClassification(FolderName, BinningSpatial, BinningTemp, b_SubROI, b_IgnoreStim)

%%%%DEFINES -> THESE CONSTANTS ARE HARDCODED!!!! DO NOT CHANGE THEM.
DEF_VISUEL = 0;

%%%%%%%%%%%%%%%%%%%%%
% Acq. Info file:
%%%%%%%%%%%%%%%%%%%%%
AcqInfoStream = ReadInfoFile(FolderName);
if( ~isfield(AcqInfoStream, 'Camera_Model') )
   AcqInfoStream.Camera_Model = 'CS2100M'; 
end

disp('Recovering stimulation parameters')
disp('**************************');

tAIChan = AcqInfoStream.AINChannels;

%%%%%%%%%%%%%%%%%%%%%
% Stimulation detected
%%%%%%%%%%%%%%%%%%%%%
if( ~b_IgnoreStim )
    Fields = fieldnames(AcqInfoStream);
    idx = contains(Fields, 'stimulation','IgnoreCase',true);
    Fields = Fields(idx);
    cnt = 0;
    for indS = 1:length(Fields)
        if(isnumeric(AcqInfoStream.(Fields{indS})) && ...
                AcqInfoStream.(Fields{indS}) > 0)
            AcqInfoStream.Stimulation = 1;
            ReadAnalogsIn(FolderName, AcqInfoStream);
            cnt = cnt + 1;
            break;
        end
    end
    if( cnt == 0 )
        fprintf('Stimulation not detected. \n');
    end
else
    fprintf('Stimulation ignored. \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Images sequence validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgFilesList = dir([FolderName filesep 'img_0*.bin']);
if( AcqInfoStream.MultiCam )
    imgFilesList2 = dir([FolderName filesep 'imgCam2_0*.bin']);
end 
aiFilesList = dir([FolderName filesep 'ai_*.bin']);

hWima = 5;
hWai = 5;
header = memmapfile([FolderName filesep imgFilesList(1).name], ...
    'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);

Version = header.Data.header(1);
nx=header.Data.header(2);
ny=header.Data.header(3);
FrameSz = header.Data.header(4);
NbImsPefFile = single(header.Data.header(5));

frameFormat = {'uint64', 3, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
ImRes_XY = [nx, ny];
SizeImage = nx*ny*2 + 3*8;

NombreImage = 0;
for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    NombreImage = NombreImage+size(data.Data,1);
end
idImg = zeros(NombreImage, 3);

for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    idImg((NbImsPefFile*(ind-1)+1):(NbImsPefFile*(ind-1)+size(data.Data,1)),:) = cell2mat(arrayfun(@(x) data.Data(x).framej, 1:size(data.Data,1),'UniformOutput',false))';
end

if( AcqInfoStream.MultiCam )
    NombreImage2 = 0;
    for ind = 1:size(imgFilesList2,1)
        data = memmapfile([FolderName filesep imgFilesList2(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
        NombreImage2 = NombreImage2 + size(data.Data,1);
    end
    idImg2 = zeros(NombreImage2, 3);

    for ind = 1:size(imgFilesList2,1)
        data = memmapfile([FolderName filesep imgFilesList2(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
        idImg2((NbImsPefFile*(ind-1)+1):(NbImsPefFile*(ind-1)+size(data.Data,1)),:) = cell2mat(arrayfun(@(x) data.Data(x).framej, 1:size(data.Data,1),'UniformOutput',false))';
    end
end
clear nx ny data header ind;

% Verbose
disp(['Opening of: ' FolderName]);
disp(['Number of Frames acquired: ' int2str(NombreImage)]);
disp(['Frames'' resolution: ' int2str(ImRes_XY(1)) ' pix X ' int2str(ImRes_XY(2)) ' pix']);
if( AcqInfoStream.MultiCam )
    disp('Two Cameras acquisition!');
end
% end of Verbose

%%%%%%%%%%%%%%%%%%%%%
% Sub ROI
%%%%%%%%%%%%%%%%%%%%%
if( b_SubROI )
    ButtonName = questdlg('Would you like to use a pre-defined ROI?', ...
        'ROI', ...
        'Pre-defined', 'Draw', 'Cancel', 'Draw');
    switch ButtonName
        case 'Pre-defined'
             [filename, pathname] = uigetfile('*.mat', 'Select ROI file');
            if isequal(filename,0) || isequal(pathname,0)
                disp('User pressed cancel')
                Pos = [1 1 1023 1023];
            else
                load([pathname filesep filename]);
            end
        case 'Draw'
            dat = memmapfile([FolderName filesep...
                imgFilesList(1).name],...
                'Offset', hWima*4 + 5*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
            dat = dat.Data.imgj;
            fig = figure; imagesc(dat);
            h = imrect;
            wait(h);
            Pos = h.getPosition;
            close(fig);
        case 'Cancel'
            disp('User pressed cancel')
            Pos = [1 1 1023 1023];
    end
   
   LimX = [round(Pos(1)) round(Pos(1)+Pos(3))];
   LimY = [round(Pos(2)) round(Pos(2)+Pos(4))];
   save([FolderName filesep 'ROI.mat'],'Pos');
else
   LimX = [1 ImRes_XY(1)];
   LimY = [1 ImRes_XY(2)];
end
Rx = LimX(2) - LimX(1) + 1;
Ry = LimY(2) - LimY(1) + 1;

%%%%%%%%%%%%%%%%%%%%%
% Binning
%%%%%%%%%%%%%%%%%%%%%
if( BinningSpatial )
    %Verbose
    disp('Binning Spatial option is ON');
    %end of Verbose
    
    Rx = round(Rx/BinningSpatial);
    Ry = round(Ry/BinningSpatial);
else
    Rx = Rx;
    Ry = Ry;
end

%%%%%%%%%%%%%%%
% Analog inputs
%%%%%%%%%%%%%%%
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', hWai*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, AcqInfoStream.AISampleRate, tAIChan, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],tAIChan);
    AnalogIN = [AnalogIN; tmp];
end
clear tmp ind data;

%%%%%%%%%%%%%%%%%%%
%Stimulation Params
%%%%%%%%%%%%%%%%%%%
Str = [];
NbStim = 0;
StimLength = 0;
InterStim_min = 0;
InterStim_max = 0;
Stim = 0;
CamTrig = 0;
if( exist([FolderName filesep 'StimParameters.mat'], 'file') )
    load([FolderName filesep 'StimParameters.mat']);
    if( NbStim > 0 )
        Str = sprintf('%s\r%s%s\r%s%s%s',...
            'Stim detected: yes',...
            'Number of events: ', int2str(NbStim),...
            'Length of each event: ', int2str(StimLength), ' sec');
        bStim = 1;
    else
        Str = sprintf('%s',...
            'Stim detected: no');
        bStim = 0;
    end
    fprintf(Str);
    fprintf('\n');
    
else
    bStim = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Camera Trigs
%%%%%%%%%%%%%%%%%%%%%%%
CamTrig = find((AnalogIN(1:(end-1),1) < 2.5) & (AnalogIN(2:end,1) >= 2.5))+1;
StartDelay = round(CamTrig(1)/10);
EndDelay = round((length(AnalogIN(:,1)) - CamTrig(end))/10);

% Verbose
Str = sprintf('%s\r%s\r%s',...
    ['Camera Trigs detected: ' int2str(length(CamTrig))],...
    ['Recording of analog inputs starts ' int2str(StartDelay) ' ms before the first trigger.'],...
    ['Recording of analog inputs ends ' int2str(EndDelay) ' ms after the last trigger.']);

fprintf(Str);
fprintf('\n');
% end of Verbose
clear StartDelay EndDelay

%Less trig than images... something's wrong!
% if( length(CamTrig) < NombreImage  )
%     disp('IOI Error: Analog recordings and Image files don''t match. Impossible to continue further.');
%     out = 'Error';
%     return
% end

%%%%%%%%%%%%%%%%%%%%%%%
% Illumination Sequence
%%%%%%%%%%%%%%%%%%%%%%%
nbColors = sum(cellfun(@(X) contains(X,'Illumination'),fieldnames(AcqInfoStream)));
nbFluo = 0;
for indC = 1:nbColors
    Tag = eval(['AcqInfoStream.Illumination' int2str(indC) ';']);
    if( contains(Tag.Color, 'Fluo') )
        nbFluo = nbFluo + 1;
    end
end
bFluo = zeros(nbFluo,1); nFluo = zeros(nbFluo,1); tFluo = {}; iFluo = 1; cFluo = zeros(nbFluo,1);
bGreen = 0; nGreen = 0; cGreen = 0;
bRed = 0; nRed = 0; cRed = 0;
bYellow = 0; nYellow = 0; cYellow = 0;
bSpeckle = 0;nSpeckle = 0; cSpeckle = 0;
for indC = 1:nbColors
   Tag = eval(['AcqInfoStream.Illumination' int2str(indC) ';']);
   if( contains(Tag.Color, 'Fluo') )
       tFluo{iFluo} = Tag.Color;
       Tag.Color = 'Fluo';
   end
   switch(Tag.Color)
       case 'Fluo'
           bFluo(iFluo) = 1;
           if( AcqInfoStream.MultiCam == 0 )
               nFluo(iFluo) = indC;
               cFluo(iFluo) = 1;
           else
               cFluo(iFluo) = Tag.CamIdx;
               nFluo(iFluo) = Tag.FrameIdx;
           end
           iFluo = iFluo + 1;
       case 'Green'
           bGreen = 1;
           if( AcqInfoStream.MultiCam == 0 )
               cGreen = 1;
               nGreen = indC;
           else
               cGreen = Tag.CamIdx;
               nGreen = Tag.FrameIdx;
           end
           
       case 'Amber'
           bYellow = 1;
           if( AcqInfoStream.MultiCam == 0 )
               cYellow = 1;
               nYellow = indC;
           else
               cYellow = Tag.CamIdx;
               nYellow = Tag.FrameIdx;
           end
       case 'Red'
           bRed = 1;
           if( AcqInfoStream.MultiCam == 0 )
               cRed = 1;
               nRed = indC;
           else
               cRed = Tag.CamIdx;
               nRed = Tag.FrameIdx;
           end
       case 'Speckle'
           bSpeckle = 1;
           if( AcqInfoStream.MultiCam == 0 )
               cSpeckle = 1;
               nSpeckle = indC;
           else
               cSpeckle = Tag.CamIdx;
               nSpeckle = Tag.FrameIdx;
           end
   end
end

Freq = AcqInfoStream.FrameRateHz;
if( AcqInfoStream.MultiCam )
    nbCam = 2;
else
    nbCam = 1;
end
if( any(bFluo) )
    Flist = dir([FolderName filesep 'Data_Fluo*.mat']);
    if( ~isempty(Flist) )
        for ind = 1:size(Flist,1)
            delete([FolderName filesep Flist(ind).name]);
        end
    end
    if( nbFluo > 1 )
        for ind = 1:nbFluo
            waveTag = regexp(tFluo{ind}, '[0-9]{3}','match');
            
            fFluo{ind} = matfile([FolderName filesep 'Data_Fluo_' waveTag{:} '.mat'],'Writable',true);
            fFluo{ind}.datFile = [FolderName filesep 'fChan_' waveTag{:} '.dat'];
            fFluo{ind}.datSize = [Rx, Ry];
            fFluo{ind}.Stim = zeros(floor(nbCam*NombreImage/(nbColors*BinningTemp)),1, 'single');
            fFluo{ind}.Freq = nbCam*Freq/(nbColors*BinningTemp);
            fFluo{ind}.Wavelength = waveTag;
            fidF(ind) = fopen([FolderName filesep 'fChan_' waveTag{:} '.dat'],'w');
        end
    else
        fFluo = matfile([FolderName filesep 'Data_Fluo.mat'],'Writable',true);
        fFluo.datFile = [FolderName filesep 'fChan.dat'];
        fFluo.datSize = [Rx, Ry];
        fFluo.Stim = zeros(floor(nbCam*NombreImage/(nbColors*BinningTemp)),1, 'single');
        fFluo.Freq = nbCam*Freq/(nbColors*BinningTemp);
        fidF = fopen([FolderName filesep 'fChan.dat'],'w');
        fFluo.Wavelength = 1;
    end
end
if( bSpeckle )
    if( exist([FolderName filesep 'Data_speckle.mat'],'file') )
        delete([FolderName filesep 'Data_speckle.mat']);
    end
    fSpeckle = matfile([FolderName filesep 'Data_speckle.mat'],'Writable',true);
    fSpeckle.datFile = [FolderName filesep 'sChan.dat'];
    fSpeckle.datSize = [Rx, Ry];
    fSpeckle.Stim = zeros(floor(nbCam*NombreImage/(nbColors*BinningTemp)),1, 'single');
    fSpeckle.Freq = nbCam*Freq/(nbColors*BinningTemp);
    fidS = fopen([FolderName filesep 'sChan.dat'],'w');
end
if( bRed )
    if( exist([FolderName filesep 'Data_red.mat'],'file') )
        delete([FolderName filesep 'Data_red.mat']);
    end
    fRed = matfile([FolderName filesep 'Data_red.mat'],'Writable',true);
    fRed.datFile = [FolderName filesep 'rChan.dat'];
    fRed.datSize = [Rx, Ry];
    fRed.Stim = zeros(floor(nbCam*NombreImage/(nbColors*BinningTemp)), 1, 'single');
    fRed.Freq = nbCam*Freq/(nbColors*BinningTemp);
    fidR = fopen([FolderName filesep 'rChan.dat'],'w');
end
if( bYellow )
    if( exist([FolderName filesep 'Data_yellow.mat'],'file') )
        delete([FolderName filesep 'Data_yellow.mat']);
    end
    fYellow = matfile([FolderName filesep 'Data_yellow.mat'],'Writable',true);
    fYellow.datFile = [FolderName filesep 'yChan.dat'];
    fYellow.datSize = [Rx, Ry];
    fYellow.Stim = zeros(floor(nbCam*NombreImage/(nbColors*BinningTemp)),1, 'single');
    fYellow.Freq = nbCam*Freq/(nbColors*BinningTemp);
    fidY = fopen([FolderName filesep 'yChan.dat'],'w');
end
if( bGreen )
    if( exist([FolderName filesep 'Data_green.mat'],'file') )
        delete([FolderName filesep 'Data_green.mat']);
    end
    fGreen = matfile([FolderName filesep 'Data_green.mat'],'Writable',true);
    fGreen.datFile = [FolderName filesep 'gChan.dat'];
    fGreen.datSize = [Rx, Ry];
    fGreen.Stim = zeros(floor(nbCam*NombreImage/(nbColors*BinningTemp)), 1, 'single');
    fGreen.Freq = nbCam*Freq/(nbColors*BinningTemp);
    fidG = fopen([FolderName filesep 'gChan.dat'],'w');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolation for bad or missing frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImAddBook1 = ImagesCheck(AcqInfoStream.Camera_Model, idImg, nbColors/nbCam, imgFilesList, 1, NombreImage );
if( AcqInfoStream.MultiCam )
    ImAddBook2 = ImagesCheck(AcqInfoStream.Camera_Model,idImg2, nbColors/nbCam, imgFilesList2, 2, NombreImage2 );   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Images Classification, filtering and writing on disk:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImsMin = floor(NombreImage/(nbColors));
ImsMin = ImsMin*nbColors;
if( bFluo )
    if( nbFluo > 1 )
        for ind = 1:nbFluo
            waveTag = fFluo{ind}.Wavelength;
            disp(['Fluorescence ' waveTag{:} 'nm channel classification:']);
            tags = nFluo:nbColors/nbCam:ImsMin;
            if( cFluo(ind) == 1 )
                WriteChannel(tags, fFluo{ind}, 1, imgFilesList, fidF(ind), ImAddBook1, 1);
            else
                WriteChannel(tags, fFluo{ind}, 1, imgFilesList2, fidF(ind), ImAddBook2, 2);
            end
            disp('done');
        end
    else
        waveTag = fFluo.Wavelength;
        disp(['Fluorescence channel classification:']);
        tags = nFluo:nbColors/nbCam:ImsMin;
        WriteChannel(tags, fFluo, 1, imgFilesList, fidF, ImAddBook1, 1);
        disp('done');
    end
    
end
if( bSpeckle )
    disp('Speckle channel classification:');
    tags = nSpeckle:nbColors/nbCam:ImsMin;
    if( cSpeckle == 1)
        WriteChannel(tags, fSpeckle, 1, imgFilesList, fidS, ImAddBook1, 1);
    else
        WriteChannel(tags, fSpeckle, 1, imgFilesList2, fidS, ImAddBook2, 2);
    end
    disp('done');
end
if( bRed )
    disp('Red channel classification:');
    tags = nRed:nbColors/nbCam:ImsMin;
    if( cRed == 1)
        WriteChannel(tags, fRed, 1, imgFilesList, fidR, ImAddBook1, 1);
    else
        WriteChannel(tags, fRed, 1, imgFilesList2, fidR, ImAddBook2, 2);
    end
    disp('Done');
end
if( bYellow )
    disp('Yellow channel classification:');
    tags = nYellow:nbColors/nbCam:ImsMin;
    if( cYellow == 1)
        WriteChannel(tags, fYellow, 1, imgFilesList, fidY, ImAddBook1, 1);
    else
        WriteChannel(tags, fYellow, 1, imgFilesList2, fidY, ImAddBook2, 2);
    end
    disp('Done.');
end
if( bGreen )
    disp('Green channel classification:');
    tags = nGreen:nbColors/nbCam:ImsMin;
    if( cGreen == 1)
        WriteChannel(tags, fGreen, 1, imgFilesList, fidG, ImAddBook1, 1);
    else
        WriteChannel(tags, fGreen, 1, imgFilesList2, fidG, ImAddBook2, 2);
    end
    disp('done');
end
fprintf('\n');
%Verbose
fprintf(['Done!']);
fprintf('\n');
%end of Verbose

    function WriteChannel(Tags, fPtr, cPtr, iFList, fidPtr, imAddr, iCam)
        
        PrcTag = round(linspace(0, length(Tags), 20));
        
        indT = 1;
        if( rem(length(Tags),BinningTemp) > 0 )
            Tags = Tags(1:(end - rem(length(Tags),BinningTemp)));
        end
        for indI = 1:BinningTemp:length(Tags)
            Images = zeros(ImRes_XY(1), ImRes_XY(2), BinningTemp, 'single');
            for indB = 0:(BinningTemp-1)
                indF = Tags(indI + indB);
                
                if( imAddr(indF,1) <= size(iFList,1) )
                    datloc =   memmapfile([FolderName filesep...
                        iFList(imAddr(indF,1)).name],...
                        'Offset', hWima*4 + (imAddr(indF,2)-1)*SizeImage,...
                        'Format', frameFormat, 'repeat', 1);
                else
                    datloc =   memmapfile([FolderName filesep 'img_interp_' int2str(iCam) '.bin'],...
                        'Offset', (imAddr(indF,2)-1)*SizeImage,...
                        'Format', frameFormat, 'repeat', 1);
                end
                Images(:,:,indB+1) = single(datloc.Data.imgj);
                if( bStim )
                    fPtr.Stim(cPtr,1) = single(Stim(indF));
                else
                    fPtr.Stim(cPtr,1) = fPtr.Stim(cPtr,1);
                end
            end
            img = mean(Images,3);
            
            if( b_SubROI )
                img = img(round(LimY(1)):round(LimY(2)),round(LimX(1)):round(LimX(2)));
            end
            
            if( BinningSpatial )
                img = imresize(img,1/BinningSpatial);
            end
            fwrite(fidPtr, img, 'single');
            cPtr = cPtr + 1;
            
            if( indI >= PrcTag(indT) )
                P = round((100*PrcTag(indT))/length(Tags));
                
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                indT = indT + 1;
            end
        end
        
        fPtr.datLength = cPtr - 1;
        fPtr.FirstDim = 'y';
        fclose(fidPtr);
    end

    function ImAddBook = ImagesCheck(cModel, idI, nColors, ifList, camID, nIms )
        badFrames = [];
        if( strcmp(cModel, 'CS2100M') )
            if( idI(1,1) > 1 )
                idI(1,1) = 0;
                idI(1,2) = 0;
            end
            SkipNFirst = sum(idI(:,1) == 0);
            MissingOffset = cumsum(idI(:,2));
            idI(:,1) = idI(:,1) + MissingOffset;
            goodFrames = find(accumarray(idI((SkipNFirst+1):end,1),1) >= 1)';
            badFrames = 1:max(goodFrames(:));
            badFrames = badFrames(~ismember(badFrames, goodFrames));
        elseif( strcmp(cModel, 'D1024') ||...
                strcmp(cModel, 'D1312'))
            if( idI(1,1) > 1 )
                idI(1,1) = 0;
                idI(1,2) = 0;
            end
            SkipNFirst = sum(idI(:,1) == 0);
            MissingOffset = cumsum(idI(:,2));
            idI(:,1) = idI(:,1) + MissingOffset;
            goodFrames = find(accumarray(idI((SkipNFirst+1):end,1),1)==1)';
            ConseqFromLeft = [1 diff(goodFrames,1,2)==1];
            ConseqFromRight = fliplr([true diff(fliplr(goodFrames),1,2)==-1]);
            goodFrames = goodFrames(ConseqFromLeft|ConseqFromRight);
            badFrames = 1:max(goodFrames(:));
            badFrames = badFrames(~ismember(badFrames, goodFrames));
        end
        
        %%% Lookup Table For missing frames
        InterpLUT = zeros(8,1);
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
                tmpBefore = tmpID - (nColors:nColors:(tmpID-1));
                idx = find(ismember(tmpBefore,idI(:,1))&~ismember(tmpBefore,badFrames),1,'first');
                tmpBefore = tmpBefore(idx);
                if( ~isempty(tmpBefore) )
                    InterpLUT(1,ind) = tmpBefore;
                    idx = find(tmpBefore == idI(:,1),1,'first');
                    InterpLUT(2,ind) = floor((idx-1)/NbImsPefFile) + 1;
                    InterpLUT(3,ind) = rem((idx-1),NbImsPefFile) + 1;
                end
                
                tmpAfter = tmpID + (nColors:nColors:(nIms));
                idx = find(ismember(tmpAfter,idI(:,1))&~ismember(tmpAfter,badFrames),1,'first');
                tmpAfter = tmpAfter(idx);
                if( isempty(tmpBefore) )
                    InterpLUT(1,ind) = tmpAfter;
                    idx = find(tmpAfter == idI(:,1),1,'first');
                    InterpLUT(2,ind) = floor((idx-1)/NbImsPefFile) + 1;
                    InterpLUT(3,ind) = rem((idx-1),NbImsPefFile) + 1;
                end
                InterpLUT(4,ind) = tmpAfter;
                idx = find(tmpAfter == idI(:,1), 1, 'first');
                InterpLUT(5,ind) = floor((idx-1)/NbImsPefFile) + 1;
                InterpLUT(6,ind) = rem((idx-1),NbImsPefFile) + 1;
                
                tmpRatio =  (tmpID - InterpLUT(1,ind))./...
                    (InterpLUT(4,ind) - InterpLUT(1,ind));
                InterpLUT(7,ind) = tmpRatio;
                InterpLUT(8,ind) = badFrames(ind);
            end
            clear tmpRatio tmpAfter tmpBefore;
            %%% Interpolation of missing frames
            TmpFrames = struct('framej',[], 'imgj',[]);
            for ind = 1:size(InterpLUT,2)
                dBefore = memmapfile([FolderName filesep...
                    ifList(InterpLUT(2,ind)).name],...
                    'Offset', hWima*4 + (InterpLUT(3,ind) - 1)*SizeImage,...
                    'Format', frameFormat, 'repeat', 1);
                dAfter = memmapfile([FolderName filesep...
                    ifList(InterpLUT(5,ind)).name],...
                    'Offset', hWima*4 + (InterpLUT(6,ind) - 1)*SizeImage,...
                    'Format', frameFormat, 'repeat', 1);
                TmpFrames(ind).imgj = uint16(round(InterpLUT(7,ind)*double(dAfter.Data.imgj - dBefore.Data.imgj))) + dBefore.Data.imgj;
                TmpFrames(ind).framej = uint64([InterpLUT(8,ind), 1, 1]);
            end
            
            fid = fopen([FolderName filesep 'img_interp_' int2str(camID) '.bin'],'w');
            for ind = 1:size(InterpLUT,2)
                fwrite(fid, TmpFrames(ind).framej, 'uint64');
                fwrite(fid, TmpFrames(ind).imgj, 'uint16');
            end
            fclose(fid);
        end
        
        %Rebuilding addresses for each frames...
        NbImage = max(goodFrames);
        ImAddressBook = zeros(NbImage,2);
        for ind = 1:NbImage
            if( ismember(ind, badFrames) )
                fidx = find( ind == InterpLUT(8,:), 1, 'first');
                ImAddressBook(ind,1) = size(ifList,1) + 1;
                ImAddressBook(ind,2) = fidx;
            elseif( ismember(ind, goodFrames) )
                fidx = find( ind == idI, 1, 'first');
                ImAddressBook(ind,1) = floor((fidx-1)/NbImsPefFile) + 1;
                ImAddressBook(ind,2) = rem(fidx-1, NbImsPefFile) + 1;
            end
        end
        
        %Saving infos
        if( ~strcmp(FolderName(end), filesep) )
            FolderName = [FolderName filesep];
        end
        save([FolderName 'ImagesLUT_' int2str(camID) '.mat'], 'ImAddressBook');
        ImAddBook = ImAddressBook;
    end
end
