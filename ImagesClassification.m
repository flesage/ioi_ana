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
if( ~b_IgnoreStim && (AcqInfoStream.Stimulation > 0) )
    ReadAnalogsIn(FolderName, AcqInfoStream);
else
    fprintf('No stimulation detected. \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Images sequence validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgFilesList = dir([FolderName filesep 'img_0*.bin']);
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
clear nx ny data header ind;

% Verbose
disp(['Opening of: ' FolderName]);
disp(['Number of Frames acquired: ' int2str(NombreImage)]);
disp(['Frames'' resolution: ' int2str(ImRes_XY(1)) ' pix X ' int2str(ImRes_XY(2)) ' pix']);
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

%%%%%%%%%%%%%%%%%%%%%
% Analog inputs
%%%%%%%%%%%%%%%%%%%%%
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

%%%%
%Stimulation Params
%%%%
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
if( length(CamTrig) < NombreImage  )
    disp('IOI Error: Analog recordings and Image files don''t match. Impossible to continue further.');
%     out = 'Error';
%     return
end

%%%%%%%%%%%%%%%%%%%%%%%
% Illumination Sequence
%%%%%%%%%%%%%%%%%%%%%%%
nbColors = sum(cellfun(@(X) contains(X,'Illumination'),fieldnames(AcqInfoStream)));

bFluo = 0; nFluo = 0;
bGreen = 0; nGreen = 0;
bRed = 0; nRed = 0;
bYellow = 0; nYellow = 0;
bSpeckle = 0;nSpeckle = 0;
for indC = 1:nbColors
   Tag = eval(['AcqInfoStream.Illumination' int2str(indC) ';']);
   if( contains(Tag, 'Fluo') )
       Tag = 'Fluo';
   end
   switch(Tag)
       case 'Fluo'
           bFluo = 1;
           nFluo = indC;
       case 'Green'
           bGreen = 1;
           nGreen = indC;
       case 'Amber'
           bYellow = 1;
           nYellow = indC;
       case 'Red'
           bRed = 1;
           nRed = indC;
       case 'Speckle'
           bSpeckle = 1;
           nSpeckle = indC;
   end
end

Freq = AcqInfoStream.FrameRateHz;

if( bFluo )
       if( exist([FolderName filesep 'Data_Fluo.mat'],'file') )
            delete([FolderName filesep 'Data_Fluo.mat']);
        end
        fFluo = matfile([FolderName filesep 'Data_Fluo.mat'],'Writable',true);
        fFluo.datFile = [FolderName filesep 'fChan.dat'];
        fFluo.datSize = [Rx, Ry];
        fFluo.Stim = zeros(floor(NombreImage/(nbColors*BinningTemp)),1, 'single');
        fFluo.Freq = Freq/(nbColors*BinningTemp);
        cFluo = 1;
        fidF = fopen([FolderName filesep 'fChan.dat'],'w');
end    
if( bSpeckle )
        if( exist([FolderName filesep 'Data_speckle.mat'],'file') )
            delete([FolderName filesep 'Data_speckle.mat']);
        end
        fSpeckle = matfile([FolderName filesep 'Data_speckle.mat'],'Writable',true);
        fSpeckle.datFile = [FolderName filesep 'sChan.dat'];
        fSpeckle.datSize = [Rx, Ry];
        fSpeckle.Stim = zeros(floor(NombreImage/(nbColors*BinningTemp)),1, 'single');
        fSpeckle.Freq = Freq/(nbColors*BinningTemp);
        cSpeckle = 1;
        fidS = fopen([FolderName filesep 'sChan.dat'],'w');

end
if( bRed )
    if( exist([FolderName filesep 'Data_red.mat'],'file') )
        delete([FolderName filesep 'Data_red.mat']);
    end
    fRed = matfile([FolderName filesep 'Data_red.mat'],'Writable',true);
    fRed.datFile = [FolderName filesep 'rChan.dat'];
    fRed.datSize = [Rx, Ry];
    fRed.Stim = zeros(floor(NombreImage/(nbColors*BinningTemp)), 1, 'single');
    fRed.Freq = Freq/(nbColors*BinningTemp);
    cRed = 1;
    fidR = fopen([FolderName filesep 'rChan.dat'],'w');
end
if( bYellow )
    if( exist([FolderName filesep 'Data_yellow.mat'],'file') )
        delete([FolderName filesep 'Data_yellow.mat']);
    end
    fYellow = matfile([FolderName filesep 'Data_yellow.mat'],'Writable',true);
    fYellow.datFile = [FolderName filesep 'yChan.dat'];
    fYellow.datSize = [Rx, Ry];
    fYellow.Stim = zeros(floor(NombreImage/(nbColors*BinningTemp)),1, 'single');
    fYellow.Freq = Freq/(nbColors*BinningTemp);
    cYellow = 1;
    fidY = fopen([FolderName filesep 'yChan.dat'],'w');
end
if( bGreen )    
    if( exist([FolderName filesep 'Data_green.mat'],'file') )
        delete([FolderName filesep 'Data_green.mat']);
    end
    fGreen = matfile([FolderName filesep 'Data_green.mat'],'Writable',true);
    fGreen.datFile = [FolderName filesep 'gChan.dat'];
    fGreen.datSize = [Rx, Ry];
    fGreen.Stim = zeros(floor(NombreImage/(nbColors*BinningTemp)), 1, 'single');
    fGreen.Freq = Freq/(nbColors*BinningTemp);
    cGreen = 1;
    fidG = fopen([FolderName filesep 'gChan.dat'],'w');
end

%Interpolation for bad or missing frames
if( strcmp(AcqInfoStream.Camera_Model, 'CS2100M') )
    if( idImg(1,1) > 1 )
        idImg(1,1) = 0;
        idImg(1,2) = 0;
    end
    SkipNFirst = sum(idImg(:,1) == 0);
    MissingOffset = cumsum(idImg(:,2));
    idImg(:,1) = idImg(:,1) + MissingOffset;
    goodFrames = find(accumarray(idImg((SkipNFirst+1):end,1),1) >= 1)';
    badFrames = 1:max(goodFrames(:));
    badFrames = badFrames(~ismember(badFrames, goodFrames));
elseif( strcmp(AcqInfoStream.Camera_Model, 'D1024') ||...
        strcmp(AcqInfoStream.Camera_Model, 'D1312'))
    if( idImg(1,1) > 1 )
        idImg(1,1) = 0;
        idImg(1,2) = 0;
    end
    SkipNFirst = sum(idImg(:,1) == 0);
    MissingOffset = cumsum(idImg(:,2));
    idImg(:,1) = idImg(:,1) + MissingOffset;
    goodFrames = find(accumarray(idImg((SkipNFirst+1):end,1),1)==1)';
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
        tmpBefore = tmpID - (nbColors:nbColors:(tmpID-1));
        idx = find(ismember(tmpBefore,idImg(:,1))&~ismember(tmpBefore,badFrames),1,'first');
        tmpBefore = tmpBefore(idx);
        if( ~isempty(tmpBefore) )
            InterpLUT(1,ind) = tmpBefore;
            idx = find(tmpBefore == idImg(:,1),1,'first');
            InterpLUT(2,ind) = floor((idx-1)/NbImsPefFile) + 1;
            InterpLUT(3,ind) = rem((idx-1),NbImsPefFile) + 1;
        end
        
        tmpAfter = tmpID + (nbColors:nbColors:(NombreImage));
        idx = find(ismember(tmpAfter,idImg(:,1))&~ismember(tmpAfter,badFrames),1,'first');
        tmpAfter = tmpAfter(idx);
        if( isempty(tmpBefore) )
            InterpLUT(1,ind) = tmpAfter;
            idx = find(tmpAfter == idImg(:,1),1,'first');
            InterpLUT(2,ind) = floor((idx-1)/NbImsPefFile) + 1;
            InterpLUT(3,ind) = rem((idx-1),NbImsPefFile) + 1;
        end
        InterpLUT(4,ind) = tmpAfter;
        idx = find(tmpAfter == idImg(:,1), 1, 'first');
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
NombreImage = max(goodFrames);
ImAddressBook = zeros(NombreImage,2);
for ind = 1:NombreImage
    if( ismember(ind, badFrames) )
        fidx = find( ind == InterpLUT(8,:), 1, 'first');
        ImAddressBook(ind,1) = size(imgFilesList,1) + 1;
        ImAddressBook(ind,2) = fidx;
    elseif( ismember(ind, goodFrames) )
        fidx = find( ind == idImg, 1, 'first');
        ImAddressBook(ind,1) = floor((fidx-1)/NbImsPefFile) + 1;
        ImAddressBook(ind,2) = rem(fidx-1, NbImsPefFile) + 1;
    end
end

%Saving infos...
if( ~strcmp(FolderName(end), filesep) )
    FolderName = [FolderName filesep];
end
save([FolderName 'ImagesLUT.mat'], 'ImAddressBook');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Images Classification, filtering and writing on disk:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( bFluo )
    disp('Fluorescence channel classification:');
    tags = nFluo:nbColors:NombreImage;
    WriteChannel(tags, fFluo, cFluo, fidF);
    disp('done');
end
if( bSpeckle )
    disp('Speckle channel classification:');
    tags = nSpeckle:nbColors:NombreImage;
    WriteChannel(tags, fSpeckle, cSpeckle, fidS);
    disp('done');
end
if( bRed )
    disp('Red channel classification:');
    tags = nRed:nbColors:NombreImage;
    WriteChannel(tags, fRed, cRed, fidR);
    disp('Done');
end
if( bYellow )
    disp('Yellow channel classification:');
    tags = nYellow:nbColors:NombreImage;
    WriteChannel(tags, fYellow, cYellow, fidY);
    disp('Done.');
end
if( bGreen )
    disp('Green channel classification:');
    tags = nGreen:nbColors:NombreImage;
    WriteChannel(tags, fGreen, cGreen, fidG);
    disp('done');
end
fprintf('\n');
%Verbose
fprintf(['Done!']);
fprintf('\n');
%end of Verbose

    function WriteChannel(Tags, fPtr, cPtr, fidPtr)
                
        PrcTag = round(linspace(0, length(Tags), 20));
        
        indT = 1;
        if( rem(length(Tags),BinningTemp) > 0 )
            Tags = Tags(1:(end - rem(length(Tags),BinningTemp)));
        end
        for indI = 1:BinningTemp:length(Tags)
            Images = zeros(ImRes_XY(1), ImRes_XY(2), BinningTemp, 'single');
            for indB = 0:(BinningTemp-1)
                indF = Tags(indI + indB);
                
                if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
                    datloc =   memmapfile([FolderName filesep...
                        imgFilesList(ImAddressBook(indF,1)).name],...
                        'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                        'Format', frameFormat, 'repeat', 1);
                else
                    datloc =   memmapfile([FolderName filesep 'img_interp.bin'],...
                        'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
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

end
