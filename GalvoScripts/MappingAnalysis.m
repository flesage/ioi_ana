function Data = MappingAnalysis(FolderPath)
%Infos files opening
if( ~strcmp(FolderPath(end), filesep) )
    FolderPath = strcat(FolderPath, filesep);
end
InfoStim = Read_OptoGenParams_File(FolderPath);
%InfoAcq = ReadInfoFile(FolderPath);
InfoAcq.AISampleRate = 10000;
InfoAcq.AINChannels = 11;

%Analog Inputs reading:
aiFiles = dir([FolderPath 'ai_*.bin']);
AnalogIn = [];
for indF = 1:length(aiFiles)
    fid = fopen([FolderPath aiFiles(indF).name]);
    fseek(fid,5*4,'bof');
    dat = fread(fid,inf,'double');
    dat = reshape(dat, InfoAcq.AISampleRate, InfoAcq.AINChannels, []);
    dat = permute(dat,[2 1 3]);
    dat = reshape(dat,11,[]);
    AnalogIn = [AnalogIn, dat];
end
clear dat indF aiFiles fid

%Stimulation detection: 
Stim = AnalogIn(2,:);
iStim = find(diff(Stim,1,2)>2.50);
Xpos = InfoStim.Positions(:,1);
Ypos = InfoStim.Positions(:,2);
sX = length(unique(Xpos));
sY = length(unique(Ypos));
iX = round((sX-1)*(Xpos - min(Xpos(:)))/(max(Xpos(:)) - min(Xpos(:))) + 1);
iY = round((sY-1)*(Ypos - min(Ypos(:)))/(max(Ypos(:)) - min(Ypos(:))) + 1);
Xpos = unique(Xpos);
Ypos = unique(Ypos);

%Analog Inputs reshaping based on stimulations:
%Channels for Paw 1: AI #1(4 in Matlab), AI #2(5 in Matlab), AI #3(6 in Matlab)
%Channels for Paw 2: AI #5(8 in Matlab), AI 6#(9 in Matlab), AI #7(10 in Matlab)
f = fdesign.lowpass('N,F3dB', 4, 250, 10000);
lpass = design(f,'butter');
clear f
AnalogIn = filtfilt(lpass.sosMatrix, lpass.ScaleValues, AnalogIn');
clear lpass

Mp1 = zeros(sY, sX, 3, 4000);
Mp2 = zeros(sY, sX, 3, 4000);
for ind = 1:size(iStim,2)
    idx = iStim(ind) + (-500:3499);
    Mp1(iY(ind), iX(ind), 1, :) = AnalogIn(idx,4);
    Mp1(iY(ind), iX(ind), 2, :) = AnalogIn(idx,5);
    Mp1(iY(ind), iX(ind), 3, :) = AnalogIn(idx,6);
    
    Mp2(iY(ind), iX(ind), 1, :) = AnalogIn(idx,8);
    Mp2(iY(ind), iX(ind), 2, :) = AnalogIn(idx,9);
    Mp2(iY(ind), iX(ind), 3, :) = AnalogIn(idx,10);
end
clear aIn dat idx ind iStim iX iY Stim sX sY
dims = size(Mp1);

Mp1 = bsxfun(@minus, Mp1, mean(Mp1,4));
Mp2 = bsxfun(@minus, Mp2, mean(Mp2,4));

Vp1 = squeeze(sqrt(sum((Mp1).^2,3)/3));
Vp2 = squeeze(sqrt(sum((Mp2).^2,3)/3));

mVp1 = mean(reshape(Vp1(:,:,1:500),[],1));
sVp1 = std(reshape(Vp1(:,:,1:500),[],1),0,1);
zVp1 = bsxfun(@rdivide, bsxfun(@minus, Vp1, mVp1), sVp1);
%zVp1 = imfilter(zVp1, fspecial('gaussian',3,1.5),'same', 'symmetric');
mVp2 = mean(Vp2(:,:,1:500),3);
sVp2 = std(Vp2(:,:,1:500),0,3);
zVp2 = bsxfun(@rdivide, bsxfun(@minus, Vp2, mVp2), sVp2);
%zVp2 = imfilter(zVp2, fspecial('gaussian',3,1.5),'same', 'symmetric');
clear mVp* sVp*

[vM1, tM1] = max(zVp1,[],3);
MapActiv1 = vM1 >= 12;
MapActiv1(tM1<500) = false;
MapActiv1(tM1>1000) = false;
tM1 = tM1 - 500;
tM1(~MapActiv1) = 0;
tM1 = tM1/10;

[vM2, tM2] = max(zVp2,[],3);
MapActiv2 = vM2 >= 12;
MapActiv2(tM2<500) = false;
MapActiv2(tM2>1000) = false;
tM2 = tM2 - 500;
tM2(~MapActiv2) = 0;
tM2 = tM2/10;

MapActiv = MapActiv1 | MapActiv2;

tZ1 = zeros(dims(1), dims(2));
tZ2 = zeros(dims(1), dims(2));
for indX = 1:dims(2)
    for indY = 1:dims(1)
        tmp1 = find(zVp1(indY, indX,:)>=3, 1,'first');
        tmp2 = find(zVp2(indY, indX,:)>=3, 1,'first');
        if( isempty(tmp1) )
            tZ1(indY, indX) = 0;
        else
            tZ1(indY, indX) = tmp1;
        end
        if( isempty(tmp2) )
            tZ2(indY, indX) = 0;
        else
            tZ2(indY, indX) = tmp2;
        end
    end
end
tZ1 = tZ1 - 500;
tZ1(tZ1 < 0) = 0;
tZ1(tZ1 > 500) = 0;
tZ1 = tZ1/10;
tZ1(~MapActiv) = 0;

tZ2 = tZ2 - 500;
tZ2(tZ2 < 0) = 0;
tZ2(tZ2 > 500) = 0;
tZ2 = tZ2/10;
tZ2(~MapActiv) = 0;

clear tmp* ind*

Data.Vp123 = Vp1;
Data.Vp567 = Vp2;
Data.timeOnset_123 = tZ1;
Data.timeOnset_567 = tZ2;
Data.timeMax_123 = tM1;
Data.timeMax_567 = tM2;
Data.valueMax_123 = max(Vp1,[],3);
Data.valueMax_567 = max(Vp2,[],3);
Data.valueMin_123 = min(Vp1,[],3);
Data.valueMin_567 = min(Vp2,[],3);
Data.sagittal_axis = Xpos;
Data.coronal_axis = Ypos;
Data.Infos = InfoStim;

save([FolderPath 'DataMapping.mat'], 'Data');
%Generate Figures:
Xpix = round(InfoStim.RefX + Xpos/InfoStim.MMpPix);
Ypix = round(InfoStim.RefY + Ypos/InfoStim.MMpPix);
%Max: 
hfig = figure;
axRef = axes('Parent', hfig);
axOvr = axes('Parent', hfig);
axis(axOvr, 'off');

IList = dir([FolderPath '*.png']);
if( length(IList) > 1 )
    str = {IList.name};
    [v,c] = listdlg('PromptString','Select a file for Reference:',...
                'SelectionMode','single',...
                'ListString',str);
    if( c < 1 )
        IRef = zeros(1024,1024,3);
    else
        IList = IList(v);
        IRef = imread([FolderPath IList.name]);
    end
elseif( isempty(IList) )
    IRef = zeros(1024,1024,3);
else
    IRef = imread([FolderPath IList.name]);
end

Xaxmm = ((1:1024) - InfoStim.RefX)*InfoStim.MMpPix;
Yaxmm = ((1:1024) - InfoStim.RefY)*InfoStim.MMpPix;
imagesc(axRef, Xaxmm, Yaxmm, IRef);
hold(axRef,'on');
plot(axRef, 0, 0, 'or');
text(axRef, 0.1, -0.1, '{\beta}','FontSize',16,'FontWeight', 'bold', 'Color', 'r')

%MAX 3-4-5
imagesc(axOvr, Ypos, Xpos, flipud(rot90(vM1)),'AlphaData', flipud(rot90(MapActiv))*0.5);
title('Max Amplitude Channels 3-4-5')
xlabel(axRef, 'Coronal axis (mm)');
ylabel(axRef, 'Sagital axis (mm)');
axis(axRef, 'image');
axis(axOvr, 'image');
axis(axOvr,'off');
colorbar('AxisLocation','in');
linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
    'CameraPosition', 'XLim', 'YLim'});
saveas(hfig, [FolderPath 'MaxAmpChan345.png']);

%MAX 7-8-9
imagesc(axOvr, Ypos, Xpos, flipud(rot90(vM2)),'AlphaData', flipud(rot90(MapActiv))*0.5);
title('Max Amplitude Channels 7-8-9')
xlabel(axRef, 'Coronal axis (mm)');
ylabel(axRef, 'Sagital axis (mm)');
axis(axRef, 'image');
axis(axOvr, 'image');
axis(axOvr,'off');
colorbar('AxisLocation','in');
linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
    'CameraPosition', 'XLim', 'YLim'});
saveas(hfig, [FolderPath 'MaxAmpChan789.png']);

%Tmax 3-4-5
imagesc(axOvr, Ypos, Xpos, flipud(rot90(tM1)),'AlphaData', flipud(rot90(MapActiv))*0.5);
title('Rising Time to Maximum Channels 3-4-5')
xlabel(axRef, 'Coronal axis (mm)');
ylabel(axRef, 'Sagital axis (mm)');
axis(axRef, 'image');
axis(axOvr, 'image');
axis(axOvr,'off');
colorbar('AxisLocation','in');
linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
    'CameraPosition', 'XLim', 'YLim'});
saveas(hfig, [FolderPath 'TmaxChan345.png']);

%Tmax 7-8-9
imagesc(axOvr, Ypos, Xpos, flipud(rot90(tM2)),'AlphaData', flipud(rot90(MapActiv))*0.5);
title('Rising Time to Maximum Channels 7-8-9')
xlabel(axRef, 'Coronal axis (mm)');
ylabel(axRef, 'Sagital axis (mm)');
axis(axRef, 'image');
axis(axOvr, 'image');
axis(axOvr,'off');
colorbar('AxisLocation','in');
linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
    'CameraPosition', 'XLim', 'YLim'});
saveas(hfig, [FolderPath 'TmaxChan789.png']);

%TOnset 3-4-5
imagesc(axOvr, Ypos, Xpos, flipud(rot90(tZ1)),'AlphaData', flipud(rot90(MapActiv))*0.5);
title('Onset Time to Maximum Channels 3-4-5')
xlabel(axRef, 'Coronal axis (mm)');
ylabel(axRef, 'Sagital axis (mm)');
axis(axRef, 'image');
axis(axOvr, 'image');
axis(axOvr,'off');
caxis(axOvr,[0 50]);
colorbar('AxisLocation','in');
linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
    'CameraPosition', 'XLim', 'YLim'});
saveas(hfig, [FolderPath 'TonChan345.png']);

%TOnset 7-8-9
imagesc(axOvr, Ypos, Xpos, flipud(rot90(tZ1)),'AlphaData', flipud(rot90(MapActiv))*0.5);
title('Onset Time to Maximum Channels 7-8-9')
xlabel(axRef, 'Coronal axis (mm)');
ylabel(axRef, 'Sagital axis (mm)');
axis(axRef, 'image');
axis(axOvr, 'image');
axis(axOvr,'off');
caxis(axOvr,[0 50]);
colorbar('AxisLocation','in');
linkprop([axRef axOvr],{'Position', 'Units','OuterPosition'...
    'CameraPosition', 'XLim', 'YLim'});
saveas(hfig, [FolderPath 'TonChan789.png']);

close(hfig);
end