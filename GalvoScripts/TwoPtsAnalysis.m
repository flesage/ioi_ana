function TwoPtsAnalysis(FolderPath)

disp(['Two Points Analysis of: ' FolderPath]);
if( ~strcmp(FolderPath(end), filesep) )
    FolderPath = strcat(FolderPath, filesep);
end

InfoStim = Read_OptoGenParams_File(FolderPath);
AISampleRate = 10000;
AINChannels = 11;

disp('Reading Analog In files.');
aiFiles = dir([FolderPath 'ai_*.bin']);
AnalogIn = [];
for indF = 1:length(aiFiles)
    fid = fopen([FolderPath aiFiles(indF).name]);
    fseek(fid,5*4,'bof');
    dat = fread(fid,inf,'double');
    dat = reshape(dat, AISampleRate, AINChannels, []);
    dat = permute(dat,[2 1 3]);
    dat = reshape(dat,11,[]);
    AnalogIn = [AnalogIn, dat];
end
clear dat indF aiFiles fid

disp('Events detection');
StimE = find((diff(AnalogIn(2,:),1,2) > 2.5) & (AnalogIn(2,2:(end))>0));
[~,index] = sort(InfoStim.EOrder);
StimOrder = StimE(index);

ChanIdxPaw1 = 4:6;
ChanIdxPaw2 = 8:10;
hf1 = figure; plot(AnalogIn(ChanIdxPaw1,:)'); saveas(hf1, [FolderPath 'RawDat_123.png']);
hf2 = figure; plot(AnalogIn(ChanIdxPaw2,:)');saveas(hf2, [FolderPath 'RawDat_123.png']);
close(hf1); close(hf2);

disp('Signal segmentation');
Sig = zeros(2, InfoStim.NbEvnts*InfoStim.NbReps, 4500,'single');
for indS = 1:size(StimE,2)
    Sig(1, indS, :) = sum(AnalogIn(ChanIdxPaw1, StimOrder(indS) + (0:4499)).^2,1);
    Sig(2, indS, :) = sum(AnalogIn(ChanIdxPaw2, StimOrder(indS) + (0:4499)).^2,1);
end
Sig = bsxfun(@rdivide, Sig, mean(Sig,3));
Sig = reshape(Sig, 2, InfoStim.NbReps, InfoStim.NbEvnts,[]);

mSig = squeeze(mean(Sig,2));
sSig = squeeze(std(Sig,0,2));
if( InfoStim.NbEvnts == 1 )
    mSig = permute(mSig, [1 3 2]);
    sSig = permute(sSig, [1 3 2]);
end
zSig = bsxfun(@rdivide,bsxfun(@minus, mSig, mean(mSig(:,:,1:500),3)),...
    std(mSig(:,:,1:500),0,3));

TZeros = InfoStim.TOffset;
timeVect = linspace(-TZeros,450-TZeros, 4500); %in msec 
Results = matfile([FolderPath 'Results.mat'], 'Writable', true);
Results.timeVect = timeVect';

%@Maximum time
disp('Timing Parameters measurement')
maxTiming = zeros(2, length(InfoStim.IStimDelay), length(InfoStim.PreCondPwr));
onsetTiming = zeros(2, length(InfoStim.IStimDelay), length(InfoStim.PreCondPwr));
for indP = 1:length(InfoStim.PreCondPwr)
    pTag = int2str(InfoStim.PreCondPwr(indP));
    for indI = 1:length(InfoStim.IStimDelay)
        iTag = int2str(InfoStim.IStimDelay(indI));
        Cnt = (indP-1)*length(InfoStim.IStimDelay) + indI;
        if( InfoStim.NullCond )
            Cnt = Cnt + 1;
        end
        timeCourse = squeeze(mSig(1,Cnt,:)) - 1;
        eval(['Results.Ch123_Pwr_' pTag '_Delay_' iTag ' = timeCourse;']);
       
        timeCourse = (timeCourse - mean(timeCourse(4001:4500)))./(std(timeCourse(4001:4500)));
        timeCourse = medfilt1(timeCourse,5);
        idx = find(abs(timeCourse) > 6,1,'first');
        if( isempty( idx ) )
            disp(['No signal detected for Channels 1-2-3 @ ' pTag ' mW with ' iTag ' ms of interstim.']);
            onsetTiming(1,indI,indP) = -1;
            maxTiming(1,indI,indP) = -1;
        else
            idx = idx - find(abs(timeCourse(idx:-1:1)) < 3,1,'first');
            onsetTiming(1,indI,indP) = (idx - TZeros*10)/10;
            [~,idx] = max(abs(timeCourse(TZeros*10 + (1:500))));
            maxTiming(1,indI,indP) = idx/10;
        end
        timeCourse = squeeze(mSig(2,Cnt,:)) - 1;
        eval(['Results.Ch567_Pwr_' pTag '_Delay_' iTag ' = timeCourse;']);
        
        
        timeCourse = (timeCourse - mean(timeCourse(4001:4500)))./(std(timeCourse(4001:4500)));
        timeCourse = medfilt1(timeCourse,5);
        idx = find(abs(timeCourse) > 6,1,'first');
        if( isempty( idx ) )
            disp(['No signal detected for Channels 5-6-7 @ ' pTag ' mW with ' iTag ' ms of interstim.']);
            onsetTiming(1,indI,indP) = -1;
            maxTiming(2,indI,indP) = -1;
        else
            idx = idx - find(abs(timeCourse(idx:-1:1)) < 3,1,'first');
            onsetTiming(2,indI,indP) = (idx - TZeros*10)/10;
            [~,idx] = max(abs(timeCourse(TZeros*10 + (1:500))));
            maxTiming(2,indI,indP) = idx/10;
        end
    end
end

if( InfoStim.NbEvnts > 1 )
    disp('Multiple conditions graph generation');
    for indP = 1:length(InfoStim.PreCondPwr)
        figure;
        plot(InfoStim.IStimDelay, onsetTiming(:,:,indP));
        ylim([0 50]);
        title(['Onset time (msec) with pre-cond power @ ' int2str(InfoStim.PreCondPwr(indP)) '%']);
        legend('Channels 1-2-3', 'Channels 6-7-8');
        ylabel('Response Delay (ms)')
        xlabel('ISI (ms)');
        savefig(gcf, [FolderPath 'Onset_Pwr_' int2str(InfoStim.PreCondPwr(indP)) '.fig']);
        saveas(gcf, [FolderPath 'Onset_Pwr_' int2str(InfoStim.PreCondPwr(indP)) '.png']);
        
        close(gcf);
        
        figure;
        plot(InfoStim.IStimDelay, maxTiming(:,:,indP));
        ylim([0 50]);
        title(['Max amplitude time (msec) with pre-cond power @ ' int2str(InfoStim.PreCondPwr(indP)) '%']);
        legend('Channels 1-2-3', 'Channels 6-7-8');
        ylabel('Maximum Delay (ms)')
        xlabel('ISI (ms)');
        savefig(gcf, [FolderPath 'Max_Pwr_' int2str(InfoStim.PreCondPwr(indP)) '.fig']);
        saveas(gcf, [FolderPath 'Max_Pwr_' int2str(InfoStim.PreCondPwr(indP)) '.png']);
        close(gcf);
    end
end

figure;
disp('Zscore Graph');
plot(timeVect, squeeze(zSig(1,:,:)));
xlim([-TZeros, 450-TZeros]);
LegendStrings = {};
for indE = 1:size(InfoStim.EDesc,1)
    LegendStrings{indE} = ['Pre-Cond Power: ' int2str(InfoStim.EDesc(indE,1))...
        ' ISI: ' int2str(InfoStim.EDesc(indE,2))];
end
legend(LegendStrings);
xlabel('Time (ms)');
ylabel('Amplitude (zscore)');
title('Channels 3-4-5 Raw');
savefig(gcf, [FolderPath 'Raw_Chan345.fig']);
saveas(gcf, [FolderPath 'Raw_Chan345.png']);
close(gcf);

figure;
plot(timeVect, squeeze(zSig(2,:,:)));
xlim([-TZeros, 450-TZeros]);
LegendStrings = {};
for indE = 1:size(InfoStim.EDesc,1)
    LegendStrings{indE} = ['Pre-Cond Power: ' int2str(InfoStim.EDesc(indE,1))...
        ' ISI: ' int2str(InfoStim.EDesc(indE,2))];
end
legend(LegendStrings);
xlabel('Time (ms)');
ylabel('Amplitude (zscore)');
title('Channels 7-8-9 Raw');
savefig(gcf, [FolderPath 'Raw_Chan789.fig']);
saveas(gcf, [FolderPath 'Raw_Chan789.png']);
close(gcf);

end
