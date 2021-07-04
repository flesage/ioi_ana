function out = ReadAnalogsIn(FolderPath, Infos)

if( ~strcmp(FolderPath, filesep) )
    FolderPath = strcat(FolderPath, filesep);
end

%List of analog files containing raw data:
aiFilesList = dir([FolderPath 'ai_*.bin']);

%Opening of the files:
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([FolderPath aiFilesList(ind).name],...
        'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, 1e4, Infos.AINChannels, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],Infos.AINChannels);
    AnalogIN = [AnalogIN; tmp];
end
clear tmp ind data aiFilesList;

%CamTrig is on the first channel:
CamSig = AnalogIN(:,1);
CamTrig = find((AnalogIN(1:(end-1),1) < 2.5) & (AnalogIN(2:end,1) >= 2.5))+1;

%StimTrig is on the second channel (except if slave):
StimTrig = find((AnalogIN(1:(end-1), 2) < 2.5) & (AnalogIN(2:end, 2) >= 2.5))+1;

if( ~isempty(StimTrig) && Infos.Stimulation == 1 )
    Period = median(StimTrig(2:end)-StimTrig(1:(end-1)))/Infos.AISampleRate;
    Freq = 1/Period;
    Width = sum(AnalogIN(StimTrig(1):StimTrig(2),2) > 2.5)...
        /(Period*Infos.AISampleRate);
    
    StimLim = find(diff(StimTrig)>20000);
    NbStim = length(StimLim)+1;
    if( NbStim == length(StimTrig) ) %Single Pulse trigged Stims
        StimLim = find((AnalogIN(1:(end-1), 2) > 2.5) &...
            (AnalogIN(2:end, 2) <= 2.5))+1;
        StimLength = mean(StimLim - StimTrig)./Infos.AISampleRate;
        if StimLength < (CamTrig(2) - CamTrig(1))/Infos.AISampleRate
            StimLength = 3*(CamTrig(2) - CamTrig(1))/Infos.AISampleRate;
            StimLim = StimLim + 3*(CamTrig(2) - CamTrig(1));
        end
        InterStim_min = min((StimTrig(2:end) - StimLim(1:(end-1)))./10000);
        InterStim_max = max((StimTrig(2:end) - StimLim(1:(end-1)))./10000);
        Stim = zeros(length(AnalogIN(:,2)),1);
        for indS = 1:NbStim
           Stim(StimTrig(indS):StimLim(indS)) = 1; 
        end
    else %Pulses train Stim
        StimLength = round(length(StimTrig)/(NbStim*Freq));
        InterStim_min = min((StimTrig(StimLim + 1) - StimTrig(StimLim))./10000);
        InterStim_max = max((StimTrig(StimLim + 1) - StimTrig(StimLim))./10000);
        InterStim_min = InterStim_min - StimLength;
        InterStim_max = InterStim_max - StimLength;
    
        Stim = zeros(length(AnalogIN(:,2)),1);
        Stim(StimTrig(1):StimTrig(StimLim(1))) = 1;
        for indS = 2:length(StimLim)
            Stim(StimTrig(StimLim(indS-1)+1):StimTrig((StimLim(indS)))) = 1;
        end
        Stim(StimTrig(StimLim(end)+1):StimTrig(end)) = 1;
    end
    
    Stim = Stim(CamTrig);
    save([FolderPath filesep 'StimParameters.mat'],'CamSig', 'CamTrig', 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
elseif( ~isempty(StimTrig) && Infos.Stimulation == 2 )
    NbStimAI = length(StimTrig);
    NbStimCycle = Infos.StimulationRepeat;
    NbStim = sum(contains(fieldnames(Infos), 'Stim')) - 3;
    NbColIll = sum(contains(fieldnames(Infos), 'Illumination'));
    InterFrame = mean(diff(CamTrig));
    Extenseur = [zeros(1,ceil((NbColIll-1)*InterFrame*1.1)) ones(1,ceil((NbColIll-1)*InterFrame*1.1))]./ceil((NbColIll-1)*InterFrame*1.1);
    Stim = conv(AnalogIN(:,2),Extenseur,'same')>0.1;
    Stim = Stim(CamTrig);
    
    if( NbStimAI ~= NbStimCycle*NbStim )
        disp('Acquisition might have been stoped before the end. Not all stimulations were acquired!');
    end
    
    StimTrig = find(diff(Stim)>0.5)+1;
    StimIDs = [];
    StimDurations = [];
    for indS = 1:NbStim
        eval(['StimIDs = cat(1, StimIDs, Infos.Stim' int2str(indS) '.code);']);
        eval(['StimDurations = cat(1, StimDurations, Infos.Stim' int2str(indS) '.Duration);']);
    end
    
    Stim = zeros(size(CamTrig),'single');
    for ind = 1:length(StimTrig)
        ID = mod(ind-1, NbStim) + 1;
        St = StimTrig(ind);
        En = StimDurations(ID).*Infos.FrameRateHz + St;
        Stim(St:En) = StimIDs(ID);
    end
        
    StimLength = 2*InterFrame;
    NbStim = NbStimAI;
    StStart = find(diff(AnalogIN(:,2)) > 2.5);
    InterStim_min = mean(StStart(2:end) - StStart(1:(end-1)))/1e4;
    InterStim_max = InterStim_min;
    
    save([FolderPath filesep 'StimParameters.mat'],'CamTrig', 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
else
    disp('No Stimulations detected. Resting State experiment?');
    Stim = 0;
    StimLength = 0;
    NbStim = 0;
    InterStim_min = 0;
    InterStim_max = 0;
    save([FolderPath filesep 'StimParameters.mat'], 'CamTrig', 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
end

end

