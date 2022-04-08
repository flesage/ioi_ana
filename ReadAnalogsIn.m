function out = ReadAnalogsIn(FolderPath, SaveFolder, Infos)

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
if( ~isfield(Infos, 'Stimulation1_Amplitude') )
    Infos.Stimulation1_Amplitude = 5;
end
StimTrig = find((AnalogIN(1:(end-1), 2) < Infos.Stimulation1_Amplitude/2) &...
    (AnalogIN(2:end, 2) >= Infos.Stimulation1_Amplitude/2))+1;

if( ~isempty(StimTrig) && Infos.Stimulation == 1 )
    Period = median(StimTrig(2:end)-StimTrig(1:(end-1)))/Infos.AISampleRate;
    Freq = 1/Period;
    Width = sum(AnalogIN(StimTrig(1):StimTrig(2),2) > 2.5)...
        /(Period*Infos.AISampleRate);
    
    indx_stimLim = find(diff(StimTrig)>20000); 
    NbStim = length(indx_stimLim)+1;
    StimLim = find((AnalogIN(1:(end-1), 2) > Infos.Stimulation1_Amplitude/2) &...
            (AnalogIN(2:end, 2) <= Infos.Stimulation1_Amplitude/2))+1;
    if( NbStim == length(StimTrig) ) %Single Pulse trigged Stims        
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
        % Calculate length of 1st Burst trial in seconds:
        StimLength = (StimLim(indx_stimLim(1)) - StimTrig(1) + StimLim(1) - StimTrig(1))/10000;
        % Get last falling edge of bursts + off time of a burst.
        StimLim = [StimLim(indx_stimLim); StimLim(end)] + StimTrig(2)-StimLim(1); 
        StimTrig = [StimTrig(1); StimTrig(indx_stimLim+1)]; % Get first rising edge of burst.        
        %
        InterStim_min = min((StimTrig(2:end) - StimLim(1:end-1))./10000);
        InterStim_max = max((StimTrig(2:end) - StimLim(1:end-1))./10000);
        %  
        Stim = zeros(length(AnalogIN(:,2)),1);
        if( NbStim > 1 )           
            for indS = 1:length(StimLim)
                Stim(StimTrig(indS):StimLim(indS)) = 1;
            end            
        else
            Stim(StimTrig(1):StimTrig(end)) = 1;
        end
    end
    
    Stim = Stim(CamTrig);
    save([SaveFolder filesep 'StimParameters.mat'],'CamSig', 'CamTrig', 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
elseif( ~isempty(StimTrig) && Infos.Stimulation == 2 )
    NbStimAI = length(StimTrig);
    NbStimCycle = Infos.StimulationRepeat;
    NbStim = sum(contains(fieldnames(Infos), 'Stim')) - 3;
    NbColIll = sum(contains(fieldnames(Infos), 'Illumination'));
    InterFrame = mean(diff(CamTrig));
    if (NbColIll == 1)
        expander = [zeros(1,ceil(NbColIll*InterFrame*1.1)) ones(1,ceil(NbColIll*InterFrame*1.1))]./ceil(NbColIll*InterFrame*1.1);
    else
        expander = [zeros(1,ceil((NbColIll-1)*InterFrame*1.1)) ones(1,ceil((NbColIll-1)*InterFrame*1.1))]./ceil((NbColIll-1)*InterFrame*1.1);
    end
    Stim = conv(AnalogIN(:,2),expander,'same')>0.1;
    Stim = Stim(CamTrig);
    
    if( NbStimAI ~= NbStimCycle*NbStim )
        disp('Acquisition might have been stopped before the end. Not all stimulations were acquired!');
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
    StStart = find(diff(AnalogIN(:,2)) > Infos.Stimulation1_Amplitude/2);
    InterStim_min = mean(StStart(2:end) - StStart(1:(end-1)))/1e4;
    InterStim_max = InterStim_min;
    
    save([SaveFolder filesep 'StimParameters.mat'],'CamTrig', 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
else
    disp('No Stimulations detected. Resting State experiment?');
    Stim = 0;
    StimLength = 0;
    NbStim = 0;
    InterStim_min = 0;
    InterStim_max = 0;
    save([SaveFolder filesep 'StimParameters.mat'], 'CamTrig', 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
end

end

