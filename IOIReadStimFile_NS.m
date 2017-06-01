function IOIReadStimFile_NS(FolderName)

aiFilesList = dir([FolderName filesep 'ai_*.bin']);

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
CamTrig = find((AnalogIN(1:(end-1),1) < 2.5) & (AnalogIN(2:end,1) >= 2.5))+1;

Stim = zeros(size(CamTrig));
StimLength = 0;
NbStim = 0;

StimON = find((AnalogIN(1:(end-1), 2) < 2.5) & (AnalogIN(2:end, 2) >= 2.5))+1;
if( ~isempty(StimON) )
    Period = (StimON(2)-StimON(1))/10000;
    Freq = 1/Period;
    Width = sum(AnalogIN(StimON(1):StimON(2),2) >2.5)/(Period*10000);
    
    StimLim = find(diff(StimON)>20000);
    NbStim = length(StimLim)+1;
    StimLength = round(length(StimON)/(NbStim*Freq));
    InterStim_min = min((StimON(StimLim + 1) - StimON(StimLim))./10000);
    InterStim_max = max((StimON(StimLim + 1) - StimON(StimLim))./10000);
    
    Stim = zeros(length(AnalogIN(:,2)),1);
    Stim(StimON(1):StimON(StimLim(1))) = 1;
    for indS = 2:length(StimLim)
        Stim(StimON(StimLim(indS-1)+1):StimON((StimLim(indS)))) = 1;
    end
    Stim(StimON(StimLim(end)+1):StimON(end)) = 1;
    
    Stim = Stim(CamTrig);
    save([FolderName filesep 'StimParameters.mat'], 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
else 
    disp('No Stimulations detected. Resting State experiment?');
    Stim = 0;
    StimLength = 0;
    NbStim = 0;
    InterStim_min = 0;
    InterStim_max = 0;
    save([FolderName filesep 'StimParameters.mat'], 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
end

end