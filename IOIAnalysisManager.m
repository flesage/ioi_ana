function out = IOIAnalysisManager()

fig = figure('Name', 'IOI Manager',...
    'NumberTitle', 'off',...
    'ToolBar', 'none',...
    'MenuBar', 'none',...
    'Position',[250 100 400 600]);

m_SelectPB = uicontrol('Style', 'pushbutton', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.025 0.925 0.175 0.05],...
    'Callback', @FilesSelection, 'String', 'Select Folder');

m_AnaList = uicontrol('Style', 'listbox', 'Parent', fig,...
    'Units', 'Normalized', 'Position', [0.025 0.025 0.75 0.875],...
    'Background', 'w', 'max', 1);

m_RemovePB = uicontrol('Style', 'pushbutton', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.225 0.925 0.175 0.05],...
    'Callback', @RemoveSelected, 'String', 'Remove');

uicontrol('Style', 'text', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.43 0.9125 0.1 0.05],...
    'String', 'CPUs');
m_CPUsToUse = uicontrol('Style', 'edit', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.53 0.925 0.1 0.05],...
    'String', '1', 'Callback', @CPUChange);
set(m_CPUsToUse, 'Enable', 'off');

m_StartPB = uicontrol('Style', 'pushbutton', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.75 0.925 0.175 0.05],...
    'Callback', @NewStart, 'String', 'Start');

m_binningCB = uicontrol('Style', 'checkbox', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.7875 0.8625 0.20 0.05],...
    'String', 'Bining');

m_EraseCB = uicontrol('Style', 'checkbox', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.7875 0.825 0.20 0.05],...
    'String', 'Clean Folder');

uicontrol('Style', 'text', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.7875 0.775 0.20 0.05],...
    'String', 'Pre-Stim');

m_PreStimDuration = uicontrol('Style', 'listbox', 'Parent', fig,...
    'Units', 'normalized', 'Position', [0.7875 0.75 0.20 0.05],...
    'String', '5 s|10 s');

    function FilesSelection(~, ~, ~)
        RootFolder = uigetdir();
        if( RootFolder == 0 )
            return;
        end
        DirsUnderTest = dir(RootFolder);
        DirsUnderTest(end+1).name = '';
        DirsUnderTest(end).isdir = 1;
        ExpToProcess = {};
        Limit = size(DirsUnderTest, 1);
        indD = 3;
        while( indD <= Limit )
            if( DirsUnderTest(indD).isdir )
                if( ~isempty(DirsUnderTest(indD).name) )
                    FilesCheck = dir([RootFolder filesep DirsUnderTest(indD).name]);
                    SeqFile = 0; AuxFile = 0; InfoFile = 0;
                    for indF = 3:size(FilesCheck,1)
                        if(FilesCheck(indF).isdir)
                            DirsUnderTest(end+1).name = [DirsUnderTest(indD).name filesep FilesCheck(indF).name];
                            DirsUnderTest(end).isdir = 1;
                        elseif( strcmp(FilesCheck(indF).name, 'IOI_scan.seq') )
                            SeqFile = 1;
                        elseif( strcmp(FilesCheck(indF).name, 'IOI_aux.mat') )
                            AuxFile = 1;
                        elseif( strcmp(FilesCheck(indF).name, 'IOI_scaninfo.mat') )
                            InfoFile = 1;
                        end
                    end
                    if( SeqFile*AuxFile*InfoFile )
                        ExpToProcess{end+1} = [RootFolder filesep DirsUnderTest(indD).name];
                    end
                else
                    FilesCheck = dir(RootFolder);
                    SeqFile = 0; AuxFile = 0; InfoFile = 0;
                    for indF = 3:size(FilesCheck,1)
                        if( strcmp(FilesCheck(indF).name, 'IOI_scan.seq') )
                            SeqFile = 1;
                        elseif( strcmp(FilesCheck(indF).name, 'IOI_aux.mat') )
                            AuxFile = 1;
                        elseif( strcmp(FilesCheck(indF).name, 'IOI_scaninfo.mat') )
                            InfoFile = 1;
                        end
                    end
                    if( SeqFile*AuxFile*InfoFile )
                        ExpToProcess{end+1} = [RootFolder filesep DirsUnderTest(indD).name];
                    end
                end
            end
            Limit = size(DirsUnderTest,1); indD = indD + 1;
        end
        set(m_AnaList,'String', ExpToProcess);
    end

    function RemoveSelected(~, ~, ~)
        toBeRemoved = get(m_AnaList, 'Value');
        A = get(m_AnaList, 'String');
        A(toBeRemoved) = [];
        set(m_AnaList, 'String', A);
    end

    function CPUChange(~,~,~)
       if( str2double(get(m_CPUsToUse, 'String')) > 7 )
           set(m_CPUsToUse, 'String', 4);
       elseif( str2double(get(m_CPUsToUse, 'String')) < 1 )
           set(m_CPUsToUse, 'String', 1);
       end
    end

    function NewStart(~,~,~)
        BinData = get(m_binningCB, 'value');
        EraseOldFiles = get(m_EraseCB, 'value');
        PreStimD = get(m_PreStimDuration, 'Value') == 1;
        
        List = get(m_AnaList,'String');
        ToOpen = ones(size(List,1),1);
        
        %%%%%%%%%
        %Expe Folder cleaning...
        %%%%%%%%%
        if( EraseOldFiles )
            for indE = 1:size(List,1)
                Files = dir([List{indE} 'Data*.mat']);
                Files = [Files;  dir([List{indE} 'ROIs.mat'])];
                Files = [Files;  dir([List{indE} 'Hb_*.mat'])];
                arrayfun(@(X) delete([List{indE} X.name]), Files);
            end
        else
            for indE = 1:size(List,1)
                Files = dir([List{indE} 'Data*.mat']);
                Files = [Files;  dir([List{indE} 'Hb_Concentrations.mat'])];
                arrayfun(@(X) delete([List{indE} X.name]), Files);
            end
        end
               
        %%%%%%%%%
        %Speckle?
        %%%%%%%%%
        ToSpeckle = ones(size(List,1),1);
        for indE = 1:size(List,1)
            load([List{indE} filesep 'IOI_scaninfo.mat'],'Signaux');
            ToSpeckle(indE) = any(Signaux(:,5));
            clear Signaux;
        end
        
        %%%%%%%%%
        %Stimulations
        %%%%%%%%%
        for indE = 1:size(List,1)
            if( exist([List{indE} filesep 'StimParameters.mat'],'file') )
                delete([List{indE} filesep 'StimParameters.mat']);
            end
            IOIReadStimFile(List{indE});
            if( PreStimD )
                PreStimLength = 5;
            else
                PreStimLength = 10;
            end
            save([List{indE} filesep 'StimParameters.mat'], 'PreStimLength', '-append');
        end

        %%%%%%%%%
        %Main Loop
        %%%%%%%%%
        NbCPUs = 1; %str2double(get(m_CPUsToUse, 'String'));
        
        %Single Threaded loop:
        if( NbCPUs == 1 )
            h = waitbar(0, 'Starting ...');
            for indR = 1:sum(ToOpen)
                disp('Step 1: Opening data files')
                disp('**************************');
                Fast_OpenIOI( List{indR}, BinData );
                disp('Step 2: Hb Computations')
                disp('**************************');
                Ana_IOI_FullFrame( List{indR}, 0 );
                if( ToSpeckle(indR) )
                    disp('Step 3: Speckle')
                    disp('**************************');
                    Ana_Speckle( List{indR} );
                end
                %IOIFiguresGen(List{indR});
                disp(['Done for:'  List{indR}])
                disp('**************************');
                if( ishghandle(h) )
                    waitbar(indR/sum(ToOpen), h, ['Task ' int2str(indR) ' of ' int2str(sum(ToOpen))]);
                else
                    h = waitbar(indR/sum(ToOpen), ['Task ' int2str(indR) ' of ' int2str(sum(ToOpen))]);
                end
            end
            if( ishghandle(h) )
                close(h);
            end
        else
            %Multi-Threaded loop:
            if( ispc )
                [~, result] = system('tasklist');
                Id = strfind(result,'MATLAB.exe');
            elseif( isunix )
                [~, result] = system('ps -a');
                Id = strfind(result,'MATLAB');
            end
            [p, ~, ~] = fileparts(mfilename('fullpath'));
            
            IdSub = Id;
            indR = 1;
            TaskRunning = 0;
            TaskDone = 0;
            
            h = waitbar(0, 'Starting ...');
            while( any(ToOpen) )
                
                ExecStr = ['matlab -nodisplay -nojvm -nosplash -nodesktop -logfile '...
                    List{indR} filesep 'Log.txt -r "addpath(''' p '''); Fast_OpenIOI('''...
                    List{indR} ''',''' int2str(BinData) '''); Ana_IOI_FullFrame(''' List{indR} ...
                    ''', 0);'];
                
                if( ToSpeckle(indR) )
                    ExecStr = strcat(ExecStr, [' Ana_Speckle(''' List{indR} ''');']);
                end
                
                %ExecStr = strcat(ExecStr, [' IOIFiguresGen(''' List{indR} '''); exit;"  &']);
                
                ExecStr = strcat(ExecStr, [' exit;"  &']);
                
                if( ispc )
                    [~, result] = system('tasklist');
                    NewIds = strfind(result,'MATLAB.exe');
                elseif( isunix )
                    [~, result] = system('ps -a');
                    NewIds = strfind(result,'MATLAB');
                end
                IdSub = NewIds;
                
                system(ExecStr);
                TaskRunning = TaskRunning + 1;
                
                while( size(NewIds,2) == size(IdSub,2) )
                    pause(10);
                    if( ispc )
                        [~, result] = system('tasklist');
                        NewIds = strfind(result,'MATLAB.exe');
                    elseif( isunix )
                        [~, result] = system('ps -a');
                        NewIds = strfind(result,'MATLAB');
                    end
                end
                IdSub = NewIds;
                ToOpen(indR) = 0;
                indR = indR + 1;
                
                CPUsInUse = size(IdSub,2);
                while( CPUsInUse == NbCPUs + 1 )
                    pause(10);
                    if( ispc )
                        [~, result] = system('tasklist');
                        NewIds = strfind(result,'MATLAB.exe');
                    elseif( isunix )
                        [~, result] = system('ps -a');
                        NewIds = strfind(result,'MATLAB');
                    end
                    CPUsInUse = size(NewIds,2);
                end
                IdSub = NewIds;
                TaskDone = TaskDone + (TaskRunning - (CPUsInUse - 1));
                TaskRunning = (CPUsInUse - 1);
                
                if( ishghandle(h) )
                    waitbar(TaskDone/sum(ToOpen), h, ['Exec: ' int2str(TaskRunning) ', Done: ' int2str(TaskDone) ', Waiting: ' int2str(sum(ToOpen) - TaskDone)]);
                else
                    h = waitbar(TaskDone/sum(ToOpen), ['Exec: ' int2str(TaskRunning) ', Done: ' int2str(TaskDone) ', Waiting: ' int2str(sum(ToOpen) - TaskDone)]);
                end
                
            end
            if( ishghandle(h) )
                close(h);
            end
        end
    end

    function Start(~, ~, ~)
        
        BinData = get(m_binningCB, 'value');
        EraseOldFiles = get(m_EraseCB, 'value');
        PreStimD = get(m_PreStimDuration, 'Value') == 1;
        
        List = get(m_AnaList,'String');
        ToOpen = ones(size(List,1),1);
        
        %%%%%%%%%
        %Expe Folder cleaning...
        %%%%%%%%%
        if( EraseOldFiles )
            for indE = 1:size(List,1)
                Files = dir([List{indE} 'Data*.mat']);
                Files = [Files;  dir([List{indE} 'ROIs.mat'])];
                Files = [Files;  dir([List{indE} 'Hb_*.mat'])];
                arrayfun(@(X) delete([List{indE} X.name]), Files);
            end
        else
            for indE = 1:size(List,1)
                Files = dir([List{indE} 'Data*.mat']);
                Files = [Files;  dir([List{indE} 'Hb_Concentrations.mat'])];
                arrayfun(@(X) delete([List{indE} X.name]), Files);
            end
        end
        
        %%%%%%%%%
        %ROI missing?
        %%%%%%%%%
        ToROI = ones(size(List,1),1);
        for indE = 1:size(List,1)
            Files = dir(List{indE});
            ToROI(indE) = isempty(strfind([Files(:).name], 'ROIs'));
            
            %Dimension test:
            if( ~ToROI(indE) )
                NomSeq = [List{indE} filesep 'IOI_scan.seq'];
                ImRes_XY = memmapfile(NomSeq,'Offset',548,'Format','uint32','Repeat',2);
                ImRes_XY = double(ImRes_XY.Data);
                
                if( BinData )
                    ImRes_XY = ImRes_XY./2;
                end
                load([List{indE} filesep 'ROIs.mat'],'ROIs');
                [rY, rX] = size(ROIs{1}.mask);
                if( (ImRes_XY(2) ~= rX) || (ImRes_XY(1) ~= rY) )
                    delete([List{indE} filesep 'ROIs.mat']); 
                    ToROI(indE) = 1;
                end
                clear datFiles matObj sY sX sT ROIs rY rX;
            end            
        end
        
        %%%%%%%%%
        %Is ROI size compatible with data?
        %%%%%%%%%
        idx = find(ToROI);
        if( ~isempty(idx) )
            for indR = 1:length(idx)
                NomSeq = [List{idx(indR)} filesep 'IOI_scan.seq'];
                SizeImage = memmapfile(NomSeq,'Offset',580,'Format','uint32','Repeat',1);
                ImRes_XY = memmapfile(NomSeq,'Offset',548,'Format','uint32','Repeat',2);
                SizeImage = double(SizeImage.Data)/2;
                ImRes_XY = double(ImRes_XY.Data);
                data = memmapfile(NomSeq,'Offset',1024,'Format',{'uint16', [ImRes_XY(1) ImRes_XY(2)], 'framej';'int32', 1, 'datej';'uint16', 2, 'timej';'uint16', SizeImage-4-ImRes_XY(1)*ImRes_XY(2), 'junkj'},'repeat',inf);
                
                incr = 1;
                while( sum(std(single(data.Data(incr).framej),0,2)) == 0 )
                    incr = incr + 1;
                end
                incr = incr + 1;
                
                Map = data.Data(incr).framej;
                if( BinData )
                    Map = imresize(Map, 0.5);
                end                
                
                Tmp = ROI_dialog(Map, List{idx(indR)});
                waitfor(Tmp.m_dialogFig);
                clear Tmp;
            end
        end
        
        %%%%%%%%%
        %Speckle?
        %%%%%%%%%
        ToSpeckle = ones(size(List,1),1);
        for indE = 1:size(List,1)
            load([List{indE} filesep 'IOI_scaninfo.mat'],'Signaux');
            ToSpeckle(indE) = any(Signaux(:,5));
            clear Signaux;
        end
        
        %%%%%%%%%
        %Stimulations
        %%%%%%%%%
        for indE = 1:size(List,1)
            if( exist([List{indE} filesep 'StimParameters.mat'],'file') )
                delete([List{indE} filesep 'StimParameters.mat']);
            end
            IOIReadStimFile(List{indE});
            if( PreStimD )
                PreStimLength = 5;
            else
                PreStimLength = 10;
            end
            save([List{indE} filesep 'StimParameters.mat'], 'PreStimLength', '-append');
        end

        %%%%%%%%%
        %Main Loop
        %%%%%%%%%
        NbCPUs = str2double(get(m_CPUsToUse, 'String'));
        
        %Single Threaded loop:
        if( NbCPUs == 1 )
            h = waitbar(0, 'Starting ...');
            for indR = 1:sum(ToOpen)
                OpenIOI( List{indR}, BinData);
                Ana_IOI( List{indR});
                if( ToSpeckle(indR) )
                    Ana_Speckle( List{indR});
                end
                IOIFiguresGen(List{indR});
                
                if( ishghandle(h) )
                    waitbar(indR/sum(ToOpen), h, ['Task ' int2str(indR) ' of ' int2str(sum(ToOpen))]);
                else
                    h = waitbar(indR/sum(ToOpen), ['Task ' int2str(indR) ' of ' int2str(sum(ToOpen))]);
                end
            end
            if( ishghandle(h) )
                close(h);
            end
        else
            %Multi-Threaded loop:
            [~, result] = system('tasklist');
            Id = strfind(result,'MATLAB.exe');
            [p, ~, ~] = fileparts(mfilename('fullpath'));
            IdSub = Id;
            indR = 1;
            TaskRunning = 0;
            TaskDone = 0;
            h = waitbar(0, 'Starting ...');
            while( any(ToOpen) )
                 ExecStr = ['matlab -nodisplay -nojvm -nosplash -nodesktop -logfile '...
                    List{indR} filesep 'Log.txt -r "addpath(''' p '''); OpenIOI('''...
                    List{indR} ''',' int2str(BinData) '); Ana_IOI(''' List{indR} ...
                    ''');'];
                
                if( ToSpeckle(indR) )
                    ExecStr = strcat(ExecStr, [' Ana_Speckle(''' List{indR} ''');']);
                end
                
                ExecStr = strcat(ExecStr, [' IOIFiguresGen(''' List{indR} '''); exit;"  &']);
                
                system(ExecStr); TaskRunning = TaskRunning + 1;
                
                [~, result] = system('tasklist');
                NewIds = strfind(result,'MATLAB.exe');
                while( size(NewIds,2) == size(IdSub,2) )
                    pause(10);
                    [~, result] = system('tasklist');
                    NewIds = strfind(result,'MATLAB.exe');
                end
                IdSub = NewIds;
                ToOpen(indR) = 0;
                indR = indR + 1;
                CPUsInUse = size(IdSub,2);
                while( CPUsInUse == NbCPUs + 1 )
                    pause(10);
                    [~, result] = system('tasklist');
                    NewIds = strfind(result,'MATLAB.exe');
                    CPUsInUse = size(NewIds,2);
                end
                IdSub = NewIds;
                TaskDone = TaskDone + (TaskRunning - (CPUsInUse - 1));
                TaskRunning = (CPUsInUse - 1);
                if( isgraphics(h) )
                    waitbar(TaskDone/sum(ToOpen), h, ['Exec: ' int2str(TaskRunning) ', Done: ' int2str(TaskDone) ', Waiting: ' int2str(sum(ToOpen) - TaskDone)]);
                else
                    h = waitbar(TaskDone/sum(ToOpen), ['Exec: ' int2str(TaskRunning) ', Done: ' int2str(TaskDone) ', Waiting: ' int2str(sum(ToOpen) - TaskDone)]);
                end
            end
            
            
            if( isgraphics(h) )
                close(h);
            end
        end
    end
end