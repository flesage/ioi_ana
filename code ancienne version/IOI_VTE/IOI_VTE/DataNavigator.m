function out = DataNavigator(varargin)
out = '0';
%%%%%%%
% Vars init
%%%%%%%
h.paths.FolderName = '';
h.flags.saveROIS = false; %flag to know if any changes were made to ROIs
h.flags.saveEvnts = false;
h.flags.IsThereHbO = false;
h.flags.IsThereHbR = false;
h.flags.IsThereHbT = false;
h.flags.IsThereFlow = false;
h.flags.IsThereFluo = false;
h.flags.IsThereGreen = false;
h.flags.IsThereYellow = false;
h.flags.IsThereRed = false;
h.flags.VideoPlaying = false;
h.flags.IsThereSpeckle = false;
%Map
h.data.Map = [];
h.data.EventBuf = 0;
h.data.MasterStim.PreStimLength = 5;
%%%%%%%
% GUI init
%%%%%%%

% Interface Generation
h.ui.fig = figure('Name', 'DataExplorer', 'Position', [100 25 1000 550], 'CloseRequestFcn', @my_closereq);

%%% ROI 
% ROIs view and management
h.ui.ROIs = uipanel('Parent', h.ui.fig, 'Title','ROI','FontSize',12,...
             'Position',[.01 .25 .495 .74]);
h.ui.AddButton = uicontrol('Style','pushbutton','Parent', h.ui.ROIs,...
    'Units', 'normalized', 'Position',[0.01 0.925 0.1 0.075],...
    'String','+', 'Callback', @AddNewROI);
h.ui.RemButton = uicontrol('Style','pushbutton','Parent', h.ui.ROIs,...
    'Units', 'normalized', 'Position',[0.11 0.925 0.1 0.075],...
    'String','-', 'Callback', @RemoveROI);
h.ui.ROIsPanel = uipanel('Parent', h.ui.ROIs,...
    'Position',[0.01 0.09 0.21 0.83]);
h.ui.ROIsListax = axes('Parent', h.ui.ROIsPanel, 'Position',[0 0 1 1], 'XLim', [0 1], 'YLim', [0 1]);
axis(h.ui.ROIsListax,'off');
h.ui.ROIsMap = axes('Parent', h.ui.ROIs, 'Position',[0.23 0.05 0.75 0.92]);

h.ui.SaveROIpb = uicontrol('Style','pushbutton','Parent', h.ui.ROIs,...
    'Units', 'normalized', 'Position',[0.01 0.01 0.1 0.075],...
    'String','save', 'Callback', @SaveROIs);
h.ui.LoadROIpb = uicontrol('Style','pushbutton','Parent', h.ui.ROIs,...
    'Units', 'normalized', 'Position',[0.11 0.01 0.1 0.075],...
    'String', 'load', 'Callback', @LoadROIs);

set(h.ui.AddButton, 'Enable', 'off'); 
set(h.ui.RemButton, 'Enable', 'off');
set(h.ui.SaveROIpb, 'Enable', 'off');
set(h.ui.LoadROIpb, 'Enable', 'off');

%%% Events 
% Events view and management
h.ui.EventsPan = uipanel('Parent', h.ui.fig, 'Title','Events','FontSize',12,...
             'Position',[.51 .25 .485 .74]);
h.ui.EventDispCB = uicontrol('Style','checkbox','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.01 0.90 0.25 0.075],...
    'String','Show events','Value',1,'Callback', @ShowEvents);
h.ui.EventSaveCB = uicontrol('Style','checkbox','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.25 0.90 0.25 0.075],...
    'String','Save Events','Value',1,'Callback', @ShowEvents);
h.ui.EventsDispPan.Container = uipanel('Parent', h.ui.EventsPan, 'Title', 'Selection', 'FontSize', 12,...
              'Position', [0.0 0.5 1.0 0.4]);
h.ui.EventsDispPan.Ax = axes('Parent', h.ui.EventsDispPan.Container, 'Position', ...
                    [0.175 0.15 0.72 0.85]);
h.ui.EventsDispPan.Cbox = uicontrol('Style','checkbox','Parent', h.ui.EventsDispPan.Container,...
    'Units', 'normalized', 'Position', [0.01 0.4 0.1 0.2],...
    'String','1', 'Callback', @OnEditEventsClicked);
h.ui.EventsDispPan.Slider = uicontrol('Parent', h.ui.EventsDispPan.Container,...
             'Style', 'slider',  'Min', 1, 'Max', 10, 'Value', 1,...
             'Units', 'normalized', 'Position', [0.95 0.0 0.05 1], 'SliderStep', [0.01 0.3],...
             'Callback', @MoveEventsDisp); 
h.ui.ROIsSelector = uicontrol('Style','popupmenu','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.20 0.01 0.15 0.075],...
    'String', 'Empty', 'Callback', @SrcChange);
h.ui.ChannelSelector = uicontrol('Style','popupmenu','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.37 0.01 0.15 0.075],...
    'String', 'Empty', 'Callback', @SrcChange);
h.ui.EventsDispRegen = uicontrol('Style','pushbutton','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.01 0.01 0.15 0.075],...
    'String','Gen', 'Callback', @PopulateEvntsDisplay);
h.ui.FilteringOpt = uicontrol('Style','checkbox','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.575 0.01 0.25 0.075],...
    'String','Filter data', 'Callback', @FiltOptChange);
h.ui.GlobalOpt = uicontrol('Style','checkbox','Parent', h.ui.EventsPan,...
    'Units', 'normalized', 'Position',[0.775 0.01 0.25 0.075],...
    'String','Global signal', 'Callback', @GlobalSigChange);
h.ui.EventsMeanPan.Ax = axes('Parent', h.ui.EventsPan, 'Position', ...
                    [0.15 0.15 0.825 0.325]);
h.ui.EventsMeanPan.Tag = axes('Parent', h.ui.EventsPan, 'Position', ...
                    [0.0 0.15 0.10 0.18]); axis(h.ui.EventsMeanPan.Tag, 'off');
text(h.ui.EventsMeanPan.Tag, 0.225, 0.5, 'Mean', 'Rotation', 90); 

set(h.ui.EventsDispPan.Slider, 'Enable', 'off'); 
set(h.ui.ROIsSelector, 'Enable', 'off'); 
set(h.ui.ChannelSelector, 'Enable', 'off'); 
set(h.ui.FilteringOpt, 'Enable', 'off'); 
set(h.ui.GlobalOpt, 'Enable', 'off'); 
set(h.ui.EventsDispRegen, 'Enable', 'off'); 
%%% Videos 
% Display animated sequences 
h.ui.AnimatedPan = uipanel('Parent', h.ui.fig, 'Title','Sequence','FontSize',12,...
             'Position',[.01 .01 .125 .105]);
h.ui.StartVideo = uicontrol('Style','pushbutton','Parent', h.ui.AnimatedPan,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Start', 'Callback', @StartVideo);
h.ui.Overlay = 0;

%%% Data Loading:
h.ui.NewData = uipanel('Parent', h.ui.fig, 'Title','Load DataSet','FontSize',12,...
             'Position',[.01 .125 .125 .105]);
h.ui.LoadData = uicontrol('Style','pushbutton','Parent', h.ui.NewData,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Load', 'Callback', @OpenFolder);

%%% Data graphs:
h.ui.Graph = uipanel('Parent', h.ui.fig, 'Title','Figures','FontSize',12,...
             'Position',[.150 .125 .125 .105]);
h.ui.GGraph = uicontrol('Style','pushbutton','Parent', h.ui.Graph,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Generate', 'Callback', @GenerateGraphs);

%%% Data Export:
h.ui.XLS = uipanel('Parent', h.ui.fig, 'Title','Spreadsheet','FontSize',12,...
             'Position',[.150 .01 .125 .105]);
h.ui.Eport = uicontrol('Style','pushbutton','Parent', h.ui.XLS,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Export', 'Callback', @exportXLS);

%%% Data Export:
h.ui.PStimL = uipanel('Parent', h.ui.fig, 'Title','PreStim','FontSize',12,...
             'Position',[.290 .125 .125 .105]);
h.ui.SetPS = uicontrol('Style','popupmenu','Parent', h.ui.PStimL,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String',{'0s', '2s', '5s', '10s'}, 'Value', 1, 'Callback', @setPreStimLength);

%%% Data path display
h.ui.Pdisp = uipanel('Parent', h.ui.fig, 'Title','Current Path','FontSize',12,...
             'Position',[.430 .125 .200 .105]);
h.ui.dispString = uicontrol('Style','text','Parent', h.ui.Pdisp,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],'string','No data loaded.');

%%% temporal speckle
h.ui.Speckdisp = uipanel('Parent', h.ui.fig, 'Title','Speckle','FontSize',12,...
             'Position',[.430 .01 .200 .105]);
h.ui.SpeckleButton = uicontrol('Style','pushbutton','Parent', h.ui.Speckdisp,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],'string','Save'...
    ,'Callback',@Temporal_Speckle);


%%% Intensity checkup:
h.ui.Icheck = uipanel('Parent', h.ui.fig, 'Title','Intensity','FontSize',12,...
             'Position',[.290 .01 .125 .105]);
h.ui.IChckButton = uicontrol('Style','pushbutton','Parent', h.ui.Icheck,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String', 'Check', 'Callback', @validateIntensity);

%%%%%%%%%%%%%%%
%Fonctions & Callbacks:
%%%%%%%%%%%%%%%
    function validateIntensity(~,~,~)
        cmap = gray(4096);
        cmap(1:512,1) = 1;
        cmap(1:512,2:3) = 0;
        cmap((end-511):end,1:2) = 0;
        cmap((end-511):end,3) = 1;
        if( h.flags.IsThereGreen )
            Datptr = matfile([h.paths.FolderName filesep 'Data_green.mat']);
            d = memmapfile(Datptr.datFile, 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Green Channel');
            colormap(cmap);
            axis image; axis off; colorbar;
        end
        if( h.flags.IsThereYellow )
            Datptr = matfile([h.paths.FolderName filesep 'Data_yellow.mat']);
            d = memmapfile(Datptr.datFile, 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Yellow Channel');
            colormap(cmap)
            axis image; axis off; colorbar;
        end
        if( h.flags.IsThereRed )
            Datptr = matfile([h.paths.FolderName filesep 'Data_red.mat']);
            d = memmapfile(Datptr.datFile, 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Red Channel');
            colormap(cmap)
            axis image; axis off; colorbar;
        end
        if( h.flags.IsThereFlow )
            d = memmapfile([h.paths.FolderName filesep 'sChan.dat'], 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Flow Channel');
            colormap(cmap)
            axis image; axis off; colorbar;
        end 
        if( h.flags.IsThereFluo )
            d = memmapfile([h.paths.FolderName filesep 'fChan.dat'], 'Format', 'single');
            d = d.Data(1:(5*length(h.data.Map(:))));
            d = reshape(d, size(h.data.Map,1), size(h.data.Map,2), []);
            figure;
            imagesc(mean(d,3),[0 4095]);
            title('Fluo Channel');
            colormap(cmap)
            axis image; axis off; colorbar;
        end 
    end

    function setPreStimLength(Src,~,~)
        sID = get(h.ui.SetPS, 'Value');
        if( sID == 1 )
            h.data.MasterStim.PreStimLength = 0.11; 
        elseif( sID == 2 )
            h.data.MasterStim.PreStimLength = 2; 
        elseif( sID == 3 )
            h.data.MasterStim.PreStimLength = 5; 
        elseif( sID == 4 )
            h.data.MasterStim.PreStimLength = 10;
        end
        OpenFolder(Src);
    end

    function GenerateGraphs(~,~,~)
            prompt = {'Green Colormap minimum:','Green Colormap maximum:',...
                'Red Colormap minimum:','Red Colormap maximum:',...
                'Yellow Colormap minimum:','Yellow Colormap maximum:',...
                'HbO Colormap minimum:','HbO Colormap maximum:',...
                'HbR Colormap minimum:','HbR Colormap maximum:',...
                'HbT Colormap minimum:','HbT Colormap maximum:',...
                'Flow Colormap minimum:','Flow Colormap maximum:'};
            dlg_title = 'Colormap Limits';
            num_lines = 1;
            defaultans = {'0.99','1.01','0.99','1.01','0.99','1.01','-5','5','-5','5','-5','5','0.75','1.25'};
            cm_answers = inputdlg(prompt,dlg_title,num_lines,defaultans);
        %Waiting Dlg...                
        GraphsDlg = dialog('Position',[500 500 250 150],'Name','Graphs');
        GraphsStr = uicontrol('Parent', GraphsDlg, 'Style','text',...
            'Position',[20 80 210 40], 'String', 'Saving Map...');
        pause(0.1);
        Map = zeros(size(h.data.Map,1), size(h.data.Map,2), 3);
        G = Map(:,:,1);
        if( h.flags.IsThereGreen )
            G = reshape(h.data.gDatPtr.Data(1:length(h.data.Map(:))), size(h.data.Map));
        end
        Y = Map(:,:,1);
        if( h.flags.IsThereYellow )
            Y = reshape(h.data.yDatPtr.Data(1:length(h.data.Map(:))), size(h.data.Map));
        end
        R = Map(:,:,1);
        if( h.flags.IsThereRed )
            R = reshape(h.data.rDatPtr.Data(1:length(h.data.Map(:))), size(h.data.Map));
        end
        Map(:,:,1) = (84.46*R + 77.5326*Y + 27.6758*G);
        Map(:,:,2) = (25.07*R + 52.1235*Y + 143.7679*G);
        Map(:,:,3) = (17.95*R + 21.7241*Y + 66.4394*G);
        
        Map = bsxfun(@rdivide, Map, permute(max(reshape(Map,[],3), [], 1), [3 1 2]));
        Map = bsxfun(@rdivide, Map, permute([1 1.275 1.35], [3 1 2]));
        fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
        imshow(Map);
        axis image; axis off;
        title('Imaged Area');
        figName = 'RawMap';
        saveas(fig, [h.paths.Graphs figName], 'png');
        close(fig);
        
       
        eLen = floor(h.data.AcqFreq*...
                (h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min));
        T = linspace(-h.data.MasterStim.PreStimLength, h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min - h.data.MasterStim.PreStimLength, eLen);
        %for each ROI
        Map = zeros(size(h.data.Map));
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
            Map = Map + mask;
            set(GraphsStr, 'String', ['ROI #' int2str(indR) ' Map saving ...']);
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            image(Map); hold 'on';
            image(~repmat(mask,1,1,3), 'AlphaData', ...
                mask.*0.10 + (~imerode(mask,strel('diamond',1))&mask)*0.9);
            title([h.data.ROIs{indR}.name ' Map']);
            axis image; axis off;
            figName = [h.data.ROIs{indR}.name ' Map'];
            saveas(fig, [h.paths.Graphs figName], 'png');
            close(fig);
            
            AccumGreen = zeros(eLen, length(h.data.EvntList),'single');
            AccumYellow = zeros(eLen, length(h.data.EvntList),'single');
            AccumRed = zeros(eLen, length(h.data.EvntList),'single');
            AccumHbO = zeros(eLen, length(h.data.EvntList),'single');
            AccumHbR = zeros(eLen, length(h.data.EvntList),'single');
            AccumFlow = zeros(eLen, length(h.data.EvntList),'single');
            for indE = 1:length(h.data.EvntList)
                set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ', Colours...']);
                
                fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                ax = axes('Parent', fig);
                hold(ax,'on');
                maxi = 1.01;
                mini = 0.99;

                if( h.flags.IsThereGreen )
                    %Open
                    Datptr = matfile([h.paths.FolderName filesep 'Data_green.mat']);
                    d = memmapfile(Datptr.datFile, 'Format', 'single');
                    d = d.Data((length(h.data.Map(:))*(h.data.G_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.G_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:floor(5*h.data.AcqFreq)));
                     Pend = median(d((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.MasterStim.PreStimLength/2)));
                     d = d./L;
                
                     if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     %Plot
                     plot(ax, T, d, 'Color', [0.0 0.75 0.0], 'LineWidth', 2);
                     AccumGreen(:,indE) = d;
                end
                if( h.flags.IsThereYellow )
                    %Open
                    Datptr = matfile([h.paths.FolderName filesep 'Data_yellow.mat']);
                    d = memmapfile(Datptr.datFile, 'Format', 'single');
                    d = d.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:floor(5*h.data.AcqFreq)));
                     Pend = median(d((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.MasterStim.PreStimLength/2)));
                     d = d./L;
                
                      if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     
                     %Plot
                     plot(ax, T, d, 'Color', [0.75 0.75 0.0], 'LineWidth', 2);
                     AccumYellow(:,indE) = d;
                end
                if( h.flags.IsThereRed )
                    %Open
                    Datptr = matfile([h.paths.FolderName filesep 'Data_red.mat']);
                    d = memmapfile(Datptr.datFile, 'Format', 'single');
                    d = d.Data((length(h.data.Map(:))*(h.data.R_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.R_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:floor(5*h.data.AcqFreq)));
                     Pend = median(d((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.MasterStim.PreStimLength/2)));
                     d = d./L;
                
                      if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     
                     %Plot
                     plot(ax, T, d, 'Color', [0.75 0.0 0.0], 'LineWidth', 2);
                     AccumRed(:,indE) = d;
                end
                clear d;
                box(ax,'on');
                set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                title(['{\Delta}Reflectance over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                ylabel('Normalized Reflectance');
                xlabel('Time (sec)');
                xlim([T(1), T(end)]);
                ylim([mini, maxi]);
                line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                line(ax, [0 0], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
                line(ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
                %Save figure for colours:
                figName = ['Colours_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                saveas(fig, [h.paths.Graphs figName], 'png');
                close(fig);
                
                %Figure for Flow:
                if( h.flags.IsThereFlow )
                    set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ', Flow...']);
                    fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                    ax = axes('Parent', fig);
                    hold(ax,'on');
                    maxi = 1.25;
                    mini = 0.75;
                    %Open
                    dF = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
                    dF = dF.Data((length(h.data.Map(:))*(h.data.F_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.F_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                    
                     %Detrend
                     Pstart = median(dF(1:floor(5*h.data.AcqFreq)));
                     Pend = median(dF((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.MasterStim.PreStimLength/2)));
                     dF = dF./L;
                                    
                     if( max(dF) > maxi) 
                         maxi = max(dF);
                     end
                     if( min(dF) < mini) 
                         mini = min(dF);
                     end
                     
                     %Plot
                     plot(ax, T, dF, 'Color', [0.0 0.0 0.0], 'LineWidth', 2);
                     box(ax,'on');
                     set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                     title(['{\Delta}Flow over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                     ylabel('{\Delta}Flow');
                     xlabel('Time (sec)');
                     
                     line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                     line(ax, [0 0], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     line(ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     xlim([T(1), T(end)]);
                     ylim([mini, maxi]);
                     %Save figure for colours:
                     figName = ['Flow_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                     saveas(fig, [h.paths.Graphs figName], 'png');
                     close(fig);
                     AccumFlow(:,indE) = dF;
                     clear dF;
                end
                
                if( h.flags.IsThereFluo )
                    set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ', Fluo...']);
                    fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                    ax = axes('Parent', fig);
                    hold(ax,'on');
                    maxi = 1.25;
                    mini = 0.75;
                    %Open
                    dF = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
                    dF = dF.Data((length(h.data.Map(:))*(h.data.F_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.F_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                    
                     %Detrend
                     Pstart = median(dF(1:floor(5*h.data.AcqFreq)));
                     Pend = median(dF((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.MasterStim.PreStimLength/2)));
                     dF = dF./L;
                                    
                     if( max(dF) > maxi) 
                         maxi = max(dF);
                     end
                     if( min(dF) < mini) 
                         mini = min(dF);
                     end
                     
                     %Plot
                     plot(ax, T, dF, 'Color', [0.0 0.0 0.0], 'LineWidth', 2);
                     box(ax,'on');
                     set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                     title(['{\Delta}Flow over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                     ylabel('{\Delta}Flow');
                     xlabel('Time (sec)');
                     
                     line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                     line(ax, [0 0], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     line(ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     xlim([T(1), T(end)]);
                     ylim([mini, maxi]);
                     %Save figure for colours:
                     figName = ['Flow_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                     saveas(fig, [h.paths.Graphs figName], 'png');
                     close(fig);
                     AccumFlow(:,indE) = dF;
                     clear dF;
                end
                
                %Figure for Hbs:
                if( h.flags.IsThereHbO && h.flags.IsThereHbR )
                    set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ' Hbs...']);
                    fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                    ax = axes('Parent', fig);
                    hold(ax,'on');
                    maxi = 5;
                    mini = -5;
                    %Open
                    dO = memmapfile(h.data.HBinfos.datFileHbO, 'Format', 'single');
                    dO = dO.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                    dR = memmapfile(h.data.HBinfos.datFileHbR, 'Format', 'single');
                    dR = dR.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                    dO = reshape(dO, [], eLen);
                    dO = mean(dO(mask(:) == 1, :), 1);
                    dR = reshape(dR, [], eLen);
                    dR = mean(dR(mask(:) == 1, :), 1);
                    
                     %Detrend
                     Pstart = median(dO(1:floor(5*h.data.AcqFreq)));
                     Pend = median(dO((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.MasterStim.PreStimLength/2)));
                     dO = dO - L;
                     Pstart = median(dR(1:floor(5*h.data.AcqFreq)));
                     Pend = median(dR((end-floor(5*h.data.AcqFreq)):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.MasterStim.PreStimLength/2)));
                     dR = dR - L;
                
                     if( max(dO) > maxi) 
                         maxi = max(dO);
                     end
                     if( min(dO) < mini) 
                         mini = min(dO);
                     end
                     if( max(dR) > maxi)
                         maxi = max(dR);
                     end
                     if( min(dR) < mini)
                         mini = min(dR);
                     end
                     
                     %Plot
                     plot(ax, T, dR, 'Color', [0.0 0.0 1.0], 'LineWidth', 2);
                     plot(ax, T, dO, 'Color', [1.0 0.0 0.0], 'LineWidth', 2);
                     plot(ax, T, dR + dO, 'Color', [0.0 1.0 0.0], 'LineWidth', 2);
                     box(ax,'on');
                     set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                     title(['{\Delta}Hb over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                     ylabel('{\Delta}Hb Concentration');
                     xlabel('Time (sec)');
                     
                     line(ax, [T(1) T(end)], [0 0], 'Color', 'k', 'LineStyle',':');
                     line(ax, [0 0], [-5 5], 'Color', 'k', 'LineStyle','--');
                     line(ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength], [-5 5], 'Color', 'k', 'LineStyle','--');
                     xlim([T(1), T(end)]);
                     ylim([mini, maxi]);
                     %Save figure for colours:
                     figName = ['Hb_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                     saveas(fig, [h.paths.Graphs figName], 'png');
                     close(fig);
                     AccumHbO(:,indE) = dO;
                     AccumHbR(:,indE) = dR;
                     clear dO dR;
                end
                
            end
            set(GraphsStr, 'String', ['ROI #' int2str(indR) ' average graphs...']);
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            hold(ax,'on');
            if(h.flags.IsThereGreen)
                errorbar(ax, T, mean(AccumGreen,2), std(AccumGreen,1,2)/sqrt(indE),...
                    'Color', [0.0 0.75 0.0], 'LineWidth', 2);
            end
            if(h.flags.IsThereYellow)
                errorbar(ax, T, mean(AccumYellow,2), std(AccumYellow,1,2)/sqrt(indE),...
                    'Color', [0.75 0.75 0.], 'LineWidth', 2);
            end
            if(h.flags.IsThereRed)
                errorbar(ax, T, mean(AccumRed,2), std(AccumRed,1,2)/sqrt(indE),...
                    'Color', [0.75 0.0 0.0], 'LineWidth', 2);
            end
            box(ax,'on');
            set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
            title(['Mean reflectance over ' h.data.ROIs{indR}.name]);
            ylabel('Reflectance');
            xlabel('Time (sec)');
            line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
            line(ax, [0 0], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
            line(ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
            axis tight;
            %Save figure for colours:
            figName = ['MeanColours_' h.data.ROIs{indR}.name];
            saveas(fig, [h.paths.Graphs figName], 'png');
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            hold(ax,'on');
            if(h.flags.IsThereHbO)
                errorbar(ax, T, mean(AccumHbO,2), std(AccumHbO,1,2)/sqrt(indE),...
                    'Color', [1. 0.0 0.0], 'LineWidth', 2);
                errorbar(ax, T, mean(AccumHbR,2), std(AccumHbR,1,2)/sqrt(indE),...
                    'Color', [0.0 0.0 1.], 'LineWidth', 2);
                errorbar(ax, T, mean(AccumHbR+AccumHbO,2), std(AccumHbR+AccumHbO,1,2)/sqrt(indE),...
                    'Color', [0.0 1. 0.0], 'LineWidth', 2);
            end
            box(ax,'on');
            set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
            title(['{\Delta}Hb over ' h.data.ROIs{indR}.name]);
            ylabel('{\Delta}Hb Concentration');
            xlabel('Time (sec)');
            line(ax, [T(1) T(end)], [0 0], 'Color', 'k', 'LineStyle',':');
            line(ax, [0 0], [-5 5], 'Color', 'k', 'LineStyle','--');
            line(ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength], [-5 5], 'Color', 'k', 'LineStyle','--');
            axis tight;
            %Save figure for colours:
            figName = ['MeanHbs_' h.data.ROIs{indR}.name];
            saveas(fig, [h.paths.Graphs figName], 'png');
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            hold(ax,'on');
            if( h.flags.IsThereFlow )
                errorbar(ax, T, mean(AccumFlow,2), std(AccumFlow,1,2)/sqrt(indE),...
                    'Color', 'k', 'LineWidth', 2);
                title(['Mean {\Delta}Flow over ' h.data.ROIs{indR}.name]);
                ylabel('{\Delta}Flow');
                xlabel('Time (sec)');
                
                line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                line(ax, [0 0], [0.9 1.1], 'Color', 'k', 'LineStyle','--');
                line(ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength], [0.9 1.1], 'Color', 'k', 'LineStyle','--');
                axis tight;
                %Save figure for colours:
                figName = ['MeanFlow_' h.data.ROIs{indR}.name];
                saveas(fig, [h.paths.Graphs figName], 'png');
                close(fig);
            end
        end
        set(GraphsStr, 'String', 'Video Sequences...');
        
        Map(Map > 1) = 1;
        Map = Map ==1;
        % getting ROIs boundaries
        se = strel('disk',3,4);
        for indk = 1:size(h.data.ROIs,2)
            er = imerode(h.data.ROIs{indk}.mask,se);
            diff = h.data.ROIs{indk}.mask - er;
            nb_pos = 0;
            for i = 1:size(diff,1)
                for j = 1:size(diff,2)
                    if(diff(i,j))
                        nb_pos = nb_pos +1;
                        pos(nb_pos,1,indk) = i;
                        pos(nb_pos,2,indk) = j;
                    end
                end
            end
        end
        
        %Video sequence of each Colour
        if( h.flags.IsThereGreen )
            fig = figure('InvertHardcopy','off','Color',[1 1 1],'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'gChan.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                d = d./L;
                Accum = Accum + d;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'RelGreen.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [str2num(cm_answers{1}) str2num(cm_answers{2})]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['Green Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'gChan.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Accum = Accum + d;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            Accum = reshape(Accum,[],size(Accum,3));
            
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'RelGreen.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['Green Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
        end
        
        if( h.flags.IsThereYellow )
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'yChan.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                d = d./L;
                Accum = Accum + d;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum(Map(:),:),[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'RelYellow.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                 hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [str2num(cm_answers{5}) str2num(cm_answers{6})]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['Yellow Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off')
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
            
             Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'yChan.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Accum = Accum + d;
            end
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            Accum = Accum./sum(h.data.EvntList);
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum,[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'AbsYellow.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['Yellow Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
        end
        
        if( h.flags.IsThereRed )
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'rChan.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                d = d./L;
                Accum = Accum + d;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum,[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'RelRed.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                 hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [str2num(cm_answers{3}) str2num(cm_answers{4})]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['Red Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
            
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'rChan.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Accum = Accum + d;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum,[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'AbsRed.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [P(1) P(2)]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['Red Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame);
            end
            close(v);
            close(fig);
        end
        
        %Video sequence of HbO, HbR & HbT
        if( h.flags.IsThereHbO )
            %section for HbO
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            AccumO = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'HbO.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                d = d - L;
                AccumO = AccumO + d;
            end
            AccumO = AccumO./sum(h.data.EvntList);
            AccumO = imfilter(AccumO, fspecial('gaussian',5,3),'same','symmetric');
            AccumO = reshape(AccumO,[],size(AccumO,3));
            P = prctile(reshape(AccumO,[],1),[5 95]);
            AccumO = reshape(AccumO, size(Map,1),[],size(AccumO,2));
            v = VideoWriter([h.paths.Graphs filesep 'HbO.avi']);
            v.FrameRate = 7.5;
            open(v);
            
            for indF = 1:size(AccumO,3)
                imagesc(ax, squeeze(AccumO(:,:,indF)));
                hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [str2num(cm_answers{7}) str2num(cm_answers{8})]);
                axis(ax, 'image', 'off');
                colormap(ax ,jet(256));
                colorbar(ax);
                title(ax,['HbO Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame); 
            end
            close(v);
            close(fig);
            
            % section for HbR
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            AccumR = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'HbR.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                d = d - L;
                AccumR = AccumR + d;
            end
            AccumR = AccumR./sum(h.data.EvntList);
            AccumR = imfilter(AccumR, fspecial('gaussian',5,3),'same','symmetric');
            AccumR = reshape(AccumR,[],size(AccumR,3));
            P = prctile(reshape(AccumR,[],1),[5 95]);
            AccumR = reshape(AccumR, size(Map,1),[],size(AccumR,2));
            v = VideoWriter([h.paths.Graphs filesep 'HbR.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(AccumR,3)
                imagesc(ax, squeeze(AccumR(:,:,indF)));
                hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [str2num(cm_answers{9}) str2num(cm_answers{10})]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['HbR Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame); 
            end
            close(v);
            close(fig);
            
            % section for HbT
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dato = memmapfile([h.paths.FolderName filesep 'HbO.dat'], 'Format', 'single');
                datr = memmapfile([h.paths.FolderName filesep 'HbR.dat'], 'Format', 'single');
                d = dato.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = d + datr.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                d = d - L;
                Accum = Accum + d;
            end
            
            Accum = Accum./sum(h.data.EvntList);
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum,[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'HbT.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [str2num(cm_answers{11}) str2num(cm_answers{12})]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['HbT Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame); 
            end
            close(v);
            close(fig);
        end
        
        %Video sequence of flow
        if( h.flags.IsThereFlow )
            fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
            ax = axes('Parent', fig);
            
            Accum = zeros([size(Map), length(T)], 'single');
            for indE = find(h.data.EvntList)
                dat = memmapfile([h.paths.FolderName filesep 'Flow.dat'], 'Format', 'single');
                d = dat.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                d = reshape(d, size(Map,1), size(Map,2), []);
                Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                d = d./L;
                Accum = Accum + d;
            end
            Accum = Accum./sum(h.data.EvntList);
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            Accum = reshape(Accum,[],size(Accum,3));
            P = prctile(reshape(Accum,[],1),[5 95]);
            Accum = reshape(Accum, size(Map,1),[],size(Accum,2));
            v = VideoWriter([h.paths.Graphs filesep 'Flow.avi']);
            v.FrameRate = 7.5;
            open(v);
            for indF = 1:size(Accum,3)
                imagesc(ax, squeeze(Accum(:,:,indF)));
                hold(ax,'on');
                for i = 1:size(h.data.ROIs,2)
                    temp_pos = squeeze(pos(:,:,i));
                    plot(ax,temp_pos(:,2),temp_pos(:,1),'k.','MarkerSize',1);
                end
                set(ax, 'CLim', [str2num(cm_answers{13}) str2num(cm_answers{14})]);
                axis(ax, 'image', 'off');
                colormap(ax,jet(256));
                colorbar(ax);
                title(ax,['Flow Intensity at: ' num2str(T(indF)) ' sec']);
                hold(ax,'off');
                frame = getframe(fig);
                writeVideo(v,frame); 
            end
            close(v);
            close(fig);
        end
        delete(GraphsDlg);
    end

    function exportXLS(~,~,~)
        %Waiting Dlg...                
        ExportDlg = dialog('Position',[500 500 250 150],'Name','Export');
        uicontrol('Parent', ExportDlg, 'Style','text',...
            'Position',[20 80 210 40], 'String', 'Exporting data...');
        pause(0.1);
        
        eLen = floor(h.data.AcqFreq*...
                (h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min));
        T = linspace(-h.data.MasterStim.PreStimLength, h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min - h.data.MasterStim.PreStimLength, eLen);
        array = zeros(eLen, size(h.data.ROIs,2)*length(h.data.EvntList)+1, 'single');
        pos = strfind(h.paths.FolderName,filesep);
        name = h.paths.FolderName(pos(end)+1:end);
        filename = [h.paths.FolderName filesep 'Graphs' filesep name '.xls'];
        array(:, 1) = T; names = {'T'};
        if(h.flags.IsThereHbO)
            for indR = 1:size(h.data.ROIs,2)
                mask = h.data.ROIs{indR}.mask;
                for indE = 1:length(h.data.EvntList)
                    dO =   h.data.hoDatPtr.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                    dO = reshape(dO, [], eLen);
                    dO = mean(dO(mask(:) == 1, :), 1);
                    array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dO;
                    names = cat(1, names, {['R' int2str(indR) 'E' int2str(indE)]});
                end
            end
            Tabl = array2table(array, 'VariableNames', names);
            writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','HbO','Range','A1');
        end
        if(h.flags.IsThereHbR)
            for indR = 1:size(h.data.ROIs,2)
                mask = h.data.ROIs{indR}.mask;
                for indE = 1:length(h.data.EvntList)
                    dR =   h.data.hrDatPtr.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                    dR = reshape(dR, [], eLen);
                    dR = mean(dR(mask(:) == 1, :), 1);
                   array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dR;
                end
            end
            Tabl = array2table(array, 'VariableNames', names);
            writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','HbR','Range','A1');
        end
        if(h.flags.IsThereFlow)
            for indR = 1:size(h.data.ROIs,2)
                mask = h.data.ROIs{indR}.mask;
                for indE = 1:length(h.data.EvntList)
                    dF =   h.data.fDatPtr.Data((length(h.data.Map(:))*(h.data.F_eflag(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(h.data.F_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                   array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dF;
                end
            end
            Tabl = array2table(array, 'VariableNames', names);
            writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','Flow','Range','A1');
        end
        %Green
        if(h.flags.IsThereGreen)
            for indR = 1:size(h.data.ROIs,2)
                mask = h.data.ROIs{indR}.mask;
                for indE = 1:length(h.data.EvntList)
                    dF =   h.data.gDatPtr.Data((length(h.data.Map(:))*(h.data.G_eflag(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(h.data.G_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                   array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dF;
                end
            end
            Tabl = array2table(array, 'VariableNames', names);
            writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','Green','Range','A1');
        end
        %Yellow
        if(h.flags.IsThereYellow)
            for indR = 1:size(h.data.ROIs,2)
                mask = h.data.ROIs{indR}.mask;
                for indE = 1:length(h.data.EvntList)
                    dF =   h.data.yDatPtr.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                   array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dF;
                end
            end
            Tabl = array2table(array, 'VariableNames', names);
            writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','Yellow','Range','A1');
        end
        if(h.flags.IsThereRed)
            for indR = 1:size(h.data.ROIs,2)
                mask = h.data.ROIs{indR}.mask;
                for indE = 1:length(h.data.EvntList)
                    dF =   h.data.rDatPtr.Data((length(h.data.Map(:))*(h.data.R_eflag(indE) - 1) + 1):...
                            (length(h.data.Map(:))*(h.data.R_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                   array(:, (indR-1)*length(h.data.EvntList) + indE+1) = dF;
                end
            end
            Tabl = array2table(array, 'VariableNames', names);
            writetable(Tabl, filename, 'FileType', 'spreadsheet', 'Sheet','Red','Range','A1');
        end
        delete(ExportDlg);
        disp('Done exporting Excel file.');
    end

    function OpenFolder(Src, ~, ~)
        if( strcmp(Src.String, 'Load') )
            if(isfield(h.data,'ROIs'))
                if(~isempty(h.data.ROIs))
                    selection = questdlg('Do you want to save modified ROIs list?',...
                        'Before loading...',...
                        'Yes','No','Yes');
                    if( strcmp(selection, 'Yes') )
                        SaveROIs();
                        SaveEvnts();
                    end 
                end
            end
            h.paths.FolderName = uigetdir();
        elseif( strcmp(h.paths.FolderName, '') )
            return;
        end
        
        h.flags.bsaveROIS = false; %flag to know if any changes were made to ROIs
        h.flags.bsaveEvnts = false;
        h.flags.MasterStim = false;
        h.flags.SlaveStim = false;
        h.flags.IsThereHbO = false;
        h.flags.IsThereHbR = false;
        h.flags.IsThereHbT = false;
        h.flags.IsThereFlow = false;
        h.flags.IsThereGreen = false;
        h.flags.IsThereRed = false;
        h.flags.IsThereYellow = false;
        h.flags.IsThereSpeckle = false;
               
        % Files Path
        h.paths.HbFile = [h.paths.FolderName filesep 'Data_Hbs.mat'];
        h.paths.ROIsFile = [h.paths.FolderName filesep 'ROIs.mat'];
        h.paths.EVNTsFile = [h.paths.FolderName filesep 'Events.mat'];
        h.paths.MasterStimProto = [h.paths.FolderName filesep 'MasterStimParameters.mat'];
        h.paths.SlaveStimProto = [h.paths.FolderName filesep 'SlaveStimParameters.mat'];
        h.paths.Graphs = [h.paths.FolderName filesep 'Graphs' filesep];
        h.paths.Flow = [h.paths.FolderName filesep 'Flow_infos.mat'];
        h.paths.Fluo = [h.paths.FolderName filesep 'Data_Fluo.mat'];
        
        if( strcmp(Src.String, 'Load') )
            if( exist(h.paths.Graphs , 'dir') )
                ButtonName = questdlg('A folder containing figures already exist. Do you want to overwrite it?', ...
                    'Figures folder', ...
                    'Yes', 'Change', 'Yes');
                switch ButtonName
                    case 'Yes'
                        disp('Erasing old figures...');
                        delete([h.paths.Graphs '*.*']);
                    case 'Change'
                        dname = uigetdir(h.paths.FolderName);
                        h.paths.Graphs = dname;
                end % switch
                
            else
                mkdir(h.paths.Graphs);
            end
        end
        
        %Waiting Dlg...                
        opendlg = dialog('Position',[500 500 250 150],'Name','Loading...');
        uicontrol('Parent', opendlg, 'Style','text',...
            'Position',[20 80 210 40], 'String', 'Loading data. Please wait...');
        pause(0.1);
        
        %Load stimulation parameters
        if(  exist(h.paths.MasterStimProto, 'file') )
            h.data.MasterStim = load(h.paths.MasterStimProto);
            sID = get(h.ui.SetPS, 'Value');
            if( sID == 1 )
                h.data.MasterStim.PreStimLength = 0.11;
            elseif( sID == 2 )
                h.data.MasterStim.PreStimLength = 2;
            elseif( sID == 3 )
                h.data.MasterStim.PreStimLength = 5;
            elseif( sID == 4 )
                h.data.MasterStim.PreStimLength = 10;
            end
            h.flags.MasterStim = true;
        else
            h.flags.MasterStim = false;
            h.data.MasterStim.NbStim = 1;
            h.data.MasterStim.PreStimLength = 0;
            h.data.MasterStim.StimLength = 0;
            h.data.MasterStim.InterStim_min = 0;
        end
        
        if(  exist(h.paths.SlaveStimProto, 'file') )
            h.data.SlaveStim = load(h.paths.SlaveStimProto);
            sID = get(h.ui.SetPS, 'Value');
            if( sID == 1 )
                h.data.SlaveStim.PreStimLength = 0.11;
            elseif( sID == 2 )
                h.data.SlaveStim.PreStimLength = 2;
            elseif( sID == 3 )
                h.data.SlaveStim.PreStimLength = 5;
            elseif( sID == 4 )
                h.data.SlaveStim.PreStimLength = 10;
            end
            h.flags.SlaveStim = true;
        else
            h.flags.SlaveStim = false;
            h.data.SlaveStim.NbStim = 1;
            h.data.SlaveStim.PreStimLength = 0;
            h.data.SlaveStim.StimLength = 0;
            h.data.SlaveStim.InterStim_min = 0;
        end
        
        h.data.AcqFreq = 0;
        Ts = 1e6;
        if( exist(h.paths.Flow, 'file') )
             h.flags.IsThereFlow = true;
             
             h.data.fInfo = matfile(h.paths.Flow);
             h.data.AcqFreq = h.data.fInfo.Freq;
             h.data.fDatPtr = memmapfile([h.paths.FolderName filesep 'Flow.dat'], 'Format', 'single');
             nframes = h.data.fInfo.datLength;
             Ts = min(nframes, Ts);
             
             stim = h.data.fInfo.Stim;
             if( size(stim,2) > size(stim,1) )
                 stim = stim';
             end
             idxS = floor(h.data.AcqFreq*h.data.MasterStim.PreStimLength);
             if( idxS < 1 )
                 idxS = 1;
             end
             Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
            
             if( isempty(Start) )
                 h.data.F_eflag = 1;
             else
                h.data.F_eflag = Start;
             end
        else
             disp('No flow measures for this experiment!');
             h.flags.IsThereFlow = false;
        end
        if( exist(h.paths.Fluo, 'file') )
             h.flags.IsThereFluo = true;
             
             h.data.fInfo = matfile(h.paths.Fluo);
             h.data.AcqFreq = h.data.fInfo.Freq;
             h.data.fDatPtr = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
             nframes = h.data.fInfo.datLength;
             Ts = min(nframes, Ts);
             
             stim = h.data.fInfo.Stim;
             if( ~h.flags.MasterStim )
                 h.data.MasterStim.StimLength = length(stim)/h.data.fInfo.Freq;
             end
             if( size(stim,2) > size(stim,1) )
                 stim = stim';
             end
            idxS = floor(h.data.AcqFreq*h.data.MasterStim.PreStimLength);
             if( idxS < 1 )
                 idxS = 1;
             end
             Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
            
             if( isempty(Start) )
                 h.data.F_eflag = 1;
             else
                h.data.F_eflag = Start;
             end
             
        else
             disp('No fluorescence measures for this experiment!');
             h.flags.IsThereFluo = false;
        end
        
        if( exist(h.paths.HbFile, 'file') )
            h.flags.IsThereHbO = true;
            h.flags.IsThereHbR = true;
            h.flags.IsThereHbT = true;
            
            h.data.HBinfos = matfile(h.paths.HbFile);
            h.data.AcqFreq = h.data.HBinfos.Freq;
            h.data.hoDatPtr = memmapfile([h.paths.FolderName filesep 'HbO.dat'], 'Format', 'single');
            h.data.hrDatPtr = memmapfile([h.paths.FolderName filesep 'HbR.dat'], 'Format', 'single');
            nframes = h.data.HBinfos.datLength;
            Ts = min(nframes, Ts);
             
            stim = h.data.HBinfos.Stim;
            if( ~h.flags.MasterStim )
                h.data.MasterStim.StimLength = length(stim)/h.data.HBinfos.Freq;
            end
            if( size(stim,2) > size(stim,1) )
                stim = stim';
            end
            idxS = floor(h.data.AcqFreq*h.data.MasterStim.PreStimLength);
             if( idxS < 1 )
                 idxS = 1;
             end
             Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
             if( isempty(Start) )
                 h.data.H_eflag = 1;
             else
                h.data.H_eflag = Start;
             end
        else            
            disp('No Hb concentrations were computed for this experiment!');
            h.flags.IsThereHbO = false;
            h.flags.IsThereHbR = false;
            h.flags.IsThereHbT = false;
        end
  
        Map = [];
        RawDatFiles = dir([h.paths.FolderName filesep 'Data_*.mat']);
        if( isempty(RawDatFiles) )
            %No compatible files were found
            h.flags.IsThereHbO = false;
            h.flags.IsThereHbR = false;
            h.flags.IsThereHbT = false;
            h.flags.MasterStim = false;
            
            disp(['No data files found in ' FolderName ' Folder.']);
            disp('There is nothing to show you... sorry!');
            return;
        end
        %Green channel:
        if( ~isempty(strfind([RawDatFiles.name],'green')) ) %#ok<*STREMP>
            h.flags.IsThereGreen = true;
            Dat_Gptr = matfile([h.paths.FolderName filesep 'Data_green.mat']);
            nrows = Dat_Gptr.datSize(1,2);
            ncols = Dat_Gptr.datSize(1,1);
            nframes = Dat_Gptr.datLength;
            Freq =  Dat_Gptr.Freq;
            h.data.AcqFreq = Freq;
            Ws = ncols;
            Hs = nrows;
            Ts = min(nframes, Ts);
            stim = Dat_Gptr.Stim;
            if( ~h.flags.MasterStim )
                h.data.MasterStim.StimLength = length(stim)/Dat_Gptr.Freq;
            end
            if( size(stim,2) > size(stim,1) )
                stim = stim';
            end
            idxS = floor(h.data.AcqFreq*h.data.MasterStim.PreStimLength);
             if( idxS < 1 )
                 idxS = 1;
             end
             Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
            
            if( isempty(Start) )
                 h.data.G_eflag = 1;
             else
                h.data.G_eflag = Start;
             end
            h.data.gDatPtr = memmapfile([h.paths.FolderName filesep 'gChan.dat'],...
                'Format', 'single');
            if( isempty(Map) )
                Map = reshape(h.data.gDatPtr.Data(1:(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        
        %Yellow channel:
        if( ~isempty(strfind([RawDatFiles.name],'yellow')) )
            h.flags.IsThereYellow = true;
            Dat_Yptr = matfile([h.paths.FolderName filesep 'Data_yellow.mat']);
            nrows = Dat_Yptr.datSize(1,2);
            ncols = Dat_Yptr.datSize(1,1);
            nframes = Dat_Yptr.datLength;
            Freq =  Dat_Yptr.Freq;
            h.data.AcqFreq = Freq;
            Ws = ncols;
            Hs = nrows;
            Ts = min(Ts, nframes);
            stim = Dat_Yptr.Stim;
            if( ~h.flags.MasterStim )
                h.data.MasterStim.StimLength = length(stim)/Dat_Yptr.Freq;
            end
            if( size(stim,2) > size(stim,1) )
                stim = stim';
            end
            idxS = floor(h.data.AcqFreq*h.data.MasterStim.PreStimLength);
            if( idxS < 1 )
                idxS = 1;
            end
            Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
            h.data.yDatPtr = memmapfile([h.paths.FolderName filesep 'yChan.dat'],...
                'Format', 'single');
            
              if( isempty(Start) )
                 h.data.Y_eflag = 1;
             else
                h.data.Y_eflag = Start;
             end
            if( isempty(Map) )
                Map = reshape(h.data.yDatPtr.Data(1:(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        
        %Red channel:
        if( ~isempty(strfind([RawDatFiles.name],'red')) )
            h.flags.IsThereRed = true;
            Dat_Rptr = matfile([h.paths.FolderName filesep 'Data_red.mat']);
            nrows = Dat_Rptr.datSize(1,2);
            ncols = Dat_Rptr.datSize(1,1);
            nframes = Dat_Rptr.datLength;
            Freq =  Dat_Rptr.Freq;
            h.data.AcqFreq = Freq;
            Ws = ncols;
            Hs = nrows;
            Ts = min(Ts, nframes);
            stim = Dat_Rptr.Stim;
            if( ~h.flags.MasterStim )
                h.data.MasterStim.StimLength = length(stim)/Dat_Rptr.Freq;
            end
            if( size(stim,2) > size(stim,1) )
                stim = stim';
            end
            idxS = floor(h.data.AcqFreq*h.data.MasterStim.PreStimLength);
            if( idxS < 1 )
                idxS = 1;
            end
            Start = find(diff(stim(idxS:end,1),1,1) > 0) + 1;
           
               if( isempty(Start) )
                 h.data.R_eflag = 1;
             else
                h.data.R_eflag = Start;
             end
            h.data.rDatPtr = memmapfile([h.paths.FolderName filesep 'rChan.dat'],...
                'Format', 'single');
            if( isempty(Map) )
                Map = reshape(h.data.rDatPtr.Data(1:(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        h.data.NCols = Ws;
        h.data.NRows = Hs;
        h.data.NFrames = Ts;
        
        if( round(h.data.AcqFreq*h.data.MasterStim.StimLength) > Ts )
            h.data.MasterStim.StimLength = Ts/h.data.AcqFreq;
        end
        
        %Map
        h.data.Map =  double(Map)./max(double(Map(:)));
        imshow(Map,[],'Parent',h.ui.ROIsMap);
        axis(h.ui.ROIsMap,'off', 'image');
        
        %ROIs file:
        ROIs = {};
        if( exist(h.paths.ROIsFile, 'file') )
            load(h.paths.ROIsFile);
            h.data.ROIs = ROIs;
            clear ROIs;
        else
            h.data.ROIs = ROIs;
        end
        
        %Events file:
        E = ones(1, h.data.MasterStim.NbStim);
        if(  exist(h.paths.EVNTsFile, 'file')  )
            load([h.paths.FolderName filesep 'Events.mat']);
        end
        h.data.EvntList = E;
        h.ui.EventsDispPan.Slider.Max = h.data.MasterStim.NbStim;
        clear E;
        
        Str = {};
        if( h.flags.IsThereGreen )
            Str{end+1} = 'Green';
        end
        if( h.flags.IsThereRed )
            Str{end+1} = 'Red';          
        end
        if( h.flags.IsThereYellow )    
            Str{end+1} = 'Yellow'; 
        end
        if( h.flags.IsThereHbO )    
            Str{end+1} = 'HbO';        
        end
        if( h.flags.IsThereHbR )    
            Str{end+1} = 'HbR';         
        end
        if( h.flags.IsThereHbT )
            Str{end+1} = 'HbT';         
        end
        if( h.flags.IsThereFlow )
            Str{end+1} = 'Flow';         
        end
        if( h.flags.IsThereFluo )
            Str{end+1} = 'Fluo';         
        end
        speckle_path = [h.paths.FolderName filesep 'Data_speckle.mat'];
        if(isfile(speckle_path))
             h.flags.IsThereSpeckle = true;
             Str{end+1} = 'Speckle';
        end
        if(h.flags.SlaveStim)
            idxslave = find(h.data.SlaveStim.Stim);
            idxmaster = find(h.data.MasterStim.Stim);
            h.data.ssdiff = round((idxslave(1)-idxmaster(1))*(h.data.SlaveStim.StimLength*h.data.SlaveStim.NbStim)/length(idxslave));
        end
        set(h.ui.ChannelSelector,'String', Str);
        
        set(h.ui.AddButton, 'Enable', 'on');
        set(h.ui.RemButton, 'Enable', 'on');
        set(h.ui.SaveROIpb, 'Enable', 'on');
        set(h.ui.LoadROIpb, 'Enable', 'on');
        set(h.ui.EventsDispPan.Slider, 'Enable', 'on');
        set(h.ui.ROIsSelector, 'Enable', 'on');
        set(h.ui.ChannelSelector, 'Enable', 'on');
        set(h.ui.FilteringOpt, 'Enable', 'on');
        set(h.ui.GlobalOpt, 'Enable', 'on');
        set(h.ui.StartVideo, 'Enable', 'on');
        set(h.ui.dispString,'string',h.paths.FolderName);

        RefreshLoop('All');
        
        delete(opendlg);
    end

    function StartVideo(~,~,~)
        
        if( ~h.flags.VideoPlaying )
            h.ui.VScreen = figure('Name', 'Video', 'Position', [200 200 500 500], 'Visible', 'off');
            colormap(h.ui.VScreen, 'jet');
            h.ui.ScreenAx = axes('Parent', h.ui.VScreen);
            
            Str = {};
            if( h.flags.IsThereGreen )    
                Str{end+1} = 'Green';   
            end
            if( h.flags.IsThereRed )    
                Str{end+1} = 'Red';          
            end            
            if( h.flags.IsThereYellow )
                Str{end+1} = 'Yellow'; 
            end
            if( h.flags.IsThereHbO )    
                Str{end+1} = 'HbO';        
            end
            if( h.flags.IsThereHbR )    
                Str{end+1} = 'HbR';         
            end
            if( h.flags.IsThereHbR && h.flags.IsThereHbO)  
                Str{end+1} = 'HbT'; 
            end
            if( h.flags.IsThereFlow )  
                Str{end+1} = 'Flow'; 
            end
            if( h.flags.IsThereFluo )
                Str{end+1} = 'Fluo';         
            end
            if( h.flags.IsThereSpeckle )
                Str{end+1} = 'Speckle';         
            end
            
            [selchan, valid] = listdlg('PromptString', 'Select channel:',...
                'SelectionMode', 'single',...
                'ListString', Str);
            if( ~valid )
                return;
            end
            SelectedSrc = Str{selchan};
            
            if( size(h.data.ROIs,2) > 0 )
                Str = cell(size(h.data.ROIs,2)+1,1);
                Str{1} = 'All';
                for indR = 1:size(h.data.ROIs,2)
                    Str{indR+1} = h.data.ROIs{indR}.name;
                end
            else
                Str = 'All';
            end
            [selroi, valid] = listdlg('PromptString', 'Select Roi:',...
                'SelectionMode', 'single',...
                'ListString', Str);
            if( ~valid )
                return;
            end
                      
            %Waiting Dlg...
            opendlg = dialog('Position',[500 500 250 150],'Name','Loading...');
            uicontrol('Parent', opendlg, 'Style','text',...
                'Position',[20 80 210 40], 'String', 'Loading data. Please wait...');
            pause(0.1);
            
            data = 0;
           
            if( isempty(strfind(SelectedSrc, 'Hb')) && isempty(strfind(SelectedSrc, 'Flow')) )
                isHb = 0;
                eval(['StartPts = h.data.' SelectedSrc(1) '_eflag;']);
                eval(['data = h.data.' lower(SelectedSrc(1)) 'DatPtr;']);
                h.data.vidClim = [0.99 1.01];
            elseif( length(SelectedSrc) > 3 )
                isHb = -1;
                StartPts = h.data.F_eflag;
                data = h.data.fDatPtr;
                h.data.vidClim = [0.75 1.25];
            elseif( SelectedSrc == 'HbT' )
                isHb = 2;
                StartPts = h.data.H_eflag;
                dO = h.data.hoDatPtr;
                dR = h.data.hrDatPtr;
                h.data.vidClim = [-5 5];
            else
                isHb = 1;
                StartPts = h.data.H_eflag;
                eval(['data = h.data.h' lower(SelectedSrc(3)) 'DatPtr;']);
                h.data.vidClim = [-5 5];
            end
            prompt = {'Colormap minimum:','Colormap maximum:'};
            dlg_title = 'Colormap Limits';
            num_lines = 1;
            defaultans = {num2str(h.data.vidClim(1)),num2str(h.data.vidClim(2))};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            h.data.vidClim(1) = str2double(answer{1});
            h.data.vidClim(2) = str2double(answer{2});
            
            if(selroi == 1)
                h.data.vidMap = ones(size(h.data.Map));
            else
                h.data.vidMap =  h.data.ROIs{selroi-1}.mask;
            end
            eLen = floor(h.data.AcqFreq*...
                (h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min));
            Accum = zeros(size(h.data.Map,1), size(h.data.Map,2), eLen);
            T = linspace(-h.data.MasterStim.PreStimLength, h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min - h.data.MasterStim.PreStimLength, eLen);
            nb_stim = 0;
            for indE = find(h.data.EvntList)
                if( isHb < 2 )
                    d = data.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                else
                     d = dO.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                     d = d + dR.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                end
                d = reshape(d, size(Accum));
                
                Pstart = median(d(:, :, 1:floor(5*h.data.AcqFreq)),3);
                Pend = median(d(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                if( isHb <= 0 )
                    d = d./L;
                else
                    d = d - L;
                end                
                Accum = Accum + d;
                nb_stim = nb_stim + 1;
            end
            Accum = Accum./nb_stim;
            Accum = imfilter(Accum, fspecial('gaussian',5,3),'same','symmetric');
            h.data.Accumulator = Accum;
            h.data.vidInd = 1;
            h.data.vidLimit = eLen;
            h.data.vidTimeVect = T;
            h.data.vidChan = SelectedSrc;
            set(h.ui.StartVideo, 'String', 'Stop');
            set(h.ui.VScreen, 'Visible', 'on');
            %timer...
            h.data.Timer = timer;
            h.data.Timer.TimerFcn = @ShowNextFrame;
            h.data.Timer.Period = 0.2;
            h.data.Timer.ExecutionMode = 'fixedRate';
            start(h.data.Timer);
            h.flags.VideoPlaying = true;
            delete(opendlg);
        else
            if( isvalid( h.data.Timer ) )
                stop(h.data.Timer);
                delete(h.data.Timer);
            end
            set(h.ui.StartVideo, 'String', 'Start');
            if( isvalid(h.ui.VScreen) )
                delete(h.ui.VScreen);
            end
            h.flags.VideoPlaying = false;
        end
    end

    function ShowNextFrame(~,~,~)
        
        if( ~isvalid(h.ui.VScreen) )
            stop(h.data.Timer);
            delete(h.data.Timer);
            set(h.ui.StartVideo, 'String', 'Start');
            h.flags.VideoPlaying = false;
            return;
        end
        figure(h.ui.VScreen);
        if( h.ui.Overlay )
            image(h.ui.ScreenAx, repmat(h.data.Map,1,1,3));
            hold on;
            i = imagesc(h.ui.ScreenAx, squeeze(h.data.Accumulator(:, :, h.data.vidInd)));
            hold off;
            title([h.data.vidChan ' channel at:' num2str(h.data.vidTimeVect(h.data.vidInd)) 's']);
            colorbar;
            
            Center = mean(h.data.vidClim);
            tMax = h.data.vidClim(2) - Center;
            tMin = h.data.vidClim(1) - Center;
            aMap = squeeze(h.data.Accumulator(:, :, h.data.vidInd));
            %  aMap = imfilter(aMap, fspecial('gaussian',32,8),'same','symmetric');
            aMap((aMap > Center + 0.25*tMin) & (aMap < Center + 0.25*tMax)) = 0;
            aMap(aMap > 0 ) = 1;
            
            
            alpha(i, aMap.*h.data.vidMap);
            set(h.ui.ScreenAx, 'CLim', h.data.vidClim);
            axis(h.ui.ScreenAx, 'image');
            h.data.vidInd = h.data.vidInd + 1;
            if( h.data.vidInd > h.data.vidLimit )
                h.data.vidInd = 1;
            end
        else
            i = imagesc(h.ui.ScreenAx, squeeze(h.data.Accumulator(:, :, h.data.vidInd)));
            set(h.ui.ScreenAx, 'CLim', h.data.vidClim);
            axis(h.ui.ScreenAx, 'image');
            title([h.data.vidChan ' channel at:' num2str(h.data.vidTimeVect(h.data.vidInd)) 's']);
            colorbar;
            h.data.vidInd = h.data.vidInd + 1;
            if( h.data.vidInd > h.data.vidLimit )
                h.data.vidInd = 1;
            end
        end
        
    end

    function GlobalSigChange(~, ~, ~)
        set(h.ui.EventsDispRegen, 'Enable', 'on'); 
    end

    function SrcChange(~,~,~)
        set(h.ui.EventsDispRegen, 'Enable', 'on'); 
    end

    function FiltOptChange(~, ~, ~)
        set(h.ui.EventsDispRegen, 'Enable', 'on'); 
    end

    function ret = FilterData(data, type)
        if( strcmp(type, 'IOI') )
            f = fdesign.lowpass('N,F3dB', 4, 0.4, h.data.AcqFreq);
            hpass = design(f,'butter');
            f = fdesign.lowpass('N,F3dB', 4, 1/60, h.data.AcqFreq);
            lpass = design(f,'butter');
            
            Hd= filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(data));
            Ld = filtfilt(lpass.sosMatrix,lpass.ScaleValues, double(data));
            ret = single(Hd./Ld);
        elseif( strcmp(type, 'Fluo') )
            
        end
    end

    function MoveEventsDisp(~,~,~)
        
        idx = round(get(h.ui.EventsDispPan.Slider, 'Value'));
        T = h.data.EventBuf(1, :);
        d = h.data.EventBuf(idx+1, :);
        h.ui.EventsDispPan.dispAx.mainline = plot(h.ui.EventsDispPan.Ax, T, d, 'k');

        if(h.flags.SlaveStim)
            if(h.data.ssdiff >= 0)
                h.ui.EventsDispPan.dispAx.zeroline = line(h.ui.EventsDispPan.Ax, [0 0], [min(d) max(d)],...
                    'Color', 'g', 'LineStyle','--');
                h.ui.EventsDispPan.dispAx.Masterstimline= line(h.ui.EventsDispPan.Ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength],...
                    [min(d) max(d)],'Color', 'r', 'LineStyle','--');
                h.ui.EventsDispPan.dispAx.Slavestimline= line(h.ui.EventsDispPan.Ax, [h.data.ssdiff+h.data.SlaveStim.StimLength h.data.ssdiff+h.data.SlaveStim.StimLength],...
                    [min(d) max(d)],'Color', 'm', 'LineStyle',':');
                 h.ui.EventsDispPan.dispAx.Slavezeroline = line(h.ui.EventsDispPan.Ax, [h.data.ssdiff h.data.ssdiff],...
                   [min(d) max(d)],'Color', 'b', 'LineStyle',':');
            else
                h.ui.EventsDispPan.dispAx.zeroline = line(h.ui.EventsDispPan.Ax, [abs(h.data.ssdiff) abs(h.data.ssdiff)],...
                    [min(d) max(d)],'Color', 'g', 'LineStyle','--');
                h.ui.EventsDispPan.dispAx.Masterstimline= line(h.ui.EventsDispPan.Ax,...
                    [(h.data.MasterStim.StimLength+abs(h.data.ssdiff)) (h.data.MasterStim.StimLength+abs(h.data.ssdiff))],...
                    [min(d) max(d)],'Color', 'r', 'LineStyle','--');
                h.ui.EventsDispPan.dispAx.Slavestimline= line(h.ui.EventsDispPan.Ax, ...
                    [h.data.SlaveStim.StimLength h.data.SlaveStim.StimLength],...
                    [min(d) max(d)],'Color', 'm', 'LineStyle',':');
                 h.ui.EventsDispPan.dispAx.Slavezeroline = line(h.ui.EventsDispPan.Ax, [0 0],...
                   [min(d) max(d)],'Color', 'b', 'LineStyle',':');
            end
        else
                h.ui.EventsDispPan.dispAx.zeroline = line(h.ui.EventsDispPan.Ax, [0 0], [min(d) max(d)],...
                    'Color', 'g', 'LineStyle','--');
                h.ui.EventsDispPan.dispAx.Masterstimline= line(h.ui.EventsDispPan.Ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength],...
                    [min(d) max(d)],'Color', 'r', 'LineStyle','--');
        end
        if( mean(d) > 0.5 )
            h.ui.EventsDispPan.dispAx.baseline = line(h.ui.EventsDispPan.Ax, [T(1) T(end)], [1 1],...
                'Color', 'k', 'LineStyle',':');
        else
            h.ui.EventsDispPan.dispAx.baseline = line(h.ui.EventsDispPan.Ax, [T(1) T(end)], [0 0],...
                'Color', 'k', 'LineStyle',':');
        end
        xlim(h.ui.EventsDispPan.Ax, [T(1), T(end)]);
        
        %Creer le checkbox associe
        set(h.ui.EventsDispPan.Cbox,'String',  int2str(idx) );        
        set(h.ui.EventsDispPan.Cbox,'Value',  h.data.EvntList(idx) );
        
        % hide events if Show events unchecked
        if(h.ui.EventDispCB.Value == 0)
            set(h.ui.EventsDispPan.Ax,'visible','off');
            set(h.ui.EventsDispPan.dispAx.mainline,'visible','off');
            set(h.ui.EventsDispPan.dispAx.zeroline,'visible','off');
            set(h.ui.EventsDispPan.dispAx.Masterstimline,'visible','off');
            set(h.ui.EventsDispPan.dispAx.baseline,'visible','off');
        end
    end

    function PopulateEvntsDisplay(~, ~, ~)
         if( isempty(h.data.ROIs) )
             return;
         end
         tc_path = [h.paths.FolderName filesep 'Events_timecourse' filesep];
         if(~exist(tc_path,'dir'))
             mkdir(tc_path);
         end
              
        sID = get(h.ui.ChannelSelector, 'Value');
        sStr = get(h.ui.ChannelSelector, 'String');
        SelectedSrc = sStr{sID};
        
        rID = get(h.ui.ROIsSelector, 'Value');
        rStr = get(h.ui.ROIsSelector, 'String');
        SelectedROI = rStr{rID};
           
        if( ValidateEvntSrc(SelectedSrc, SelectedROI) )
            
            waitdialog = dialog('Position',[500 500 250 150],'Name','Computing...');
            uicontrol('Parent', waitdialog, 'Style','text',...
                'Position',[20 80 210 40], 'String', 'Computing Events curve. Please wait...');
            pause(0.1);
            
            data = 0;
            StartPts = 1;
            isHbT = 0; isHb = 0;
            if( isempty(strfind(SelectedSrc, 'Hb')) && isempty(strfind(SelectedSrc, 'Flow')) )
                eval(['StartPts = h.data.' SelectedSrc(1) '_eflag;']);
                eval(['data = h.data.' lower(SelectedSrc(1)) 'DatPtr;']);
            elseif( ~isempty(strfind(SelectedSrc, 'Flow')) )
                StartPts = h.data.F_eflag;
                isHb = -1;
                data = h.data.fDatPtr;
            else
                StartPts = h.data.H_eflag;
                isHb = 1;
                if( SelectedSrc == 'HbR' )
                    eval('data = h.data.hrDatPtr;');
                elseif( SelectedSrc == 'HbO' )
                    eval('data = h.data.hoDatPtr;');
                else
                    isHbT = 1;
                end
            end
            eLen = floor(h.data.AcqFreq*...
                (h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min));
            
            if( strcmp(SelectedROI, 'AllPixels') )
                mask = ones(size(h.data.Map));
            else
                idx = arrayfun(@(a) strcmp(h.data.ROIs{a}.name, SelectedROI), 1:size(h.data.ROIs,2));
                mask = h.data.ROIs{idx == 1}.mask;                                                                       
            end
            
            if( h.data.MasterStim.NbStim == 1 )
                h.ui.EventsDispPan.Visible = false;
            else
                set(h.ui.EventsDispPan.Slider,'Max',h.data.MasterStim.NbStim,'SliderStep',[1/h.data.MasterStim.NbStim 1/h.data.MasterStim.NbStim]);
                h.ui.EventsDispPan.Visible = true;
                h.ui.EventsDispPan.min = 1 - (h.data.MasterStim.NbStim)*0.6;
                if( h.ui.EventsDispPan.min > 0 )
                    h.ui.EventsDispPan.min = 0;
                end
            end
            
            T = linspace(-h.data.MasterStim.PreStimLength, h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min - h.data.MasterStim.PreStimLength, eLen);
            h.data.EventBuf = zeros(h.data.MasterStim.NbStim + 1, eLen, 'single');
            h.data.EventBuf(1,:) = T;
            for indE = 1:h.data.MasterStim.NbStim
                
                if( isHbT == 0 )
                    d = data.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                else
                    d1 = h.data.hoDatPtr.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                    d2 = h.data.hrDatPtr.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                    d = d1+d2;
                end
                
                d_im = reshape(d,size(h.data.Map,1),size(h.data.Map,2),[]);
                
                if(strcmp(SelectedSrc,'Red') || strcmp(SelectedSrc,'Green') ||strcmp(SelectedSrc,'Yellow')||strcmp(SelectedSrc,'Flow') )
                    Pstart = median(d_im(:, :, 1:floor(5*h.data.AcqFreq)),3);
                    Pend = median(d_im(:,:,(end-floor(5*h.data.AcqFreq)):end),3);
                    m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                    L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                        (m*T(round(end - h.data.MasterStim.PreStimLength/2))));
                    d_im = d_im./L;
                end
                d_im = imfilter(d_im, fspecial('gaussian',5,3),'same','symmetric');
                if(get(h.ui.EventSaveCB,'Value'))
                    fig_name = char(strcat(sStr{sID}, '-', rStr{rID} ,' Event' ,num2str(indE)));
                    fig = figure('visible','off','units','normalized','outerposition',[0 0 1 1]);
                    nb_pict = 0;
                    for indT = 1:eLen
                        if(abs(mod(T(indT),2)) < 0.1)
                            nb_pict = nb_pict + 1;
                            subplot(1,12,nb_pict);
                            hold on
                            if(strcmp(SelectedSrc,'Red') || strcmp(SelectedSrc,'Green') ||strcmp(SelectedSrc,'Yellow'))
                                imagesc(d_im(:,:,indT),[0.99,1.01]);
                            elseif(strcmp(SelectedSrc,'Flow'))
                                imagesc(d_im(:,:,indT),[0.75 1.25]);
                            else
                                imagesc(d_im(:,:,indT),[-5,5])
                            end
                            axis off;
                            colormap jet;
                            title_string = char(strcat(num2str(round(T(indT)))));  
                            title(title_string);
                        end
                        if(nb_pict == 12)
                            break;
                        end
                    end
                    ha=get(gcf,'children');
                    set(ha(12),'position',[0.07 .5 .07 .1]);
                    set(ha(11),'position',[0.14 .5 .07 .1]);
                    set(ha(10),'position',[0.21 .5 .07 .1]);
                    set(ha(9),'position',[0.28 .5 .07 .1]);
                    set(ha(8),'position',[0.35 .5 .07 .1]);
                    set(ha(7),'position',[0.42 .5 .07 .1]);
                    set(ha(6),'position',[0.49 .5 .07 .1]);
                    set(ha(5),'position',[0.56 .5 .07 .1]);
                    set(ha(4),'position',[0.63 .5 .07 .1]);
                    set(ha(3),'position',[0.70 .5 .07 .1]);
                    set(ha(2),'position',[0.77 .5 .07 .1]);
                    set(ha(1),'position',[.84 .5 .09 .1]);
                    hold off
                    hp4 = get(ha(1),'Position');
                    colorbar('Position', [hp4(1)+hp4(3)+0.02  0.5  0.01  hp4(4)]);
                    file_path = [tc_path fig_name];
                    print(fig,file_path,'-djpeg');
                    delete(fig);
                end
                
                d = reshape(d, [], eLen);
                
                Glob = mean(d, 1);
                d = mean(d(mask(:)==1,:),1);
                
                if( (isHb == 0) && get(h.ui.GlobalOpt, 'Value') )
                    d = d./Glob;
                end
                if( (isHb == 0) && get(h.ui.FilteringOpt, 'Value') )
                    d = FilterData( d, 'IOI');
                end
                
                %Detrend
                Pstart = median(d(1:floor(5*h.data.AcqFreq)));
                Pend = median(d((end-floor(5*h.data.AcqFreq)):end));
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.MasterStim.PreStimLength));
                L = m*T + (Pend - m*T(round(end - h.data.MasterStim.PreStimLength/2)));
                if( isHb <= 0 )
                    d = d./L;
                else
                    d = d - L;
                end
                
                h.data.EventBuf(indE + 1, :) = d;

            end
            MoveEventsDisp();
            
        end
        delete(waitdialog);
        MeanRecalculation();
        set(h.ui.EventsDispRegen, 'Enable', 'off'); 
    end

    function OnEditEventsClicked(~, src, ~)
        V = get(src.Source,'Value');
        h.data.EvntList(str2double(get(src.Source, 'String'))) = V;
        MeanRecalculation();
        h.flags.saveEvnts = true;
    end

    function MeanRecalculation()
        d = zeros(1, floor(h.data.AcqFreq*...
                (h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min)));
        for indE = 1:length(h.data.EvntList)
            if(h.data.EvntList(indE))
                d = d + h.data.EventBuf(indE + 1,:);
            end
        end
        d = d/sum(h.data.EvntList);
        T = linspace(-h.data.MasterStim.PreStimLength, h.data.MasterStim.StimLength + h.data.MasterStim.InterStim_min - h.data.MasterStim.PreStimLength, length(d));
        h.ui.EventsMeanPan.dispAx.mainline = plot(h.ui.EventsMeanPan.Ax, T, d, 'k');

        if(h.flags.SlaveStim)
            if(h.data.ssdiff >= 0)
                h.ui.EventsMeanPan.dispAx.zeroline= line(h.ui.EventsMeanPan.Ax, [0 0],...
                 [min(d) max(d)],'Color', 'g', 'LineStyle','--');
                h.ui.EventsMeanPan.dispAx.Masterstimline = line(h.ui.EventsMeanPan.Ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength],...
                    [min(d) max(d)], 'Color', 'r', 'LineStyle','--');
                h.ui.EventsMeanPan.dispAx.Slavestimline = line(h.ui.EventsMeanPan.Ax, [(h.data.SlaveStim.StimLength+h.data.ssdiff) (h.data.SlaveStim.StimLength+h.data.ssdiff)],...
                    [min(d) max(d)], 'Color', 'm', 'LineStyle',':');
                 h.ui.EventsMeanPan.dispAx.Slavezeroline = line(h.ui.EventsMeanPan.Ax, [h.data.ssdiff h.data.ssdiff],...
                    [min(d) max(d)],'Color', 'b', 'LineStyle',':');
            else
                 h.ui.EventsMeanPan.dispAx.zeroline= line(h.ui.EventsMeanPan.Ax, [(abs(h.data.ssdiff)) (abs(h.data.ssdiff))],...
                    [min(d) max(d)],'Color', 'g', 'LineStyle','--');
                 h.ui.EventsMeanPan.dispAx.Masterstimline = line(h.ui.EventsMeanPan.Ax, [h.data.MasterStim.StimLength+abs(h.data.ssdiff) h.data.MasterStim.StimLength+abs(h.data.ssdiff)],...
                    [min(d) max(d)], 'Color', 'r', 'LineStyle','--');
                 h.ui.EventsMeanPan.dispAx.Slavestimline = line(h.ui.EventsMeanPan.Ax, [h.data.SlaveStim.StimLength h.data.SlaveStim.StimLength],...
                    [min(d) max(d)], 'Color', 'm', 'LineStyle',':');
                 h.ui.EventsMeanPan.dispAx.Slavezeroline = line(h.ui.EventsMeanPan.Ax, [0 0],...
                    [min(d) max(d)],'Color', 'b', 'LineStyle',':');
            end
        else
            h.ui.EventsMeanPan.dispAx.zeroline= line(h.ui.EventsMeanPan.Ax, [0 0],...
                [min(d) max(d)],'Color', 'g', 'LineStyle','--');
            h.ui.EventsMeanPan.dispAx.Masterstimline = line(h.ui.EventsMeanPan.Ax, [h.data.MasterStim.StimLength h.data.MasterStim.StimLength],...
                [min(d) max(d)], 'Color', 'r', 'LineStyle','--');   
        end
        
        sID = get(h.ui.ChannelSelector, 'Value');
        sStr = get(h.ui.ChannelSelector, 'String');
        SelectedSrc = sStr{sID};
        if( isempty(strfind(SelectedSrc, 'Hb')) )
        h.ui.EventsMeanPan.dispAx.baseline = line(h.ui.EventsMeanPan.Ax, [T(1) T(end)], [1 1],...
                'Color', 'k', 'LineStyle',':');
        else
        h.ui.EventsMeanPan.dispAx.baseline = line(h.ui.EventsMeanPan.Ax, [T(1) T(end)], [0 0],...
                'Color', 'k', 'LineStyle',':');
        end
        xlim(h.ui.EventsMeanPan.Ax,[T(1), T(end)]);
        
        % hide events if Show events unchecked
        if(h.ui.EventDispCB.Value == 0)
            set(h.ui.EventsMeanPan.Ax,'visible','off');
            set(h.ui.EventsMeanPan.dispAx.mainline,'visible','off');
            set(h.ui.EventsMeanPan.dispAx.zeroline,'visible','off');
            set(h.ui.EventsMeanPan.dispAx.Masterstimline,'visible','off');
            set(h.ui.EventsMeanPan.dispAx.baseline,'visible','off');
        end
    end

    function bRet = ValidateEvntSrc(cSrc, rSrc)
        if( cSrc(1) == 'G' )
            bRet = h.flags.IsThereGreen;
        elseif( cSrc(1) == 'R' )
            bRet = h.flags.IsThereRed;
        elseif( cSrc(1) == 'Y' )
            bRet = h.flags.IsThereYellow;
        elseif( cSrc(1) == 'F' )
            bRet = h.flags.IsThereFlow || h.flags.IsThereFluo;
        elseif( cSrc(3) == 'O' )
            bRet = h.flags.IsThereHbO;
        elseif( cSrc(3) == 'R' )
            bRet = h.flags.IsThereHbR;
        elseif( cSrc(3) == 'T' )
            bRet = h.flags.IsThereHbT;
        else
            bRet = false;
        end
        
        if( ~strcmp(rSrc, 'AllPixels') && ~any(arrayfun(@(r) strcmp(rSrc, h.data.ROIs{r}.name), 1:size(h.data.ROIs,2))) )
            bRet = false;
        end
    end

    function RefreshLoop(Option)
        if( strcmp(Option, 'ROIs') )
            RefreshROIsList();
            RefreshMapDisplay();
        elseif( strcmp(Option, 'Events') )
            if( h.data.MasterStim.NbStim > 1 )
                PopulateEvntsDisplay([],[],[]);
            end
        elseif( strcmp(Option, 'All') )
            RefreshROIsList();
            RefreshMapDisplay();
            if( h.data.MasterStim.NbStim > 1 )
                PopulateEvntsDisplay([],[],[]);
            end
        end
    end

    function RefreshMapDisplay(~,~,~)
        lst = get(h.ui.ROIsMap, 'Children');
        while( ~isempty(lst) )
           delete(lst(1));
           lst = get(h.ui.ROIsMap, 'Children');
        end
        
        if( isfield(h.ui, 'MskAxes') )
            for indA = 1:size(h.ui.MskAxes,2)
                cla(h.ui.MskAxes(indA));
            end
            h.ui = rmfield(h.ui, 'MskAxes');
        end
        
        cmap = gray(64);
        image(h.ui.ROIsMap, round(h.data.Map*64), 'AlphaData', 1);
        colormap(h.ui.ROIsMap, cmap);
        axis(h.ui.ROIsMap, 'off');
        
        NbRois = size(h.data.ROIs,2);
        for indR = 1:NbRois
            m = h.data.ROIs{indR}.mask;
            h.ui.MskAxes(indR) = axes('Parent', h.ui.ROIs, 'Position',[0.23 0.05 0.75 0.92]);
            image(h.ui.MskAxes(indR), m, 'AlphaData', m.*0.25);
            colormap(h.ui.MskAxes(indR), cat(1,[0 0 0], h.data.ROIs{indR}.color));
            axis(h.ui.MskAxes(indR), 'off');
        end 
    end

    function AddNewROI(~,~,~)
        h.flags.saveROIS = true;
        
        lst = get(h.ui.ROIsMap, 'Children');
        while( ~isempty(lst) )
            delete(lst(1));
            lst = get(h.ui.ROIsMap, 'Children');
        end
           
        if( isfield(h.ui, 'MskAxes') )
            for indA = 1:size(h.ui.MskAxes,2)
                cla(h.ui.MskAxes(indA));
            end
            h.ui = rmfield(h.ui, 'MskAxes');
        end
        
        image(h.ui.ROIsMap, round(h.data.Map*64));
        colormap(h.ui.ROIsMap, gray(64));
        axis(h.ui.ROIsMap, 'off');
        prompt = {'Name:'};
        dlg_title = 'New ROI';
        num_lines = 1;
        Tmp = {};
        if( ~isempty(h.data.ROIs) )
            Tmp = cell(size(h.data.ROIs,2),1);
            for indR = 1:size(h.data.ROIs,2)
                Tmp{indR} = h.data.ROIs{indR}.name;
            end
        end
        NewROI_ID = 1;
        while( sum(ismember(Tmp,['ROI_' int2str(NewROI_ID)])) )
            NewROI_ID = NewROI_ID + 1;
        end
        def = {['ROI_' int2str(NewROI_ID)]};
        
        answer = inputdlg(prompt, dlg_title,num_lines,def);
        if( ~isempty(answer) )
             type_list = {'Circle', 'Polygon', 'Surround','Rectangle','Free Hand','Point'};
             [indx,tf] = listdlg('PromptString','Select how to select ROI',...
                           'SelectionMode','single',...
                           'ListString',type_list);
            % Type = questdlg('ROI selection method:','Method selection',...
            %   'Circle', 'Polygon', 'Surround', 'Circle');
            Type = type_list{indx};
            h_im = get(h.ui.ROIsMap,'Children');
            
            switch Type
                case 'Circle'
                    e = imellipse(h.ui.ROIsMap);
                    mask = createMask(e, h_im(end));
                    delete(e);
                    clear e;
                case 'Rectangle'
                    r = imrect(h.ui.ROIsMap);
                    mask = createMask(r, h_im(end));
                    delete(r);
                    clear r;
                case 'Free Hand'
                    f = imfreehand(h.ui.ROIsMap);
                    mask = createMask(f, h_im(end));
                    delete(f);
                    clear f;
                case 'Polygon'
                    p = impoly(h.ui.ROIsMap);
                    mask = createMask(p, h_im(end));
                    delete(p);
                    clear p;
                 case 'Surround'
                     [orig, valid] = listdlg('PromptString', 'From wich ROI?',...
                         'SelectionMode', 'single',...
                         'ListString', Tmp);
                     if( ~valid )
                         return;
                     end
                     width = inputdlg( 'Width of the surround area:', 'Width', [1, 50], {'25'}); 
                     width = str2double(width{1});
                     
                     m = imfill(h.data.ROIs{orig}.mask,'holes');
                     mask = imdilate(m, strel('disk',width)) & ~m;
                 case 'Point'
                        p = impoint(h.ui.ROIsMap);
                        point = getPosition(p);
                        answer_rad = inputdlg('Enter Pixel radius');
                        radius = str2num(answer_rad{1});
                        th = 0:pi/50:2*pi;
                        roi_x = radius * cos(th) + point(1) ;
                        roi_y = radius * sin(th) + point(2);
                        mask = poly2mask(roi_x,roi_y, size(h.data.Map,1),size(h.data.Map,2));
                        delete(p);
                        clear p;
                     
            end
            h.data.ROIs{end+1} = struct('name',answer, 'mask', mask, 'color',  [0 0 1]);           
        end
        
        RefreshLoop('ROIs');
    end

    function RemoveROI(~,~,~)
        
        if( isempty(h.data.ROIs) )
            return;
        end
       
        Tmp = cell(size(h.data.ROIs,2),1);
        for indR = 1:size(h.data.ROIs,2)
            Tmp{indR} = h.data.ROIs{indR}.name;
        end
       
        [sel, valid] = listdlg('PromptString', 'Select the ROI to be removed:',...
            'SelectionMode', 'single',...
            'ListString', Tmp);
                     
        if( ~valid )
            return;
        end
        
        h.data.ROIs(sel) = [];  
        h.flags.saveROIS = true;
        RefreshLoop('ROIs');    
    end
        
    function OnChangeColorClicked(~, src, ~)
        h.flags.saveROIS = true;
        OldC = get(src.Source,'BackgroundColor');
        NewC = uisetcolor(OldC);
        ID = get(src.Source, 'UserData');
        h.data.ROIs{ID}.color = NewC;
        set(src.Source,'BackgroundColor',NewC);
        RefreshLoop('ROIs');
    end

    function RefreshROIsList(~ ,~, ~)
        lst = get(h.ui.ROIsPanel, 'Children');
        while( ~isempty(lst) )
            delete(lst(1));
            lst = get(h.ui.ROIsPanel, 'Children');
        end
        
        if( ~isempty(h.data.ROIs) )
            Str = cell(size(h.data.ROIs,2) + 1,1);
            Str{1} = 'AllPixels';
            for indR = 1:size(h.data.ROIs,2)
                uicontrol('Style','text','Parent',h.ui.ROIsPanel,...
                    'Units', 'normalized', 'Position',[0.01 0.99-indR*0.05 0.75 0.05],...
                    'String',h.data.ROIs{indR}.name, 'HorizontalAlignment','left');
                h.ui.ROIsColor(indR) = uicontrol('Style','pushbutton',...
                    'Parent',h.ui.ROIsPanel,'Units', 'Normalized',...
                    'Position',[0.80 0.99-indR*0.05 0.175 0.05], ...
                    'BackGroundColor', h.data.ROIs{indR}.color,...
                    'Callback',@OnChangeColorClicked,...
                    'UserData', indR);
                Str{indR+1} = h.data.ROIs{indR}.name;
            end
            set( h.ui.ROIsSelector, 'Value', 1 );
            set( h.ui.ROIsSelector, 'String', Str );
            set(h.ui.EventsDispRegen,'Enable','on');
        end
    end

    function SaveROIs(~,~,~)
        if( h.flags.saveROIS )
            ROIs = h.data.ROIs; %#ok<*NASGU>
            save(h.paths.ROIsFile, 'ROIs');
        end
        h.flags.saveROIS = false;
        SaveEvnts();
    end

    function SaveEvnts(~,~,~)
        if( h.flags.saveEvnts )
            E = h.data.EvntList;
            save(h.paths.EVNTsFile, 'E');
        end
        h.flags.saveEvnts = false;
    end

    function LoadROIs(~,~,~)
%         [FileName,PathName,FilterIndex] = uigetfile('ROIs.mat');
%         if( any(FileName ~= 'ROIs.mat') || FilterIndex == 0 )
%             return;
%         end
%         Tmp = load([PathName FileName]);
%         if( any(size(Tmp.ROIs{1}.mask) ~= size(h.data.Map)) )
%             msgbox('ROIs file selected does not fit data dimension.','Load ROIs');
%             return;
%         end
%         h.data.ROIs = Tmp.ROIs;
%         h.flags.saveROIS = true;
%         RefreshLoop('ROIs');
    end
    
    function Temporal_Speckle(~,~,~)
        if( h.flags.IsThereSpeckle)
            Dat_Speckle = matfile([h.paths.FolderName filesep 'Data_speckle.mat']);
            nrows = Dat_Speckle.datSize(1,2);
            ncols = Dat_Speckle.datSize(1,1);
            nframes = Dat_Speckle.datLength;
            Freq =  Dat_Speckle.Freq;
            temp_data = memmapfile([h.paths.FolderName filesep 'sChan.dat'], 'Format', 'single');
            data = reshape(temp_data.Data,nrows,ncols,[]);
            
            answer = inputdlg('Enter number of frames used','Input');
            
            if(isnumeric(str2num(answer{1})) && str2num(answer{1}) > 0)
                speckle_frames = str2num(answer{1});
                
                % getting ROIs boundaries
                se = strel('disk',3,4);
                for indk = 1:size(h.data.ROIs,2)
                    er = imerode(h.data.ROIs{indk}.mask,se);
                    diff = h.data.ROIs{indk}.mask - er;
                    nb_pos = 0;
                    for i = 1:size(diff,1)
                        for j = 1:size(diff,2)
                            if(diff(i,j))
                                nb_pos = nb_pos +1;
                                pos(nb_pos,1,indk) = i;
                                pos(nb_pos,2,indk) = j;
                            end
                        end
                    end
                end
                
                disp('Writing speckle video');
                v = VideoWriter([h.paths.Graphs filesep 'temp_speckle.avi']);
                v.FrameRate = Freq;
                open(v);
                
                spec_fig = figure('visible','off','color','k');
                ax = axes('Parent',spec_fig);
                
                for ind = 1:size(data,3)-speckle_frames
                    speckle_mean = mean(data(:,:,ind:ind+speckle_frames),3);
                    speckle_std = std(data(:,:,ind:ind+speckle_frames),[],3);
                    temp_speckle = speckle_std./speckle_mean;
                    imagesc(ax,temp_speckle);
                    hold(ax,'on');
                    for i = 1:size(h.data.ROIs,2)
                        roi_mean(ind,i) = mean2(temp_speckle(h.data.ROIs{i}.mask));
                        temp_pos = squeeze(pos(:,:,i));
                        plot(ax,temp_pos(:,2),temp_pos(:,1),'g.','MarkerSize',1);
                    end
                    axis(ax, 'image', 'off');
                    set(ax, 'CLim', [0 0.05]);
                    colormap(ax,gray(256));
                    hold(ax,'off');
                    frame = getframe(spec_fig);
                    writeVideo(v,frame);
                end
                close(v);
                delete(spec_fig);
                
               if(exist('roi_mean','var'))
                   t = 0:1/Freq:(size(roi_mean,1)-1)/Freq;
                   fig_plot = figure('visible','off');
                   hold on;
                   for i = 1:size(h.data.ROIs,2)
                       plot(t',roi_mean(:,i),'DisplayName',h.data.ROIs{i}.name);
                   end
                   ylabel('Intensity average');
                   xlabel('Time(s)');
                   legend;
                   title('Speckle course for each ROIs')
                   print(fig_plot,[h.paths.FolderName filesep 'Graphs' filesep 'temp_speckle.jpg'],'-djpeg');
                   delete(fig_plot);
               end
               disp('Video writing done!');
            else
                disp('Must be a non-null positive number.');
            end        
        else
            disp('No Speckle datas detected');
        end
    end

    function ShowEvents(~,~,~)
        if(h.ui.EventDispCB.Value ~= 0)
            set(h.ui.EventsMeanPan.Ax,'visible','on');
            set(h.ui.EventsDispPan.Ax,'visible','on');
            if(isfield(h.ui.EventsDispPan,'dispAx'))
                set(h.ui.EventsDispPan.dispAx.mainline,'visible','on');
                set(h.ui.EventsDispPan.dispAx.zeroline,'visible','on');
                set(h.ui.EventsDispPan.dispAx.Masterstimline,'visible','on');
                if(h.flags.SlaveStim)
                    set(h.ui.EventsDispPan.dispAx.Slavestimline,'visible','on');
                    if(isfield(h.ui.EventsDispPan.dispAx,'Slavezeroline'))
                        set(h.ui.EventsDispPan.dispAx.Slavezeroline,'visible','on')
                    end
                end
                set(h.ui.EventsDispPan.dispAx.baseline,'visible','on');
            end
            if(isfield(h.ui.EventsMeanPan,'dispAx'))
                set(h.ui.EventsMeanPan.dispAx.mainline,'visible','on');
                set(h.ui.EventsMeanPan.dispAx.zeroline,'visible','on');
                set(h.ui.EventsMeanPan.dispAx.Masterstimline,'visible','on');
                if(h.flags.SlaveStim)
                    set(h.ui.EventsMeanPan.dispAx.Slavestimline,'visible','on');
                    if(isfield(h.ui.EventsMeanPan.dispAx,'Slavezeroline'))
                        set(h.ui.EventsMeanPan.dispAx.Slavezeroline,'visible','on')
                    end
                end
                set(h.ui.EventsMeanPan.dispAx.baseline,'visible','on');
            end
        else
            set(h.ui.EventsMeanPan.Ax,'visible','off');
            set(h.ui.EventsDispPan.Ax,'visible','off');
            if(isfield(h.ui.EventsDispPan,'dispAx'))
                set(h.ui.EventsDispPan.dispAx.mainline,'visible','off');
                set(h.ui.EventsDispPan.dispAx.zeroline,'visible','off');
                set(h.ui.EventsDispPan.dispAx.Masterstimline,'visible','off');
                if(h.flags.SlaveStim)
                    set(h.ui.EventsDispPan.dispAx.Slavestimline,'visible','off');
                    if(isfield(h.ui.EventsDispPan.dispAx,'Slavezeroline'))
                        set(h.ui.EventsDispPan.dispAx.Slavezeroline,'visible','off')
                    end
                end
                set(h.ui.EventsDispPan.dispAx.baseline,'visible','off');
            end
            if(isfield(h.ui.EventsMeanPan,'dispAx'))
                set(h.ui.EventsMeanPan.dispAx.mainline,'visible','off');
                set(h.ui.EventsMeanPan.dispAx.zeroline,'visible','off');
                set(h.ui.EventsMeanPan.dispAx.Masterstimline,'visible','off');
                if(h.flags.SlaveStim)
                    set(h.ui.EventsMeanPan.dispAx.Slavestimline,'visible','off');
                    if(isfield(h.ui.EventsMeanPan.dispAx,'Slavezeroline'))
                        set(h.ui.EventsMeanPan.dispAx.Slavezeroline,'visible','off')
                    end
                end
                set(h.ui.EventsMeanPan.dispAx.baseline,'visible','off');
            end
        end
    end

    function my_closereq(~,~)
        if( h.flags.saveROIS )
            selection = questdlg('Do you want to save modified ROIs list?',...
                'Before closing...',...
                'Yes','No','Yes');
            if( strcmp(selection, 'Yes') )
                SaveROIs();
                SaveEvnts();
            end
        end
        
        delete(gcf);
    end

end