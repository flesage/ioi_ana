function out = DataNavigator(varargin)

%%%%%%%
% Vars init
%%%%%%%
h.flags.saveROIS = false; %flag to know if any changes were made to ROIs
h.flags.saveEvnts = false;
h.flags.IsThereHbO = false;
h.flags.IsThereHbR = false;
h.flags.IsThereHbT = false;
h.flags.IsThereFlow = false;
h.flags.IsThereGreen = false;
h.flags.IsThereYellow = false;
h.flags.IsThereRed = false;
h.flags.VideoPlaying = false;
%Map
h.data.Map = [];
h.data.EventBuf = 0;

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
h.ui.EventsDispPan.Container = uipanel('Parent', h.ui.EventsPan, 'Title', 'Selection', 'FontSize', 12,...
              'Position', [0.0 0.5 1.0 0.4]);
h.ui.EventsDispPan.Ax = axes('Parent', h.ui.EventsDispPan.Container, 'Position', ...
                    [0.175 0.15 0.72 0.85]);
h.ui.EventsDispPan.Cbox = uicontrol('Style','checkbox','Parent', h.ui.EventsDispPan.Container,...
    'Units', 'normalized', 'Position', [0.01 0.4 0.1 0.2],...
    'String','1', 'Callback', @OnEditEventsClicked);
h.ui.EventsDispPan.Slider = uicontrol('Parent', h.ui.EventsDispPan.Container,...
             'Style', 'slider',  'Min', 1, 'Max', 10, 'Value', 1,...
             'Units', 'normalized', 'Position', [0.95 0.0 0.05 1], 'SliderStep', [0.1 0.3],...
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

%%% Data graphs:
h.ui.XLS = uipanel('Parent', h.ui.fig, 'Title','Spreadsheet','FontSize',12,...
             'Position',[.150 .01 .125 .105]);
h.ui.Eport = uicontrol('Style','pushbutton','Parent', h.ui.XLS,...
    'Units', 'normalized', 'Position',[0.1 0.1 0.8 0.8],...
    'String','Export', 'Callback', @exportXLS);



%%%%%%%%%%%%%%%
%Fonctions & Callbacks:
%%%%%%%%%%%%%%%
    function GenerateGraphs(~,~,~)
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
        figName = ['RawMap'];
        saveas(fig, [h.paths.Graphs figName], 'png');
        close(fig);
        
        eLen = 5*floor(h.data.Stim.InterStim_min);
        T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, eLen);
        %for each ROI
        for indR = 1:size(h.data.ROIs,2)
            mask = h.data.ROIs{indR}.mask;
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
            
            for indE = 1:length(h.data.EvntList)
                set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ', Colours...']);
                
                fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                ax = axes('Parent', fig);
                hold(ax,'on');
                maxi = 1.01;
                mini = 0.99;

                if( h.flags.IsThereGreen )
                    %Open
                    d =   h.data.gDatPtr.Data((length(h.data.Map(:))*(h.data.G_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.G_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:25));
                     Pend = median(d((end-24):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     d = d./L;
                
                     if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     %Plot
                     plot(ax, T, d, 'Color', [0.0 0.75 0.0], 'LineWidth', 2);
                end
                if( h.flags.IsThereYellow )
                    %Open
                    d =   h.data.yDatPtr.Data((length(h.data.Map(:))*(h.data.Y_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.Y_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:25));
                     Pend = median(d((end-24):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     d = d./L;
                
                      if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     
                     %Plot
                     plot(ax, T, d, 'Color', [0.75 0.75 0.0], 'LineWidth', 2);
                end
                if( h.flags.IsThereRed )
                    %Open
                    d =   h.data.rDatPtr.Data((length(h.data.Map(:))*(h.data.R_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.R_eflag(indE) +eLen - 1)) );
                    d = reshape(d, [], eLen);
                    d = mean(d(mask(:) == 1, :), 1);
                    
                    %Filter
                    d = FilterData(d, 'IOI');
                    
                     %Detrend
                     Pstart = median(d(1:25));
                     Pend = median(d((end-24):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     d = d./L;
                
                      if( max(d) > maxi) 
                         maxi = max(d);
                     end
                     if( min(d) < mini) 
                         mini = min(d);
                     end
                     
                     %Plot
                     plot(ax, T, d, 'Color', [0.75 0.0 0.0], 'LineWidth', 2);
                end
                box(ax,'on');
                set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 2);
                title(['{\Delta}Reflectance over ' h.data.ROIs{indR}.name  ', E#' int2str(indE)]);
                ylabel('Normalized Reflectance');
                xlabel('Time (sec)');
                xlim([T(1), T(end)]);
                ylim([mini, maxi]);
                line(ax, [T(1) T(end)], [1 1], 'Color', 'k', 'LineStyle',':');
                line(ax, [0 0], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
                line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [0.99 1.01], 'Color', 'k', 'LineStyle','--');
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
                    dF =   h.data.fDatPtr.Data((length(h.data.Map(:))*(h.data.F_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.F_eflag(indE) +eLen - 1)) );
                    dF = reshape(dF, [], eLen);
                    dF = mean(dF(mask(:) == 1, :), 1);
                    
                     %Detrend
                     Pstart = median(dF(1:25));
                     Pend = median(dF((end-24):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
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
                     line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [mini maxi], 'Color', 'k', 'LineStyle','--');
                     xlim([T(1), T(end)]);
                     ylim([mini, maxi]);
                     %Save figure for colours:
                     figName = ['Flow_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                     saveas(fig, [h.paths.Graphs figName], 'png');
                     close(fig);
                end
                
                %Figure for Hbs:
                if( h.flags.IsThereHbO && h.flags.IsThereHbR )
                    set(GraphsStr, 'String', ['ROI #' int2str(indR) ', Event #' int2str(indE) ' Hbs...']);
                    fig = figure('InvertHardcopy','off','Color',[1 1 1], 'Visible', 'off');
                    ax = axes('Parent', fig);
                    hold(ax,'on');
                    maxi = 25;
                    mini = -25;
                    %Open
                    dO =   h.data.hoDatPtr.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                    dR =   h.data.hrDatPtr.Data((length(h.data.Map(:))*(h.data.H_eflag(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(h.data.H_eflag(indE) +eLen - 1)) );
                    dO = reshape(dO, [], eLen);
                    dO = mean(dO(mask(:) == 1, :), 1);
                    dR = reshape(dR, [], eLen);
                    dR = mean(dR(mask(:) == 1, :), 1);
                    
                     %Detrend
                     Pstart = median(dO(1:25));
                     Pend = median(dO((end-24):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                     dO = dO - L;
                     Pstart = median(dR(1:25));
                     Pend = median(dR((end-24):end));
                     m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                     L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
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
                     line(ax, [0 0], [-25 25], 'Color', 'k', 'LineStyle','--');
                     line(ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [-25 25], 'Color', 'k', 'LineStyle','--');
                     xlim([T(1), T(end)]);
                     ylim([mini, maxi]);
                     %Save figure for colours:
                     figName = ['Hb_' h.data.ROIs{indR}.name  '_Evnt_' int2str(indE)];
                     saveas(fig, [h.paths.Graphs figName], 'png');
                     close(fig);
                end
            end
        end
        delete(GraphsDlg);
          %for each event (h.data.R_eflag, h.data.G_eflag, ...
             %graph HbO, HbR and HbT
             %graph Colours
             %graph flow
          %end
          
          %graph mean of events (Colours & Hbs)
          %graph mean of flow
        %end
        
        %Video sequence of each Colour
        %Video sequence of HbO, HbR & HbT
        %Video sequence of flow
        
        % 4 x 4 of each colour, Hbs and flow (7 figs of 4x4)
         
        
    end

    function exportXLS(~,~,~)
         %Waiting Dlg...                
        ExportDlg = dialog('Position',[500 500 250 150],'Name','Export');
        uicontrol('Parent', ExportDlg, 'Style','text',...
            'Position',[20 80 210 40], 'String', 'Exporting data...');
        pause(0.1);
        
        eLen = 5*floor(h.data.Stim.InterStim_min);
        T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, eLen);
        array = zeros(eLen, size(h.data.ROIs,2)*length(h.data.EvntList)+1, 'single');
        filename = [h.paths.FolderName filesep 'DataExport.xls'];
        array(:, 1) = T; names = {'T'};
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
        delete(ExportDlg);
    end

    function OpenFolder(~, ~, ~)
        h.paths.FolderName = uigetdir();
        
        h.flags.bsaveROIS = false; %flag to know if any changes were made to ROIs
        h.flags.bsaveEvnts = false;
        h.flags.Stim = false;
        h.flags.IsThereHbO = false;
        h.flags.IsThereHbR = false;
        h.flags.IsThereHbT = false;
        h.flags.IsThereFlow = false;
        h.flags.IsThereGreen = false;
        h.flags.IsThereRed = false;
        h.flags.IsThereYellow = false;
               
        % Files Path
        h.paths.HbFile = [h.paths.FolderName filesep 'Data_Hbs.mat'];
        h.paths.ROIsFile = [h.paths.FolderName filesep 'ROIs.mat'];
        h.paths.EVNTsFile = [h.paths.FolderName filesep 'Events.mat'];
        h.paths.StimProto = [h.paths.FolderName filesep 'StimParameters.mat'];
        h.paths.Graphs = [h.paths.FolderName filesep 'Graphs' filesep];
        h.paths.Flow = [h.paths.FolderName filesep 'Flow_infos.mat'];
        
        if( exist(h.paths.Graphs , 'dir') )
            ButtonName = questdlg('A folder containing figures already exist. Do you want to overwrite it?', ...
                'Figures folder', ...
                'Yes', 'Change', 'Yes');
            switch ButtonName,
                case 'Yes',
                    disp('Erasing old figures...');
                    delete([h.paths.Graphs '*.*']);
                case 'Change',
                    dname = uigetdir(h.paths.FolderName);
                    h.paths.Graphs = dname;
            end % switch
            
        else
            mkdir(h.paths.Graphs);
        end
        
        %Waiting Dlg...                
        opendlg = dialog('Position',[500 500 250 150],'Name','Loading...');
        uicontrol('Parent', opendlg, 'Style','text',...
            'Position',[20 80 210 40], 'String', 'Loading data. Please wait...');
        pause(0.1);
        
        %Load stimulation parameters
        if(  exist(h.paths.StimProto, 'file') )
            h.data.Stim = load(h.paths.StimProto);
            h.flags.Stim = true;
        else
            h.flags.Stim = false;
        end
        
        if( exist(h.paths.Flow, 'file') )
             h.flags.IsThereFlow = true;
             h.data.fInfo = matfile(h.paths.Flow);
             h.data.fDatPtr = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
             h.data.fDatPtr = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
             Start = find(diff(h.data.fInfo.Stim(1,(5*h.data.Stim.PreStimLength):end),1,2) > 0);
             h.data.F_eflag = Start;
        else
             disp('No flow measures for this experiment!');
             h.flags.IsThereFlow = false;
        end
        
        if( exist(h.paths.HbFile, 'file') )
            h.flags.IsThereHbO = true;
            h.flags.IsThereHbR = true;
            h.flags.IsThereHbT = true;
            
            h.data.HBinfos = matfile(h.paths.HbFile);
            h.data.hoDatPtr = memmapfile(h.data.HBinfos.datFileHbO, 'Format', 'single');
            h.data.hrDatPtr = memmapfile(h.data.HBinfos.datFileHbR, 'Format', 'single');
            Start = find(diff(h.data.HBinfos.Stim(1,(5*h.data.Stim.PreStimLength):end),1,2) > 0);
            h.data.H_eflag = Start;
        else            
            disp('No Hb concentrations were computed for this experiment!');
            h.flags.IsThereHbO = false;
            h.flags.IsThereHbR = false;
            h.flags.IsThereHbT = false;
        end
  
        Ts = 1e6; Map = [];
        RawDatFiles = dir([h.paths.FolderName filesep 'Data_*.mat']);
        if( isempty(RawDatFiles) )
            %No compatible files were found
            h.flags.IsThereHbO = false;
            h.flags.IsThereHbR = false;
            h.flags.IsThereHbT = false;
            h.flags.Stim = false;
            
            disp(['No data files found in ' FolderName ' Folder.']);
            disp('There is nothing to show you... sorry!');
            return;
        end
        %Green channel:
        if( ~isempty(strfind([RawDatFiles.name],'green')) )
            h.flags.IsThereGreen = true;
            Dat_Gptr = matfile([h.paths.FolderName filesep 'Data_green.mat']);
            nrows = Dat_Gptr.datSize(1,1);
            ncols = Dat_Gptr.datSize(1,2);
            nframes = Dat_Gptr.datLength;
            Ws = ncols;
            Hs = nrows;
            Ts = nframes;
            Start = find(diff(Dat_Gptr.Stim(1,(5*h.data.Stim.PreStimLength):end),1,2) > 0);
            h.data.G_eflag = Start;
            h.data.gDatPtr = memmapfile(Dat_Gptr.datFile,...
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
            nrows = Dat_Yptr.datSize(1,1);
            ncols = Dat_Yptr.datSize(1,2);
            nframes = Dat_Yptr.datLength;
            Ws = ncols;
            Hs = nrows;
            Ts = min(Ts, nframes);
            Start = find(diff(Dat_Yptr.Stim(1,(5*h.data.Stim.PreStimLength):end),1,2) > 0);
            h.data.yDatPtr = memmapfile(Dat_Yptr.datFile,...
                'Format', 'single');
            h.data.Y_eflag = Start;
            if( isempty(Map) )
                Map = reshape(h.data.yDatPtr.Data(1:(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        
        %Red channel:
        if( ~isempty(strfind([RawDatFiles.name],'red')) )
            h.flags.IsThereRed = true;
            Dat_Rptr = matfile([h.paths.FolderName filesep 'Data_red.mat']);
            nrows = Dat_Rptr.datSize(1,1);
            ncols = Dat_Rptr.datSize(1,2);
            nframes = Dat_Rptr.datLength;
            Ws = ncols;
            Hs = nrows;
            Ts = min(Ts, nframes);
            Start = find(diff(Dat_Rptr.Stim(1,(5*h.data.Stim.PreStimLength):end),1,2) > 0);
            h.data.R_eflag = Start;
            h.data.rDatPtr = memmapfile(Dat_Rptr.datFile,...
                'Format', 'single');
            if( isempty(Map) )
                Map = reshape(h.data.rDatPtr.Data(1:(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        h.data.NCols = Ws;
        h.data.NRows = Hs;
        h.data.NFrames = Ts;
        
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
        E = ones(1, h.data.Stim.NbStim);
        if(  exist(h.paths.EVNTsFile, 'file')  )
            load([h.paths.FolderName filesep 'Events.mat']); 
        end
        h.data.EvntList = E;
        clear E;
        
        Str = {};
        if( h.flags.IsThereGreen )    Str{end+1} = 'Green';   end
        if( h.flags.IsThereRed )    Str{end+1} = 'Red';          end
        if( h.flags.IsThereYellow )    Str{end+1} = 'Yellow'; end
        if( h.flags.IsThereHbO )    Str{end+1} = 'HbO';        end
        if( h.flags.IsThereHbR )    Str{end+1} = 'HbR';         end
        if( h.flags.IsThereHbT )    Str{end+1} = 'HbT';         end
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

        RefreshLoop('All');
        
        delete(opendlg);
    end

    function StartVideo(~,~,~)
        
        if( ~h.flags.VideoPlaying )
            h.ui.VScreen = figure('Name', 'Video', 'Position', [200 200 500 500], 'Visible', 'off');
            colormap(h.ui.VScreen, 'jet');
            h.ui.ScreenAx = axes('Parent', h.ui.VScreen);
            
            Str = {};
            if( h.flags.IsThereGreen )    Str{end+1} = 'Green';   end
            if( h.flags.IsThereRed )    Str{end+1} = 'Red';          end
            if( h.flags.IsThereYellow )    Str{end+1} = 'Yellow'; end
            if( h.flags.IsThereHbO )    Str{end+1} = 'HbO';        end
            if( h.flags.IsThereHbR )    Str{end+1} = 'HbR';         end
            
            [sel, valid] = listdlg('PromptString', 'Select channel:',...
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
            SelectedSrc = Str{sel};
            if( isempty(strfind(SelectedSrc, 'Hb'))  )
                isHb = 0;
                eval(['StartPts = h.data.' SelectedSrc(1) '_eflag;']);
                eval(['data = h.data.' lower(SelectedSrc(1)) 'DatPtr;']);
                h.data.vidClim = [0.99 1.01];
            else
                isHb = 1;
                StartPts = h.data.H_eflag;
                eval(['data = h.data.h' lower(SelectedSrc(3)) 'DatPtr;']);
                h.data.vidClim = [-10 10];
            end
            
            eLen = 5*floor(h.data.Stim.InterStim_min);
            Accum = zeros(size(h.data.Map,1), size(h.data.Map,2), eLen);
            T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, eLen);
            for indE = 1:h.data.Stim.NbStim
                d = data.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                    (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                
                d = reshape(d, size(Accum));
                
                Pstart = median(d(:, :, 1:25),3);
                Pend = median(d(:,:,(end-24):end),3);
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                L = bsxfun(@minus, bsxfun(@plus, Pend, bsxfun(@times, m, permute(T,[1 3 2]))), ...
                    (m*T(round(end - h.data.Stim.PreStimLength/2))));
                if( isHb == 0 )
                    d = d./L;
                else
                    d = d - L;
                end
                
                Accum = Accum + d;
            end
            Accum = Accum./indE;
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
        
        imagesc(h.ui.ScreenAx, squeeze(h.data.Accumulator(:, :, h.data.vidInd)));
        title([h.data.vidChan ' channel at:' num2str(h.data.vidTimeVect(h.data.vidInd)) 's']);
        set(h.ui.ScreenAx, 'CLim', [0.99, 1.01]);
        axis(h.ui.ScreenAx, 'image');
        h.data.vidInd = h.data.vidInd + 1;
        if( h.data.vidInd > h.data.vidLimit )
            h.data.vidInd = 1;
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
            f = fdesign.lowpass('N,F3dB', 4, 0.4, 5);
            hpass = design(f,'butter');
            f = fdesign.lowpass('N,F3dB', 4, 1/60, 5);
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
        plot(h.ui.EventsDispPan.Ax, T, d, 'k');
        line(h.ui.EventsDispPan.Ax, [0 0], [min(d) max(d)],'Color', 'g', 'LineStyle','--');
        line(h.ui.EventsDispPan.Ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [min(d) max(d)],...
            'Color', 'r', 'LineStyle','--');
        if( mean(d) > 0.5 )
            line(h.ui.EventsDispPan.Ax, [T(1) T(end)], [1 1],...
                'Color', 'k', 'LineStyle',':');
        else
            line(h.ui.EventsDispPan.Ax, [T(1) T(end)], [0 0],...
                'Color', 'k', 'LineStyle',':');
        end
        xlim(h.ui.EventsDispPan.Ax, [T(1), T(end)]);
        
        %Creer le checkbox associe
        set(h.ui.EventsDispPan.Cbox,'String',  int2str(idx) );        
        set(h.ui.EventsDispPan.Cbox,'Value',  h.data.EvntList(idx) );        
    end

    function PopulateEvntsDisplay(~, ~, ~)
         if( isempty(h.data.ROIs) )
             return;
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
            StartPts = 0;
            isHbT = 0; isHb = 0;
            if( isempty(strfind(SelectedSrc, 'Hb')) )
                eval(['StartPts = h.data.' SelectedSrc(1) '_eflag;']);
                eval(['data = h.data.' lower(SelectedSrc(1)) 'DatPtr;']);
            else
                StartPts = h.data.H_eflag;
                isHb = 1;
                if( strfind(SelectedSrc, 'R') )
                    eval(['data = h.data.hrDatPtr;']);
                elseif( strfind(SelectedSrc, 'O') )
                    eval(['data = h.data.hoDatPtr;']);
                else
                    isHbT = 1;
                end
            end
            eLen = 5*floor(h.data.Stim.InterStim_min);
            
            mask = 0;
            if( strcmp(SelectedROI, 'AllPixels') )
                mask = ones(size(h.data.Map));
            else
                idx = arrayfun(@(a) strcmp(h.data.ROIs{a}.name, SelectedROI), 1:size(h.data.ROIs,2));
                mask = h.data.ROIs{idx == 1}.mask;                                                                       
            end
            
            h.ui.EventsDispPan.min = 1 - (h.data.Stim.NbStim)*0.6;
            if( h.ui.EventsDispPan.min > 0 )
                h.ui.EventsDispPan.min = 0;
            end
            
            T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, eLen);
            h.data.EventBuf = zeros(h.data.Stim.NbStim + 1, eLen, 'single');
            h.data.EventBuf(1,:) = T;
            for indE = 1:h.data.Stim.NbStim
                
                if( isHbT == 0 )
                    d = data.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                else
                    d1 = hoDatPtr.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                    d2 = hrDatPtr.Data( (length(h.data.Map(:))*(StartPts(indE) - 1) + 1):...
                        (length(h.data.Map(:))*(StartPts(indE) +eLen - 1)) );
                    d = d1+d2;
                end
                d = reshape(d, [], eLen);
                Glob = mean(d, 1);
                d = mean(d(mask(:)==1,:),1);
                
                if( get(h.ui.GlobalOpt, 'Value') )
                    d = d./Glob;
                end
                if( ~isHb && get(h.ui.FilteringOpt, 'Value') )
                    d = FilterData( d, 'IOI');
                end
                
                %Detrend
                Pstart = median(d(1:25));
                Pend = median(d((end-24):end));
                m = ((Pend - Pstart)/(T(end) - T(1) - h.data.Stim.PreStimLength));
                L = m*T + (Pend - m*T(round(end - h.data.Stim.PreStimLength/2)));
                if( isHb == 0 )
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
        d = zeros(1, 5*floor(h.data.Stim.InterStim_min));
        for indE = 1:length(h.data.EvntList)
            if(h.data.EvntList(indE))
                d = d + h.data.EventBuf(indE + 1,:);
            end
        end
        d = d/sum(h.data.EvntList);
        T = linspace(-h.data.Stim.PreStimLength, h.data.Stim.InterStim_min - h.data.Stim.PreStimLength, length(d));
        plot(h.ui.EventsMeanPan.Ax, T, d, 'k');
        line(h.ui.EventsMeanPan.Ax, [0 0], [min(d) max(d)],'Color', 'g', 'LineStyle','--');
        line(h.ui.EventsMeanPan.Ax, [h.data.Stim.StimLength h.data.Stim.StimLength], [min(d) max(d)],...
            'Color', 'r', 'LineStyle','--');
        sID = get(h.ui.ChannelSelector, 'Value');
        sStr = get(h.ui.ChannelSelector, 'String');
        SelectedSrc = sStr{sID};
        if( isempty(strfind(SelectedSrc, 'Hb')) )
            line(h.ui.EventsMeanPan.Ax, [T(1) T(end)], [1 1],...
                'Color', 'k', 'LineStyle',':');
        else
            line(h.ui.EventsMeanPan.Ax, [T(1) T(end)], [0 0],...
                'Color', 'k', 'LineStyle',':');
        end
        xlim(h.ui.EventsMeanPan.Ax,[T(1), T(end)]);
    end

    function bRet = ValidateEvntSrc(cSrc, rSrc)
        if( cSrc(1) == 'G' )
            bRet = h.flags.IsThereGreen;
        elseif( cSrc(1) == 'R' )
            bRet = h.flags.IsThereRed;
        elseif( cSrc(1) == 'Y' )
            bRet = h.flags.IsThereYellow;
        elseif( cSrc(3) == 'O' )
            bRet = h.flags.IsThereHbO;
        elseif( cSrc(3) == 'R' )
            bRet = h.flags.IsThereHbR;
        elseif( cSrc(3) == 'T' )
            bRet = h.flags.IsThereHbT;
        else
            bRet = false;
        end
        
        
        if( strcmp(rSrc, 'AllPixels') )
            bRet = bRet;
        elseif( any(arrayfun(@(r) strcmp(rSrc, h.data.ROIs{r}.name), 1:size(h.data.ROIs,2))) )
            bRet = bRet;
        else
            bRet = false;
        end
    end

    function RefreshLoop(Option)
        if( strcmp(Option, 'ROIs') )
            RefreshROIsList();
            RefreshMapDisplay();
        elseif( strcmp(Option, 'Events') )
            PopulateEvntsDisplay([],[],[]);
        elseif( strcmp(Option, 'All') )
            RefreshROIsList();
            RefreshMapDisplay();
            PopulateEvntsDisplay([],[],[]);
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
            Type = questdlg('ROI selection method:','Method selection',...
                'Circle', 'Polygon', 'Surround', 'Circle');
            
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
                     
            end
            h.data.ROIs{end+1} = struct('name',answer, 'mask', mask, 'color',  [0 0 1]);           
        end
        
        RefreshLoop('ROIs');
    end

    function RemoveROI(~,~,~)
        Tmp = {};
        if( ~isempty(h.data.ROIs) )
            for indR = 1:size(h.data.ROIs,2)
                Tmp{indR} = h.data.ROIs{indR}.name;
            end
        else
            return;
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
            Str = {'AllPixels'};
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
                Str{end+1} = h.data.ROIs{indR}.name;
            end
            set( h.ui.ROIsSelector, 'String', Str );
        end
    end

    function SaveROIs(~,~,~)
        if( h.flags.saveROIS )
            ROIs = h.data.ROIs;
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
        load(h.paths.ROIsFile);
        h.flags.saveROIS = false;
        load(h.paths.EVNTsFile);
        h.EvntList = E;
        h.flags.saveEvnts = false;
        RefreshLoop('ROIs');
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
