function RestingStateCorrelation()

warning off
close all;
%% Variables

% atlas variables
atlas_raw = load('atlas.mat');
atlas = atlas_raw.map_2d;

atlas_bregma = [148.9327  114.5000];
atlas_lambda = [83.7408  114.5000];
atlas_bulb_center = [225.4388  114.5000];
atlas_points = [atlas_bregma; atlas_lambda;atlas_bulb_center];

atlas_Lm1 = [162.4020   87.1649];
atlas_Rm1 = [162.4020   139];
atlas_Lv1 = [82 59];
atlas_Rv1 = [82 169];
atlas_Ls1 = [178 40];
atlas_Rs1 = [178 183];

%handler initiation
h.paths.FolderName = '';
h.flags.IsThereHbO = false;
h.flags.IsThereHbR = false;
h.flags.IsThereHbT = false;
h.flags.IsThereFlow = false;
h.flags.IsThereGreen = false;
h.flags.IsThereYellow = false;
h.flags.IsThereRed = false;

%Map
h.data.Map = [];

%other variables
dataPaths = '';
nFolders = 0;
w = struct;
h.data.brightness = 32;


%default parameters for blend map
blend = struct;

blend.alphaRange.min = -1;
blend.alphaRange.max = 1;
blend.alphaRange.min_1 = -1;
blend.alphaRange.min_2 = 0.5;
blend.alphaRange.max_1 = -0.5;
blend.alphaRange.max_2 = 1;

blend.figIntensity = 1;

blend.rangeType = 1;

%% Interface generation

ui_fig = figure('Name', 'DataSelector','NumberTitle','off','Position', [300 200 600 600]);

uicontrol('Parent',ui_fig,'Style','text','string','add the number of folders you want and load',...
    'position',[10 550 250 50]);

ui_load_PB = uicontrol('Style','pushbutton','parent', ui_fig, 'String', 'Load',...
    'Position', [50 530 100 50],'Callback',@LoadList);

ui_state_edit = uicontrol('Style','edit','Parent',ui_fig,'Position',[40 10 200 40],'String',...
    'Waiting to load folders','Enable','off');

ui_correlate_PB = uicontrol('Style','pushbutton','Parent', ui_fig, 'String', 'Correlate',...
    'Position', [260 10 100 40], 'Callback', @Process);

loadCounter = uicontrol('Style', 'edit', 'Parent',ui_fig,...
    'Position', [160 530 50 50],'Value',0);

radius_edit = uicontrol('Style', 'edit', 'Parent',ui_fig,...
    'Position', [300 530 50 50],'String',7);
uicontrol('Style', 'text', 'Parent',ui_fig,...
    'Position', [240 520 55 50],'String','ROI radius');
h.ui.ROIsMap = axes('Parent', ui_fig, 'Position',[0.05 0.15 0.700 0.700]);

ui_blend_PB = uicontrol('Style','pushbutton','Parent', ui_fig, 'String', 'Blend setup',...
    'Position', [370 10 100 40], 'Callback', @Blend_param);

set(loadCounter,'String',int2str(loadCounter.Value));

set(ui_correlate_PB,'Enable','off');
set(ui_blend_PB,'Enable','off');

%% functions

    function LoadList(~,~,~)
        
        set(ui_load_PB,'Enable','off');
        set(loadCounter,'Enable','off');
        nFolders = str2num(get(loadCounter,'String'));
        if (nFolders <= 0)
            set(loadCounter,'Enable','on');
            set(ui_load_PB,'Enable','on');
            return
        end
        dataPaths = strings(nFolders,1);
        
        %sfor each path,
        for i = 1:nFolders
            
            dataPaths(i,1) = uigetdir();
            
            %%displaying image of the dataset
            
            channelCheck(dataPaths(i,1));
            
            imagesc(h.ui.ROIsMap, round(h.data.anatomy*h.data.brightness));
            colormap(h.ui.ROIsMap, gray(64));
            axis(h.ui.ROIsMap, 'off');
            
            AjustBrightness();
            
            set(ui_state_edit,'String','Select ROI');
            selection = questdlg('Do you want to select the ROI or use the atlas?',...
                'ROI Selection',...
                'Atlas','Manual','Manual');
            
            if (strcmp(selection, 'Manual'))
                SelectROI(dataPaths(i,1));
                
            else
                SelectAtlasROI();
            end
            SelectBrain(dataPaths(i,1));
            
            w(i).paths=h.paths;
            w(i).flags=h.flags;
            w(i).ui=h.ui;
            w(i).data=h.data;
        end
        
        set(ui_correlate_PB,'Enable','on');
        set(ui_state_edit,'String','Click correlate or blend');
        set(ui_blend_PB,'Enable','on');
    end

    function channelCheck(folderPath)
        file_name = char(strcat(folderPath,filesep,'anato.mat'));
        anat = load(file_name);
        h.data.anatomy = anat.im;
        
        h.paths.FolderName = char(folderPath);
        
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
        h.paths.Graphs = [h.paths.FolderName filesep 'Graphs' filesep];
        h.paths.Flow = [h.paths.FolderName filesep 'Flow_infos.mat'];
        
        h.data.AcqFreq = 0;
        if( exist(h.paths.Flow, 'file') )
            h.flags.IsThereFlow = true;
            
            h.data.fInfo = matfile(h.paths.Flow);
            h.data.AcqFreq = h.data.fInfo.Freq;
            h.data.fDatPtr = memmapfile(h.data.fInfo.datFile, 'Format', 'single');
            
        else
            disp('No flow measures for this experiment!');
            h.flags.IsThereFlow = false;
        end
        
        if( exist(h.paths.HbFile, 'file') )
            h.flags.IsThereHbO = true;
            h.flags.IsThereHbR = true;
            h.flags.IsThereHbT = true;
            
            h.data.HBinfos = matfile(h.paths.HbFile);
            h.data.AcqFreq = h.data.HBinfos.Freq;
            h.data.hoDatPtr = memmapfile(h.data.HBinfos.datFileHbO, 'Format', 'single');
            h.data.hrDatPtr = memmapfile(h.data.HBinfos.datFileHbR, 'Format', 'single');
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
            Ts = nframes;
            h.data.gDatPtr = memmapfile(Dat_Gptr.datFile,...
                'Format', 'single');
            if( isempty(Map) )
                Map = reshape(h.data.gDatPtr.Data(1+(ncols*nrows):2*(ncols*nrows)),nrows,[]);
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
            
            h.data.yDatPtr = memmapfile(Dat_Yptr.datFile,...
                'Format', 'single');
            if( isempty(Map) )
                Map = reshape(h.data.yDatPtr.Data(1+(ncols*nrows):2*(ncols*nrows)),nrows,[]);
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
            h.data.rDatPtr = memmapfile(Dat_Rptr.datFile,...
                'Format', 'single');
            if( isempty(Map) )
                Map = reshape(h.data.rDatPtr.Data(1+(ncols*nrows):2*(ncols*nrows)),nrows,[]);
            end
            clear nrows ncols cframes Start
        end
        h.data.NCols = Ws;
        h.data.NRows = Hs;
        h.data.NFrames = Ts;
        
        h.data.Map =  double(Map)./max(double(Map(:)));
        
    end

    function SelectAtlasROI(~,~,~)
        h.data.ROIs = {};
        
        %selecting different parameters
        set(ui_state_edit,'String','Select Lambda');
        point_1 = impoint(h.ui.ROIsMap);
        setColor(point_1,'g');
        h.ui.Lambda = getPosition(point_1);
        set(ui_state_edit,'String','Select Bregma');
        point_2 = impoint(h.ui.ROIsMap);
        setColor(point_2,'y');
        h.ui.Bregma = getPosition(point_2);
        set(ui_state_edit,'String','Select the bulb');
        point_3 = impoint(h.ui.ROIsMap);
        setColor(point_3,'b');
        h.ui.Bulb = getPosition(point_3);
        h.data.anat_points =[h.ui.Bregma;h.ui.Lambda;h.ui.Bulb];
        
        map_1 = h.data.anatomy;
        im_atlas = atlas;
        map_1_points = h.data.anat_points;
        
        size_map_1 = size(map_1);
        size_im_atlas = size(im_atlas);
        im_atlas = imresize(im_atlas,[size_map_1(1),size_map_1(2)]);
        factor_x = size_map_1(2)/size_im_atlas(2);
        factor_y = size_map_1(1)/size_im_atlas(1);
        atlas_point_x = atlas_points(:,1)*factor_x;
        atlas_point_y = atlas_points(:,2)*factor_y;
        mod_atlas_points = [atlas_point_x ,atlas_point_y];
        
        list = {'Left M1','Right M1','Left V1','Right V1','Left S1','Right S1'};
        indx = listdlg('PromptString','Select a region:',...
            'SelectionMode','single',...
            'ListString',list);
        switch indx
            case 1
                atlas_selec = atlas_Lm1;
            case 2
                atlas_selec = atlas_Rm1;
            case 3
                atlas_selec = atlas_Lv1;
            case 4
                atlas_selec = atlas_Rv1;
            case 5
                atlas_selec = atlas_Ls1;
            case 6
                atlas_selec = atlas_Rs1;
        end
        
        atlas_selec_x = atlas_selec(1)*factor_x;
        atlas_selec_y = atlas_selec(2)*factor_y;
        mod_atlas_selec = [atlas_selec_x,atlas_selec_y];
        
        tform = fitgeotrans(mod_atlas_points,map_1_points,'similarity');
        ra=imref2d(size(map_1));
        [mod_atlas,rb] = imwarp(im_atlas,tform);
        
        [t_atlas_m1_x,t_atlas_m1_y]=transformPointsForward(tform,mod_atlas_selec(1),mod_atlas_selec(2));
        [mod_atlas_sel_x,mod_atlas_sel_y]=worldToIntrinsic(rb,t_atlas_m1_x,t_atlas_m1_y);
        
        answer = inputdlg('enter the ROI radius','Radius Selection');
        answer = str2num(char(answer));
        
        if (answer ~= 0 && isinteger(answer))
            radius = answer;
        else
            radius = 5;
        end
        th = 0:pi/50:2*pi;
        roi_x = radius * cos(th) + mod_atlas_sel_x;
        roi_y = radius * sin(th) + mod_atlas_sel_y;
        
        figure('Name','Atlas and picture overlay','NumberTitle','off');
        imshowpair(map_1,ra,mod_atlas,rb);
        hold on;
        plot(h.ui.Lambda(1),h.ui.Lambda(2),'bo');
        plot(h.ui.Bregma(1),h.ui.Bregma(2),'ro');
        plot(h.ui.Bulb(1),h.ui.Bulb(2),'yo')
        colormap gray;
        hold off;
        
        imagesc(h.ui.ROIsMap, mod_atlas);
        colormap(h.ui.ROIsMap, gray(64));
        axis(h.ui.ROIsMap, 'off');
        
        mask = poly2mask(roi_x,roi_y,size_map_1(1),size_map_1(2));
        
        h.data.ROI.mask = mask;
    end

    function SelectROI(dataPath,~,~)
        h.data.ROIs = {};
        
        imagesc(h.ui.ROIsMap, round(h.data.anatomy*h.data.brightness));
        colormap(h.ui.ROIsMap, gray(64));
        axis(h.ui.ROIsMap, 'off');
        
        addMore = true;
        nbROI = 0;
        
        while (addMore)
            nbROI = nbROI +1;
            Choice = questdlg('Load or create a new ROI','ROI selection',...
                'Load','Create','Create');
            if(strcmp(Choice,'Load'))
                [filename,path] = uigetfile(h.paths.FolderName);
                file_path = char(strcat(path,filesep,filename));
                data = load(file_path);
                mask = data.ROIs{1,1}.mask;
                name = data.ROIs{1,1}.name;
            end
            
            if(strcmp(Choice,'Create'))
                Type = questdlg('ROI selection method:','Method selection',...
                    'Circle', 'Surround','Point' ,'Circle');
                
                h_im = get(h.ui.ROIsMap,'Children');
                
                switch Type
                    case 'Circle'
                        e = imellipse(h.ui.ROIsMap);
                        mask = createMask(e, h_im(end));
                        delete(e);
                        clear e;
                    case 'Point'
                        p = impoint(h.ui.ROIsMap);
                        point = getPosition(p);
                        radius = str2num(get(radius_edit,'String'));
                        th = 0:pi/50:2*pi;
                        roi_x = radius * cos(th) + point(1) ;
                        roi_y = radius * sin(th) + point(2);
                        mask = poly2mask(roi_x,roi_y, size(h.data.anatomy,1),size(h.data.anatomy,2));
                        delete(p);
                        clear p;
                    case 'Surround'
                        [orig, valid] = listdlg('PromptString', 'From which ROI?',...
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
                
                selection = questdlg('Do you want to save the ROI?',...
                    'Before continuing...',...
                    'Yes','No','No');
                if (strcmp(selection, 'Yes'))
                    list ={'left_frontal','left_motor','left_cingulate','left_somato','left_retrospin','left_visual'...
                        ,'right_frontal','right_motor','right_cingulate','right_somato','right_retrospin','right_visual'};
                    indx = listdlg('PromptString','Select a ROI name',...
                        'SelectionMode','single',...
                        'ListString',list);
                    switch(indx)
                        case 1
                            name = list(1);
                        case 2
                            name = list(2);
                        case 3
                            name = list(3);
                        case 4
                            name = list(4);
                        case 5
                            name = list(5);
                        case 6
                            name = list(6);
                        case 7
                            name = list(7);
                        case 8
                            name = list(8);
                        case 9
                            name = list(9);
                        case 10
                            name = list(10);
                        case 11
                            name = list(11);
                        case 12
                            name = list(12);
                    end
                    file_name = char(strcat(dataPath,filesep,name,'.mat'));
                    ROIp.color = [1 0 0];
                    ROIp.name = name;
                    ROIp.mask = mask;
                    ROIs = {ROIp};
                    save(file_name,'ROIs');
                end
            end
            h.data.ROI.mask(:,:,nbROI) = mask;
            h.data.ROI.name(:,:,nbROI) = name;
            addROISelection = questdlg('Do you want to add another ROI?',...
                'Before continuing...',...
                'Yes','No','No');
            if(strcmp(addROISelection,'No'))
                addMore = false;
                h.data.ROI.nbROI = nbROI;
            end
        end
    end

    function SelectBrain(dataPath,~,~)
        imagesc(h.ui.ROIsMap,round(h.data.anatomy*h.data.brightness));
        colormap(h.ui.ROIsMap, gray(64));
        axis(h.ui.ROIsMap, 'off');
        
        Choice = questdlg('Load or create a new brain selection','brain selection',...
            'Load','Create','Create');
        if(strcmp(Choice,'Load'))
            [filename,path] = uigetfile(h.paths.FolderName);
            file_path = char(strcat(path,filesep,filename));
            data = load(file_path);
            mask = data.mask;
        end
        
        if(strcmp(Choice,'Create'))
            h_im = get(h.ui.ROIsMap,'Children');
            p = impoly(h.ui.ROIsMap);
            mask = createMask(p, h_im(end));
            delete(p);
            clear p;
            
            selection = questdlg('Do you want to save the brain selection?',...
                'Before continuing...',...
                'Yes','No','No');
            if (strcmp(selection, 'Yes'))
                file_name = char(strcat(dataPath,filesep,'brain.mat'));
                save(file_name,'mask');
            end
        end
        h.data.ROI.brainmask = mask;
    end

    function Correlation(channel,dataPath)
        
        if( isempty(strfind(channel, 'Hb')) && isempty(strfind(channel, 'Flow')) )
            eval(['data = h.data.' lower(channel(1)) 'DatPtr;']);
        elseif( ~isempty(strfind(channel, 'Flow')) )
            data = h.data.fDatPtr;
        else
            if( channel == 'HbR' )
                eval('data = h.data.hrDatPtr;');
            elseif( channel == 'HbO' )
                eval('data = h.data.hoDatPtr;');
            elseif( channel == 'HbT')
                eval('data = h.data.hoDatPtr;');
                data1 = data;
                eval('data = h.data.hrDatPtr;');
                data2 = data;
            end
        end
        
        nbROI = h.data.ROI.nbROI;
        
        for t = 1:nbROI
            mask = h.data.ROI.mask(:,:,t);
            name = h.data.ROI.name(:,:,t);
            brainmask = h.data.ROI.brainmask;
            if(strcmp(channel ,'HbT'))
                d1 = reshape(data1.Data, h.data.NRows, h.data.NCols, []);
                d2 = reshape(data2.Data, h.data.NRows, h.data.NCols, []);
                d = d1 + d2;
            else
                d = reshape(data.Data, h.data.NRows, h.data.NCols, []);
            end
            % definir filtre butter
            
            fc=h.data.AcqFreq;
            fs1=0.009;
            fs2=0.08;
            [b,a] = butter(4,[fs1/(fc/2) fs2/(fc/2)],'bandpass');
            
            % étape 0 : filtre
            if(~strcmp(channel,'HbO')&& ~strcmp(channel,'HbR') && ~strcmp(channel,'HbT'))
                %                 for i = 1:h.data.NFrames-1
                %                     d(:,:,i) = homomorphic_filter(d(:,:,i),10,5,0.5,1.5);
                %                 end
                %             else
                for i = 1:h.data.NFrames-1
                    d(:,:,i) = imgaussfilt(d(:,:,i),3);
                end
            end
            % etape 1 : Definir le signal moyen du cerveau
            
            for i=1:h.data.NFrames-1
                Mask_cerveau = d(:,:,i); % h.data.ROIs de tout le cerveau
                Vmoy(i,1)= mean(Mask_cerveau(brainmask(:) == 1));
            end
            
            Vmoy_filt= filtfilt(b,a,double(Vmoy)); %filtrer Vmoy
            
            % etape 2 : definir signal moyen pour region d'interet
            
            for i=1:h.data.NFrames-1
                Masked_reg = d(:,:,i);
                Vref(i,1)= mean(Masked_reg(logical(mask(:)==1)));
            end
            Vref_filt= filtfilt(b,a,double(Vref)); % filtrer Vref
            Beta= Vref_filt\Vmoy_filt;
            Vref_chapo= Vref_filt- (Beta.*Vmoy_filt);
            
            % etape 3 : construire carte de correlation
            for j= 1:h.data.NCols
                for k=1:h.data.NRows
                    d_filt= filtfilt(b,a,double(squeeze(d(k,j,1:size(Vmoy,1))))); % filter voxel
                    beta= d_filt\Vmoy_filt;
                    Vcompare = d_filt-(beta.*Vmoy_filt);
                    StatCart(j,k) = (corr2(Vcompare,Vref_chapo));
                end
            end
            %Step 4: building correlation map with ROI
            roiboundaries = bwboundaries(mask);
            xy = roiboundaries{1};
            x = xy(:,2);
            y = xy(:,1);
            
            fig_1 = figure('visible','off');
            corr_map = brainmask.*StatCart';
            im = imagesc(corr_map,[-1,1]);
            hold on
            alphadat = brainmask;
            set(im, 'AlphaData',alphadat);
            plot(x,y,'w-','LineWidth',1);
            colormap jet
            colorbar
            axis 'off';
            
            h.data.StatCart = double(StatCart);
            path = char(strcat(dataPath,filesep,'image_correlation',filesep,name));
            if(~exist(path, 'dir'))
                mkdir(path);
            end
            
            file_name = char(strcat(path,filesep,channel));
            print(fig_1,file_name,'-djpeg');
            file_name = char(strcat(path,filesep,channel,'.mat'));
            save(file_name,'corr_map');
            delete(fig_1);
            
            %Step 6: building overlay blend with correlation map
            fig_3 = figure('visible','off');
            if( blend.rangeType == 1)
                pat_overlay_blend(round(h.data.anatomy*h.data.brightness), corr_map,brainmask,[-1 1],...
                    [blend.alphaRange.min, blend.alphaRange.max], jet(256), blend.figIntensity);
            end
            if( blend.rangeType == 2)
                pat_overlay_blend(round(h.data.anatomy*h.data.brightness), corr_map,brainmask,[-1 1],...
                    [blend.alphaRange.min_1, blend.alphaRange.max_1,blend.alphaRange.min_2, blend.alphaRange.max_2],...
                    jet(256), blend.figIntensity);
            end
            hold on;
            plot(x,y,'g--','LineWidth',1);
            colormap jet
            c=colorbar;
            set(c,'color','w');
            axis 'off';
            file_name = char(strcat(path,filesep,channel,'_blend'));
            print(fig_3,file_name,'-djpeg');
            file_name = char(strcat(path,filesep,channel,'_blend.mat'));
            save(file_name,'fig_3');
            delete(fig_3);
        end
        % step 7: Build Map of regional node degree
%         node_degree_map = zeros(h.data.NRows,h.data.NCols);
%         squee_Vref_chapo = cell(h.data.NCols*h.data.NRows,1);
%         idx_squee = 1;
%         for i = 1:h.data.NRows
%             for j =1:h.data.NCols
%                 d_filt= filtfilt(b,a,double(squeeze(d(k,j,1:size(Vmoy,1))))); % filter voxel
%                 beta= d_filt\Vmoy_filt;
%                 squee_Vref_chapo{idx_squee} = d_filt-(beta.*Vmoy_filt);
%                 idx_squee = idx_squee +1;
%             end
%         end
%         for idxVref = 1:h.data.NCols*h.data.NRows
%             if(brainmask(idxVref))
%                 test = squee_Vref_chapo{idxVref};
%                 compteur = 0;
%                 for idxVcomp= 1:h.data.NRows*h.data.NCols
%                     if(brainmask(idxVcomp))
%                         temp_corr = corr2(squee_Vref_chapo{idxVcomp},test);
%                         z_corr = 0.5*log((1+temp_corr)./(1-temp_corr));
%                         if(z_corr >= 0.4)
%                             compteur = compteur +1;
%                         end
%                     end
%                 end
%                 node_degree_map(idxVref)= compteur-1; % removing his own pixel correlation
%             end
%         end
%         figure();
%         imagesc(node_degree_map);
%         colormap jet;
    end

    function Process(~,~,~)
        delete(ui_fig);
        k = waitbar(0,'Please wait...');
        for i=1:nFolders
            folder_name = char(strcat(dataPaths(i,1),filesep,'image_correlation'));
            mkdir(folder_name);
            h = w(i);
            
            if (h.flags.IsThereHbO == true)
                Correlation('HbO',dataPaths(i,1));
                waitbar((7*(i-1)+1)/(7*nFolders),k);
            end
            if (h.flags.IsThereHbR == true)
                Correlation('HbR',dataPaths(i,1));
                waitbar((7*(i-1)+2)/(7*nFolders),k);
            end
            
            if (h.flags.IsThereHbO == true && h.flags.IsThereHbR == true)
                Correlation('HbT',dataPaths(i,1));
                waitbar((7*(i-1)+3)/(7*nFolders),k);
            end
            
            if (h.flags.IsThereFlow == true)
                Correlation('Flow',dataPaths(i,1));
                waitbar((7*(i-1)+4)/(7*nFolders),k);
            end
            if (h.flags.IsThereYellow == true)
                Correlation('Yellow',dataPaths(i,1));
                waitbar((7*(i-1)+5)/(7*nFolders),k);
            end
            if (h.flags.IsThereRed == true)
                Correlation('Red',dataPaths(i,1));
                waitbar((7*(i-1)+6)/(7*nFolders),k);
            end
            if (h.flags.IsThereGreen == true)
                Correlation('Green',dataPaths(i,1));
                waitbar((7*(i-1)+7)/(7*nFolders),k);
            end
        end
        delete(k);
    end

    function Blend_param(~,~,~)
        ui_blend = figure('Name', 'Blend parameters selector','NumberTitle','off','Position', [300 200 350 350]);
        
        ui_blend_slide_max = uicontrol('Style','slider','Min',-1,'Max',1,'SliderStep',[0.01 0.10],'Parent',ui_blend,...
            'Position',[10 220 100 20],'Value',blend.alphaRange.max,'Callback',@Display_max,'BackgroundColor','w');
        ui_blend_edit_max = uicontrol('Style','edit','Enable','off','Position',[120 220 35 20],'String',...
            get(ui_blend_slide_max,'Value'),'BackgroundColor','w');
        ui_blend_slide_min = uicontrol('Style','slider','Min',-1,'Max',1,'SliderStep',[0.01 0.10],'Parent',ui_blend,...
            'Position',[10 270 100 20],'Value',blend.alphaRange.min,'Callback',@Display_min,'BackgroundColor','w');
        ui_blend_edit_min = uicontrol('Style','edit','Enable','off','Position',[120 270 35 20],'String',...
            get(ui_blend_slide_min,'Value'),'BackgroundColor','w');
        ui_blend_slide_max_1 = uicontrol('Style','slider','Min',-1,'Max',1,'SliderStep',[0.01 0.10],'Parent',ui_blend,...
            'Position',[170 220 100 20],'Value',blend.alphaRange.max_1,'Callback',@Display_max_1,'BackgroundColor','w');
        ui_blend_edit_max_1 = uicontrol('Style','edit','Enable','off','Position',[280 220 35 20],'String',...
            get(ui_blend_slide_max_1,'Value'),'BackgroundColor','w');
        ui_blend_slide_min_1 = uicontrol('Style','slider','Min',-1,'Max',1,'SliderStep',[0.01 0.10],'Parent',ui_blend,...
            'Position',[170 270 100 20],'Value',blend.alphaRange.min_1,'Callback',@Display_min_1,'BackgroundColor','w');
        ui_blend_edit_min_1 = uicontrol('Style','edit','Enable','off','Position',[280 270 35 20],'String',...
            get(ui_blend_slide_min_1,'Value'),'BackgroundColor','w');
        ui_blend_slide_max_2 = uicontrol('Style','slider','Min',-1,'Max',1,'SliderStep',[0.01 0.10],'Parent',ui_blend,...
            'Position',[170 170 100 20],'Value',blend.alphaRange.max_2,'Callback',@Display_max_2,'BackgroundColor','w');
        ui_blend_edit_max_2 = uicontrol('Style','edit','Enable','off','Position',[280 170 35 20],'String',...
            get(ui_blend_slide_max_2,'Value'),'BackgroundColor','w');
        ui_blend_slide_min_2 = uicontrol('Style','slider','Min',-1,'Max',1,'SliderStep',[0.01 0.10],'Parent',ui_blend,...
            'Position',[170 120 100 20],'Value',blend.alphaRange.min_2,'Callback',@Display_min_2,'BackgroundColor','w');
        ui_blend_edit_min_2 = uicontrol('Style','edit','Enable','off','Position',[280 120 35 20],'String',...
            get(ui_blend_slide_min_2,'Value'),'BackgroundColor','w');
        
        ui_blend_text_min = uicontrol('Style','text','String','minimum of range','Position',...
            [10 285 100 20]);
        ui_blend_text_max = uicontrol('Style','text','String','maximum of range','Position',...
            [10 235 100 20]);
        
        ui_blend_text_1 = uicontrol('Style','text','String','1 Range','Position',...
            [10 320 100 20]);
        ui_blend_text_2 = uicontrol('Style','text','String','2 Ranges','Position',...
            [170 320 100 20]);
        
        ui_blend_text_min_1 = uicontrol('Style','text','String','minimum of range 1','Position',...
            [170 285 100 20]);
        ui_blend_text_max_1 = uicontrol('Style','text','String','maximum of range 1','Position',...
            [170 235 100 20]);
        ui_blend_text_min_2 = uicontrol('Style','text','String','minimum of range 2','Position',...
            [170 135 100 20]);
        ui_blend_text_max_2 = uicontrol('Style','text','String','maximum of range 2','Position',...
            [170 185 100 20]);
        
        ui_blend_text_op = uicontrol('Style','text','String','Other parameters','Position',...
            [10 185 100 20]);
        
        ui_blend_text_choice = uicontrol('Style','text','String','number of ranges','Position',...
            [10 150 100 20]);
        ui_blend_range_choice = uicontrol('Style','popup','String',{'1','2'},'Position',...
            [110 155 30 20]);
        
        ui_blend_figint_text = uicontrol('Style','text','String','Weight of the base layer','Position',...
            [10 125 120 20]);
        ui_blend_slide_int = uicontrol('Style','slider','Min',0,'Max',1,'SliderStep',[0.01 0.10],'Parent',ui_blend,...
            'Position',[10 110 100 20],'Value',blend.figIntensity,'Callback',@Display_intensity,'BackgroundColor','w');
        ui_blend_edit_int = uicontrol('Style','edit','Enable','off','Position',[120 110 35 20],'String',...
            get(ui_blend_slide_int,'Value'),'BackgroundColor','w');
        ui_blend_done_PB = uicontrol('Style','pushbutton','Parent',ui_blend,'String','Done','Position',...
            [125 30 100 30],'Callback',@Done);
        
        function Display_max(~,~,~)
            value = get(ui_blend_slide_max,'Value');
            set(ui_blend_edit_max,'String',value);
        end
        
        function Display_min(~,~,~)
            value = get(ui_blend_slide_min,'Value');
            set(ui_blend_edit_min,'String',value);
        end
        
        function Display_max_1(~,~,~)
            value = get(ui_blend_slide_max_1,'Value');
            set(ui_blend_edit_max_1,'String',value);
        end
        
        function Display_min_1(~,~,~)
            value = get(ui_blend_slide_min_1,'Value');
            set(ui_blend_edit_min_1,'String',value);
        end
        
        function Display_max_2(~,~,~)
            value = get(ui_blend_slide_max_2,'Value');
            set(ui_blend_edit_max_2,'String',value);
        end
        
        function Display_min_2(~,~,~)
            value = get(ui_blend_slide_min_2,'Value');
            set(ui_blend_edit_min_2,'String',value);
        end
        
        function Display_intensity(~,~,~)
            value = get(ui_blend_slide_int,'Value');
            set(ui_blend_edit_int,'String',value);
        end
        
        function Done(~,~,~)
            blend.alphaRange.min = get(ui_blend_slide_min,'Value');
            blend.alphaRange.max = get(ui_blend_slide_max,'Value');
            blend.alphaRange.min_1 = get(ui_blend_slide_min_1,'Value');
            blend.alphaRange.min_2 = get(ui_blend_slide_min_2,'Value');
            blend.alphaRange.max_1 = get(ui_blend_slide_max_1,'Value');
            blend.alphaRange.max_2 = get(ui_blend_slide_max_2,'Value');
            
            blend.figIntensity = get(ui_blend_slide_int,'Value');
            
            blend.rangeType = get(ui_blend_range_choice,'Value');
            
            delete(ui_blend);
        end
    end

    function AjustBrightness()
        prompt = {'Enter image brightness value:'};
        title = 'Input';
        dims = [1 35];
        definput = {int2str(h.data.brightness)};
        answer = inputdlg(prompt,title,dims,definput);
        h.data.brightness = str2num(cell2mat(answer));
        imagesc(h.ui.ROIsMap, round(h.data.anatomy*h.data.brightness));
        colormap(h.ui.ROIsMap, gray(64));
        axis(h.ui.ROIsMap, 'off');
    end

end