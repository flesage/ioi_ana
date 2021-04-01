function Ret = VisuData()
dParams.Folder = '';
dParams.Chan = '';
OldChan = '';
Data = [];
cData = [];
Infos = [];
currentPixel = [1 1];

% Figures:
%Principale
hParams.figP = uifigure('Name', 'Parametres', 'NumberTitle','off',...
    'Position', [20 300 250 700], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'CloseRequestFcn', @FermeTout);
%Raw data:
hParams.figR = uifigure('Name', 'Images', 'NumberTitle','off',...
    'Position', [285 510 500 500], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'Visible','off',...
    'WindowButtonDownFcn', @ChangePtPos, 'CloseRequestFcn', @NeFermePas);
%Corr data:
hParams.figC = uifigure('Name', 'Correlation', 'NumberTitle','off',...
    'Position', [285 200 250 250], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'Visible','off', ...
    'CloseRequestFcn', @NeFermePas);
%Decours Temporel:
hParams.figT = uifigure('Name', 'Signal Temporel', 'NumberTitle','off',...
    'Position', [600 200 750 250], 'Color', 'w', 'MenuBar', 'none',...
    'ToolBar', 'none', 'Resize', 'off', 'Visible','off',...
    'CloseRequestFcn', @NeFermePas);

% GUI
%Pour le chemin d'acces vers le data a visualiser:
hParams.ExpLabel = uilabel(hParams.figP, 'Text','Experience:',...
    'Position',[5, 655, 100, 35], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left');
hParams.ExpEdit = uieditfield(hParams.figP, 'text',...
    'Value', 'Choisir un dossier', 'Position',[5, 640, 175, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left');
hParams.ExpPb = uibutton(hParams.figP, 'push',...
    'Text', '', 'Position',[200, 640, 35, 29], 'BackgroundColor', 'w',...
    'Icon', 'FolderIcon.png', 'ButtonPushedFcn', @ChangeFolder);

%Pre-Analyse:
hParams.PreAPB = uibutton(hParams.figP, 'push', ...
    'Text','Pre-Analyse', 'Position',[45, 590, 150, 35],...
    'BackgroundColor','w', 'FontName', 'Calibri', 'FontSize', 12,...
    'ButtonPushedFcn', @RunPreAna, 'visible', 'off');
hParams.PreALabel = uilabel(hParams.figP, 'Text','Pre-Analyse en cours. Patientez svp...',...
    'Position',[20, 550, 200, 35], 'BackgroundColor','w', ...
    'FontName', 'Calibri', 'FontSize', 12, 'visible', 'off');

%Cannal a utiliser:
hParams.ChanLabel = uilabel(hParams.figP, 'Text','Canal d''imagerie:',...
    'Position',[5, 590, 200, 35], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.ChanPopMenu = uidropdown(hParams.figP, 'Items', {'Choisir'},...
    'Position',[5, 575, 175, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10, 'visible', 'off',...
    'ValueChangedFcn', @OuvrirData);

%Type d'experience:
hParams.TypeLabel = uilabel(hParams.figP, 'Text','Type d''enregistrement',...
    'Position',[5, 525, 200, 35], 'BackgroundColor','w', 'FontName', 'Calibri',...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'visible', 'off');
hParams.TypePopMenu = uidropdown(hParams.figP, 'Items',{'Choisir', 'RestingState', 'Episodique'},...
    'Position',[5, 505, 175, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,...
    'visible', 'off', 'ValueChangedFcn', @ChangeType);

%Interaction Communes:
hParams.dFsFPB = uibutton(hParams.figP, 'push',...
    'Text','DF/F', 'Position',[45, 450, 150, 35], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12,...
    'visible', 'off', 'ButtonPushedFcn', @DFsF);
hParams.GSRPB = uibutton(hParams.figP, 'push',...
    'Text','GSR', 'Position',[45, 400, 150, 35], 'BackgroundColor','w',...
    'FontName', 'Calibri', 'FontSize', 12,...
    'visible', 'off', 'ButtonPushedFcn', @GSR);

% Interaction RS:
% hParams.CorrPB = uibutton(hParams.figP, 'push',...
%     'Text','Carte Corr', 'Position',[45, 350, 150, 35], 'BackgroundColor','w',...
%     'FontName', 'Calibri', 'FontSize', 12,...
%     'visible', 'off', 'ButtonPushedFcn', @CorrMap);

% Visualisation des images brutes:
% Graph:
hParams.axR1 = uiaxes(hParams.figR, 'Position', [67.5 100 375 375]);
% Boutons:
hParams.CurrentImageSl = uislider( hParams.figR,... 
    'Value', 1, 'Limits', [1 2], 'MajorTicks', [1 2],...
    'Position',[50, 90, 400, 3], 'ValueChangedFcn', @ChangeImage);
hParams.CI_MinLabel = uilabel( hParams.figR, 'Text', 'Minimum:', ...
    'Position',[100, 15, 75, 25],'FontName', 'Calibri', 'FontSize', 12,...
    'HorizontalAlignment', 'left');
hParams.CI_MaxLabel = uilabel( hParams.figR, 'Text', 'Maximum:', ...
    'Position',[300, 15, 75, 25],'FontName', 'Calibri', 'FontSize', 12,...
    'HorizontalAlignment', 'left');
hParams.CI_Min_Edit = uieditfield( hParams.figR, 'numeric',...
    'Value', 1, 'Position',[175, 15, 50, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left',...
    'ValueChangedFcn', @AdjustImage);
hParams.CI_Max_Edit = uieditfield( hParams.figR, 'numeric',...
    'Value', 4096, 'Position',[375, 15, 50, 25], 'BackgroundColor', 'w',...
    'FontName', 'Calibri', 'FontSize', 10,  'HorizontalAlignment', 'left',...
    'ValueChangedFcn', @AdjustImage);

% Visualisation Corrélation
%Graph:
hParams.axC1 = uiaxes(hParams.figC, 'Position', [5 5 240 240]);

% Visualisation Décours Temporel:
%Graph:
hParams.axT1 = uiaxes(hParams.figT, 'Position', [5 5 740 240]);

% Initialisation de l'interface:
ChangeMode('Ouverture');

% Fonctions et Callbacks:
    function ChangeFolder(~,~,~)
        
        selpath = uigetdir(path);
        if( selpath == 0 )
            return;
        end
        
        if( ~strcmp(selpath(end), filesep) )
            selpath = strcat(selpath, filesep);
        end
        dParams.sFolder = selpath;
        hParams.ExpEdit.Value = selpath;
        ChargerDossier();
    end

    function ChargerDossier()
         %Validation du dossier
        list = dir([dParams.sFolder '*.dat']);
        
        Channels{1} = 'Choisir un Canal';
        for ind = 1:size(list,1) 
            Channels{ind+1} = list(ind).name;
        end
        
        hParams.ChanPopMenu.Items = Channels;
        hParams.ChanPopMenu.Value = Channels{1};
        if( size(list,1) == 0 )
            ChangeMode('PreAna');
        else
            ChangeMode('SelectParams');
        end
        
        CheckDefaultParams();         
    end

    function CheckDefaultParams()
       %S'il y a un fichier mat a loader:
       dParams = matfile([dParams.sFolder  'dParams.mat'],'Writable', true);
       if( ~exist([dParams.Properties.Source 'dParams.mat'],'file') )
          %Sinon, on le cree:
          dParams.sFolder = hParams.ExpEdit.Value;
          dParams.sExpType = 'NA';
          dParams.sStimAI = 'NA';
          dParams.sFichierStim = 'NA';
       else
          %Faut changer les parametres par defauts...
          
       end
    end

    function OuvrirData(~,~,~)
        if( iscell(hParams.ChanPopMenu.Items) )
            dParams.Chan = hParams.ChanPopMenu.Value;
            if( contains(hParams.ChanPopMenu.Items{1}, 'Choisir') )
                hParams.ChanPopMenu.Items = hParams.ChanPopMenu.Items(2:end);
            end
        else
            dParams.Chan = hParams.ChanPopMenu.Value;
        end
        
        if( ~strcmp(dParams.Chan, OldChan) )
            fid = fopen([dParams.sFolder dParams.Chan]);
            Data = fread(fid,inf, 'single=>single');
            Tmp = dir([dParams.sFolder 'Data_*.mat']);
            Infos = matfile(Tmp(1).name);
            Data = reshape(Data, Infos.datSize(1,1), Infos.datSize(1,2),[]);
            fclose(fid);
            
            hParams.GSRPB.Enable = 'on';
            hParams.dFsFPB.Enable = 'on';
        end    
        OldChan = dParams.Chan;
        hParams.dFsFPD.Enable = 'on'; 
        
        if( any(contains(hParams.TypePopMenu.Items, 'Choisir')) )
            ChangeMode('SelectParams')
        elseif( strcmp(dParams.sExpType, 'RestingState') )
            ChangeMode('RestingState');
        else
            ChangeMode('Episodique');
        end
        hParams.figR.Visible = 'on';
       
        hParams.CurrentImageSl.Limits = [1 size(Data,3)];
        hParams.CurrentImageSl.Value = 1;
        hParams.CurrentImageSl.MajorTicks = 1:1000:size(Data,3);
        ChangeImage();
    end

    function ChangeImage(~,~,~)
        Id = round(hParams.CurrentImageSl.Value);
        Im = imresize(squeeze(Data(:,:,Id)),[256 256]);
        imagesc(hParams.axR1, Im);
        caxis(hParams.axR1, [hParams.CI_Min_Edit.Value, hParams.CI_Max_Edit.Value]);
        title(hParams.axR1,['Image #: ' int2str(Id)]);
        axis(hParams.axR1, 'off', 'image');
        hold(hParams.axR1, 'on');
        plot(hParams.axR1, currentPixel(1), currentPixel(2), 'or');
        hold(hParams.axR1, 'off');
        DecoursTemp();    
    end

    function AdjustImage(~,~,~)
        caxis(hParams.axR1, [hParams.CI_Min_Edit.Value, hParams.CI_Max_Edit.Value]);
        DecoursTemp();
    end

    function RunPreAna(~,~,~)
        prompt = {'Binning Spatial (1 si aucun binning):',...
            'Binning Temporel(1 si aucun binning):',...
            'Redifinir une Region d''interet? (0:non; 1:oui)',...
            'Ignorer le signal de stimulation interne du systeme? (0:non; 1:oui)'};
        dlgtitle = 'Pre-Analyse';
        dims = [1 50];
        definput = {'1','1', '0', '0'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        
        if( isempty(answer) )
            return;
        end
        hParams.PreALabel.Visible = 'on';
        pause(0.01);
        try
        ImagesClassification(dParams.sFolder, str2double(answer{1}),...
            str2double(answer{2}), str2double(answer{3}), str2double(answer{4}));
        catch e
            disp(e);
            hParams.PreALabel.Value = 'Une erreur est survenue durant la pre-analyse.';
        end
        hParams.PreALabel.Visible = 'off';
        ChargerDossier();
    end

    function ChangeType(~,~,~)
        dParams.sExpType = hParams.TypePopMenu.Value;
        if( contains(hParams.TypePopMenu.Items{1}, 'Choisir') )
           hParams.TypePopMenu.Items = hParams.TypePopMenu.Items(2:end);
        end
        if( any(contains(hParams.ChanPopMenu.Items, 'Choisir')) )
            ChangeMode('SelectParams')
        elseif( strcmp(dParams.sExpType, 'RestingState') )
            ChangeMode('RestingState');
            CorrMap();
        else
            ChangeMode('Episodique');
        end
    end
        
    function ChangeMode(NewMode)
       
        switch( NewMode )
            case 'Ouverture'
                hParams.figR.Visible = 'off';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.TypeLabel.Visible = 'off';
                hParams.TypePopMenu.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.ChanLabel.Visible = 'off';
                hParams.ChanPopMenu.Visible = 'off';
                hParams.dFsFPB.Visible = 'off';
                hParams.GSRPB.Visible = 'off';
                hParams.CorrPB.Visible = 'off';
                
            case 'PreAna'
                hParams.figR.Visible = 'off';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                
                hParams.PreAPB.Visible = 'on';
                hParams.TypeLabel.Visible = 'off';
                hParams.TypePopMenu.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.ChanLabel.Visible = 'off';
                hParams.ChanPopMenu.Visible = 'off';
                hParams.dFsFPB.Visible = 'off';
                hParams.GSRPB.Visible = 'off';
                hParams.CorrPB.Visible = 'off';
                
            case 'SelectParams'
                hParams.figR.Visible = 'on';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.TypeLabel.Visible = 'on';
                hParams.TypePopMenu.Visible = 'on';
                hParams.PreALabel.Visible = 'off';
                hParams.ChanLabel.Visible = 'on';
                hParams.ChanPopMenu.Visible = 'on';
                hParams.dFsFPB.Visible = 'off';
                hParams.GSRPB.Visible = 'off';
                hParams.CorrPB.Visible = 'off';
                
            case 'RestingState'
                hParams.figR.Visible = 'on';
                hParams.figC.Visible = 'on';
                hParams.figT.Visible = 'on';
                
                hParams.PreAPB.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.TypeLabel.Visible = 'on';
                hParams.TypePopMenu.Visible = 'on';
                hParams.ChanLabel.Visible = 'on';
                hParams.ChanPopMenu.Visible = 'on';
                hParams.dFsFPB.Visible = 'on';
                hParams.GSRPB.Visible = 'on';
                hParams.CorrPB.Visible = 'on';
                
            otherwise
                hParams.figR.Visible = 'off';
                hParams.figC.Visible = 'off';
                hParams.figT.Visible = 'off';
                
                hParams.PreAPB.Visible = 'off';
                hParams.TypeLabel.Visible = 'off';
                hParams.TypePopMenu.Visible = 'off';
                hParams.PreALabel.Visible = 'off';
                hParams.ChanLabel.Visible = 'off';
                hParams.ChanPopMenu.Visible = 'off';
                hParams.dFsFPB.Visible = 'off';
                hParams.GSRPB.Visible = 'off';
                hParams.CorrPB.Visible = 'off';
        end
        
    end

    function DFsF(~,~,~)
        dims = size(Data);
        Data = reshape(Data, [], dims(3));
        f = fdesign.lowpass('N,F3dB', 4, 1/120, Infos.Freq);
        lpass = design(f,'butter');
        lpData = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(Data)'))';
        Data = Data./lpData;
        clear lpData;
        f = fdesign.lowpass('N,F3dB', 4, 2.5, Infos.Freq);
        lpass = design(f,'butter');
        Data = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(Data)'))';        
        
        Data = reshape(Data,dims);
        hParams.CI_Min_Edit.Value = 0.95;
        hParams.CI_Max_Edit.Value = 1.05;
        ChangeImage();
        hParams.dFsFPB.Enable = 'off';
        
        CorrMap();
        DecoursTemp();
    end

    function GSR(~,~,~)
        dims = size(Data);
        Data = reshape(Data, [], dims(3));
        Signal = mean(Data,1);
        Signal = Signal./mean(Signal);
        X = [ones(1,dims(3)); Signal];
        B = X'\Data';
        A = (X'*B)';
        Data = Data - A;
        Data = reshape(Data, dims);
        hParams.CI_Min_Edit.Value = -0.05;
        hParams.CI_Max_Edit.Value = 0.05;
        ChangeImage();
        hParams.GSRPB.Enable = 'off';
        
        CorrMap();
        DecoursTemp();
    end

    function CorrMap(~,~,~)
       cData = imresize(Data,[64 64]);
       cData = cData - mean(cData,3);
       cData = corr(reshape(cData,[],size(cData,3))');
       
       RefreshCorrMap();
    end

    function RefreshCorrMap(~,~,~)
        if( isempty(cData) )
            return;
        end
        
        Facteur = 64/256;
        Id = floor((currentPixel(1)-1)*Facteur)*64 + floor(currentPixel(2)*Facteur);
        imagesc(hParams.axC1, reshape(cData(Id,:),64,64),[-1 1]);
        axis(hParams.axC1, 'off', 'image');
    end

    function ChangePtPos(Obj, Evnt)
       Pos = round(Obj.Children(6).CurrentPoint);
       Pos = Pos(1,1:2);
       if( any(Pos < 1) | any(Pos > 255) )
           return;
       end
       currentPixel = Pos;
       ChangeImage();
       RefreshCorrMap();
       DecoursTemp();
    end

    function DecoursTemp()
        hold(hParams.axT1, 'off');
        plot(hParams.axT1, squeeze(Data(currentPixel(2),currentPixel(1),:)));
        hold(hParams.axT1, 'on');
        line(hParams.axT1,[hParams.CurrentImageSl.Value, hParams.CurrentImageSl.Value],...
            [hParams.CI_Min_Edit.Value, hParams.CI_Max_Edit.Value],...
            'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        ylim(hParams.axT1, [hParams.CI_Min_Edit.Value, hParams.CI_Max_Edit.Value]);
    end

    function NeFermePas(~,~,~)
        
    end

    function FermeTout(~,~,~)
        delete(hParams.figT);
        delete(hParams.figR);
        delete(hParams.figC);
        delete(hParams.figP);
    end
end