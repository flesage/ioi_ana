% T-test analysis et 2-way ANOVA script for bilateral correlation
% Author: Brent Fortin-Boily
% Updated : 24-06-2020

clear all;

% analysis type selection
answer = questdlg('Choose the analysis you want to do.', ...
	'Analysis', ...
	'T-test','2-way ANOVA','T-test');
if(strcmp(answer,'T-test'))
    disp('T-test selected.')
    t_test();
else
     disp('two-way ANOVA selected.')
    anova_test();
end


function t_test(~,~,~)
    channel_list = {'Flow','Green','HbO','HbR','HbT','Red','Yellow'};
    [indx,tf] = listdlg('PromptString','Select a or many channel.','ListString',channel_list);
    name_g1 = inputdlg('Add first group dataset name');
    fprintf(char(strcat('Add data path for :',name_g1,'\n')));
    main_path_g1 = uigetdir();
    datas_paths_g1 = dir(main_path_g1);
    nb_path_g1 = 0;
    for indi = 3:size(datas_paths_g1,1)
        nb_path_g1 = nb_path_g1 + 1;
        paths_g1(nb_path_g1) = strcat(cellstr(main_path_g1),filesep,datas_paths_g1(indi).name);
        nb_chan = 0;
        isGoodPath_g1(nb_path_g1) = 0;
        for i = 1:length(channel_list)
            pos = strfind(char(paths_g1(nb_path_g1)),filesep);
            temp_path = char(paths_g1(nb_path_g1));
            complete_path = char(strcat(temp_path,filesep,'bilateral_correlation_',channel_list(i),'_',temp_path(pos(end)+4:end),'.mat'));
            if(isfile(complete_path) && sum(ismember(indx,i)))
                isGoodPath_g1(nb_path_g1) = 1;
                nb_chan = nb_chan + 1;
                load(complete_path);
                bi_corr_gr1(:,:,nb_path_g1,nb_chan)= bi_corr;
            end
        end
    end

    fprintf('All paths added! \n');
    %%
    % enter group 2 datas paths

    name_g2 = inputdlg('Add second group dataset name');
    fprintf(char(strcat('Add data path for :',name_g2,'\n')));
    main_path_g2 = uigetdir();
    datas_paths_g2 = dir(main_path_g2);
    nb_path_g2 = 0;

    for indi = 3:size(datas_paths_g2,1)
        nb_path_g2 = nb_path_g2 + 1;
        paths_g2(nb_path_g2) = strcat(cellstr(main_path_g2),filesep,datas_paths_g2(indi).name);
        nb_chan = 0;
        isGoodPath_g2(nb_path_g2) = 0;
        for i = 1:length(channel_list)
            pos = strfind(char(paths_g2(nb_path_g2)),filesep);
            temp_path = char(paths_g2(nb_path_g2));
            complete_path = char(strcat(temp_path,filesep,'bilateral_correlation_',channel_list(i),'_',temp_path(pos(end)+4:end),'.mat'));
            if(isfile(complete_path) && sum(ismember(indx,i)))
                isGoodPath_g2(nb_path_g2) = 1;
                nb_chan = nb_chan + 1;
                load(complete_path);
                bi_corr_gr2(:,:,nb_path_g2,nb_chan)= bi_corr;
            end
        end 
    end


        fprintf('All paths added! \n');

    % t-test between groups
    %%
    roi_list = {'frontal','motor','cingulate','somato','retrospin','visual'};

    % getting file name and folder path for saving
    answer = inputdlg('Add a file name');
    save_path = cellstr(uigetdir());
    complete_path = char(strcat(save_path,filesep,answer,'.txt'));
    folder_save_path = char(strcat(save_path,filesep,answer,'_figures'));
    % create path
    if(~exist(folder_save_path, 'dir'))
        mkdir(folder_save_path);
    end
    %saving mean bicorrelation map for each group
    disp('Starting analysis...');
        ROI_displaylist ={'left frontal','left motor','left cingulate','left somato','left retrospin','left visual'...
                ,'right frontal','right motor','right cingulate','right somato','right retrospin','right visual'};
        step = 1:1:12;
        nb_chan = 0;
    for i = 1:length(channel_list)
        if(sum(ismember(indx,i)))
            nb_chan = nb_chan + 1;
            temp = squeeze(bi_corr_gr1(:,:,:,nb_chan));
            mean_temp = mean(temp,3);
            fig5 = figure('Units','normalized','position',[0 0 0.90 0.90],'Visible','off');
            imagesc(mean_temp,[-1,1]);
            set(gca,'XTick',step);
            set(gca,'YTick',step);
            xticklabels(ROI_displaylist);
            yticklabels(ROI_displaylist);
            xlabel('compared roi');
            ylabel('selected map');
            str_title = char(strcat('mean_gr_,',name_g1,'_',channel_list(i)));
            title(str_title);
            colormap jet
            colorbar
            filename = char(strcat(folder_save_path,filesep,'mean_bicorr_gr_',name_g1,'_',channel_list(i)));
            print(fig5,filename,'-djpeg');
            delete(fig5);  
        end
    end

    nb_chan = 0;
    for i = 1:length(channel_list)
        if(sum(ismember(indx,i)))
            nb_chan = nb_chan + 1;
            temp = squeeze(bi_corr_gr2(:,:,:,nb_chan));
            mean_temp = mean(temp,3);
            fig5 = figure('Units','normalized','position',[0 0 0.90 0.90],'Visible','off');
            imagesc(mean_temp,[-1,1]);
            set(gca,'XTick',step);
            set(gca,'YTick',step);
            xticklabels(ROI_displaylist);
            yticklabels(ROI_displaylist);
            xlabel('compared roi');
            ylabel('selected map');
            str_title = char(strcat('mean_gr_,',name_g2,'_',channel_list(i)));
            title(str_title);
            colormap jet
            colorbar
            filename = char(strcat(folder_save_path,filesep,'mean_bicorr_gr_',name_g2,'_',channel_list(i)));
            print(fig5,filename,'-djpeg');
            delete(fig5);
        end
    end
    % making .txt file for saving results
    fileID = fopen(complete_path,'w');

    % saving results in .txt file
    fprintf(fileID,'Group 1 "%s" datas paths \r\n',char(name_g1));
    for i = 1:nb_path_g1
        if(isGoodPath_g1(i))
            fprintf(fileID,'%s \r\n',char(paths_g1(i)));
        end
    end
    fprintf(fileID,'\nGroup 2 "%s" datas paths \n\r',char(name_g2));
    for i = 1:nb_path_g2
        if(isGoodPath_g2(i))
            fprintf(fileID,'%s \r\n',char(paths_g2(i)));
        end
    end

    % bi correlation, figure plots and .txt file save.
    fprintf(fileID,'\n=============== Two-sample t-tests for bicorrelation ============== \r\n');
    nb_chan = 0;
    for i = 1:length(channel_list)
        if(sum(ismember(indx,i)))
            nb_chan = nb_chan + 1;
            fprintf(fileID,'\n-------------- Tests for %s channel -------------- \r\n',char(channel_list(i)));
            for j = 1:length(ROI_displaylist)-1
                fprintf(fileID,'\n Tests for reference ROI: %s =========\r\n',char(ROI_displaylist(j)));
                for k = j+1:length(ROI_displaylist)
                    fprintf(fileID,'\n Test with ROI: %s \r\n',char(ROI_displaylist(k)));
                    g1 = squeeze(bi_corr_gr1(j,k,:,nb_chan));
                    fprintf(fileID,'\n Values for group 1 \r\n');
                    for ind = 1:length(g1)
                        fprintf(fileID,'%2.4f \n',g1(ind));
                    end
                    mean_g1 = mean(g1);
                    var_g1 = var(g1);
                    g2 = squeeze(bi_corr_gr2(j,k,:,nb_chan));
                    fprintf(fileID,'\n Values for group 2 \r\n');
                    for ind = 1:length(g2)
                        fprintf(fileID,'%2.4f \n',g2(ind));
                    end
                    mean_g2 = mean(g2);
                    var_g2 = var(g2);
                    [h,p] = ttest2(g1,g2,'Vartype','unequal','Alpha',0.05);
                    fprintf(fileID,'group:        %s        %s \r\n',char(name_g1),char(name_g2));
                    fprintf(fileID,'mean:        %2.4f     %2.4f \r\n',mean_g1,mean_g2);
                    fprintf(fileID,'variance:    %2.4f     %2.4f \r\n',var_g1,var_g2);
                    fprintf(fileID,'number of data: %d          %d \r\n',nb_path_g1,nb_path_g2);
                    fprintf(fileID,'P value: %2.4f \r\n',p);
                    if(~isnan(h))
                        if(h)
                            fprintf(fileID,'Decision: significant difference \r\n');
                        else
                            fprintf(fileID,'Decision: non-significant difference \r\n');
                        end
                    else
                         fprintf(fileID,'Decision: not a number \r\n');
                    end
                    % Making compare group figure
                    c = categorical({char(name_g1),char(name_g2)});
                    e1 = sqrt(var_g1/nb_path_g1);
                    e2 = sqrt(var_g2/nb_path_g2);
                    fig = figure('visible','off');
                    bar(c,[mean_g1 mean_g2]);
                    hold on;
                    eb = errorbar(c,[mean_g1 mean_g2],[-e1 -e2],[e1 e2]); 
                    eb.Color = [0 0 0];                            
                    eb.LineStyle = 'none'; 
                    xlabel('Groups');
                    ylabel('Average');
                    fig_save_path = char(strcat(save_path,filesep,answer,'_figures',filesep,channel_list(i),'_',ROI_displaylist(j),'_',ROI_displaylist(k)));
                    print(fig,fig_save_path,'-djpeg');
                    delete(fig);
                end
            end
        end
    end
    disp('Analysis done!');
    fclose(fileID);
end

function anova_test(~,~,~)
channel_list = {'Flow','Green','HbO','HbR','HbT','Red','Yellow'};
    [indx,tf] = listdlg('PromptString','Select a or many channel.','ListString',channel_list);
    name_g1 = inputdlg('Add first group dataset name');
    %enter group 1 paths
    fprintf(char(strcat('Add data path for :',name_g1,'\n')));
    main_path_g1 = uigetdir();
    datas_paths_g1 = dir(main_path_g1);
    nb_path_g1 = 0;
    for indi = 3:size(datas_paths_g1,1)
        nb_path_g1 = nb_path_g1 + 1;
        paths_g1(nb_path_g1) = strcat(cellstr(main_path_g1),filesep,datas_paths_g1(indi).name);
        nb_chan = 0;
        isGoodPath_g1(nb_path_g1) = 0;
        for i = 1:length(channel_list)
            pos = strfind(char(paths_g1(nb_path_g1)),filesep);
            temp_path = char(paths_g1(nb_path_g1));
            complete_path = char(strcat(temp_path,filesep,'bilateral_correlation_',channel_list(i),'_',temp_path(pos(end)+4:end),'.mat'));
            if(isfile(complete_path) && sum(ismember(indx,i)))
                isGoodPath_g1(nb_path_g1) = 1;
                nb_chan = nb_chan + 1;
                load(complete_path);
                bi_corr_gr1(:,:,nb_path_g1,nb_chan)= bi_corr;
            end
        end
    end

    fprintf('All paths added! \n');
    %%
    % enter group 2 datas paths

    name_g2 = inputdlg('Add second group dataset name');
    fprintf(char(strcat('Add data path for :',name_g2,'\n')));
    main_path_g2 = uigetdir();
    datas_paths_g2 = dir(main_path_g2);
    nb_path_g2 = 0;

    for indi = 3:size(datas_paths_g2,1)
        nb_path_g2 = nb_path_g2 + 1;
        paths_g2(nb_path_g2) = strcat(cellstr(main_path_g2),filesep,datas_paths_g2(indi).name);
        nb_chan = 0;
        isGoodPath_g2(nb_path_g2) = 0;
        for i = 1:length(channel_list)
            pos = strfind(char(paths_g2(nb_path_g2)),filesep);
            temp_path = char(paths_g2(nb_path_g2));
            complete_path = char(strcat(temp_path,filesep,'bilateral_correlation_',channel_list(i),'_',temp_path(pos(end)+4:end),'.mat'));
            if(isfile(complete_path) && sum(ismember(indx,i)))
                isGoodPath_g2(nb_path_g2) = 1;
                nb_chan = nb_chan + 1;
                load(complete_path);
                bi_corr_gr2(:,:,nb_path_g2,nb_chan)= bi_corr;
            end
        end 
    end
        fprintf('All paths added! \n');
        
    %%
    % enter group 3 datas paths

    name_g3 = inputdlg('Add third group dataset name');
    fprintf(char(strcat('Add data path for :',name_g3,'\n')));
    main_path_g3 = uigetdir();
    datas_paths_g3 = dir(main_path_g3);
    nb_path_g3 = 0;

    for indi = 3:size(datas_paths_g3,1)
        nb_path_g3 = nb_path_g3 + 1;
        paths_g3(nb_path_g3) = strcat(cellstr(main_path_g3),filesep,datas_paths_g3(indi).name);
        nb_chan = 0;
        isGoodPath_g3(nb_path_g3) = 0;
        for i = 1:length(channel_list)
            pos = strfind(char(paths_g3(nb_path_g3)),filesep);
            temp_path = char(paths_g3(nb_path_g3));
            complete_path = char(strcat(temp_path,filesep,'bilateral_correlation_',channel_list(i),'_',temp_path(pos(end)+4:end),'.mat'));
            if(isfile(complete_path) && sum(ismember(indx,i)))
                isGoodPath_g3(nb_path_g3) = 1;
                nb_chan = nb_chan + 1;
                load(complete_path);
                bi_corr_gr3(:,:,nb_path_g3,nb_chan)= bi_corr;
            end
        end 
    end
        fprintf('All paths added! \n');
        
       %%
    % enter group 4 datas paths

    name_g4 = inputdlg('Add fourth group dataset name');
    fprintf(char(strcat('Add data path for :',name_g4,'\n')));
    main_path_g4 = uigetdir();
    datas_paths_g4 = dir(main_path_g4);
    nb_path_g4 = 0;

    for indi = 3:size(datas_paths_g4,1)
        nb_path_g4 = nb_path_g4 + 1;
        paths_g4(nb_path_g4) = strcat(cellstr(main_path_g4),filesep,datas_paths_g4(indi).name);
        nb_chan = 0;
        isGoodPath_g4(nb_path_g4) = 0;
        for i = 1:length(channel_list)
            pos = strfind(char(paths_g4(nb_path_g4)),filesep);
            temp_path = char(paths_g4(nb_path_g4));
            complete_path = char(strcat(temp_path,filesep,'bilateral_correlation_',channel_list(i),'_',temp_path(pos(end)+4:end),'.mat'));
            if(isfile(complete_path) && sum(ismember(indx,i)))
                isGoodPath_g4(nb_path_g4) = 1;
                nb_chan = nb_chan + 1;
                load(complete_path);
                bi_corr_gr4(:,:,nb_path_g4,nb_chan)= bi_corr;
            end
        end 
    end
        fprintf('All paths added! \n');

    % 2-way Anova test
    %%
    roi_list = {'frontal','motor','cingulate','somato','retrospin','visual'};
    ROI_displaylist ={'left frontal','left motor','left cingulate','left somato','left retrospin','left visual'...
                 ,'right frontal','right motor','right cingulate','right somato','right retrospin','right visual'};

    % getting file name and folder path for saving
    answer = inputdlg('Add a file name');
    save_path = cellstr(uigetdir());
    complete_path = char(strcat(save_path,filesep,answer,'.txt'));
    folder_save_path = char(strcat(save_path,filesep,answer,'_figures'));
    % create path
    if(~exist(folder_save_path, 'dir'))
        mkdir(folder_save_path);
    end
%     %saving mean bicorrelation map for each group
%     disp('Starting analysis...');
%         ROI_displaylist ={'left frontal','left motor','left cingulate','left somato','left retrospin','left visual'...
%                 ,'right frontal','right motor','right cingulate','right somato','right retrospin','right visual'};
%         step = 1:1:12;
%         nb_chan = 0;
%     for i = 1:length(channel_list)
%         if(sum(ismember(indx,i)))
%             nb_chan = nb_chan + 1;
%             temp = squeeze(bi_corr_gr1(:,:,:,nb_chan));
%             mean_temp = mean(temp,3);
%             fig5 = figure('Units','normalized','position',[0 0 0.90 0.90],'Visible','off');
%             imagesc(mean_temp,[-1,1]);
%             set(gca,'XTick',step);
%             set(gca,'YTick',step);
%             xticklabels(ROI_displaylist);
%             yticklabels(ROI_displaylist);
%             xlabel('compared roi');
%             ylabel('selected map');
%             str_title = char(strcat('mean_gr_,',name_g1,'_',channel_list(i)));
%             title(str_title);
%             colormap jet
%             colorbar
%             filename = char(strcat(folder_save_path,filesep,'mean_bicorr_gr_',name_g1,'_',channel_list(i)));
%             print(fig5,filename,'-djpeg');
%             delete(fig5);  
%         end
%     end
% 
%     nb_chan = 0;
%     for i = 1:length(channel_list)
%         if(sum(ismember(indx,i)))
%             nb_chan = nb_chan + 1;
%             temp = squeeze(bi_corr_gr2(:,:,:,nb_chan));
%             mean_temp = mean(temp,3);
%             fig5 = figure('Units','normalized','position',[0 0 0.90 0.90],'Visible','off');
%             imagesc(mean_temp,[-1,1]);
%             set(gca,'XTick',step);
%             set(gca,'YTick',step);
%             xticklabels(ROI_displaylist);
%             yticklabels(ROI_displaylist);
%             xlabel('compared roi');
%             ylabel('selected map');
%             str_title = char(strcat('mean_gr_,',name_g2,'_',channel_list(i)));
%             title(str_title);
%             colormap jet
%             colorbar
%             filename = char(strcat(folder_save_path,filesep,'mean_bicorr_gr_',name_g2,'_',channel_list(i)));
%             print(fig5,filename,'-djpeg');
%             delete(fig5);
%         end
%     end
    % making .txt file for saving results
    fileID = fopen(complete_path,'w');

    % saving results in .txt file
    fprintf(fileID,'Group 1 "%s" datas paths \r\n',char(name_g1));
    for i = 1:nb_path_g1
        if(isGoodPath_g1(i))
            fprintf(fileID,'%s \r\n',char(paths_g1(i)));
        end
    end
    fprintf(fileID,'\nGroup 2 "%s" datas paths \n\r',char(name_g2));
    for i = 1:nb_path_g2
        if(isGoodPath_g2(i))
            fprintf(fileID,'%s \r\n',char(paths_g2(i)));
        end
    end
    
    fprintf(fileID,'\nGroup 3 "%s" datas paths \n\r',char(name_g3));
    for i = 1:nb_path_g3
        if(isGoodPath_g3(i))
            fprintf(fileID,'%s \r\n',char(paths_g3(i)));
        end
    end
    
    fprintf(fileID,'\nGroup 4 "%s" datas paths \n\r',char(name_g4));
    for i = 1:nb_path_g4
        if(isGoodPath_g4(i))
            fprintf(fileID,'%s \r\n',char(paths_g4(i)));
        end
    end
    
    % bi correlation, figure plots and .txt file save.
    disp('Processing data...')
    fprintf(fileID,'\n=============== Two-way Anova for bicorrelation ============== \r\n');
    nb_chan = 0;
    for i = 1:length(channel_list)
        if(sum(ismember(indx,i)))
            nb_chan = nb_chan + 1;
            fprintf(fileID,'\n-------------- Tests for %s channel -------------- \r\n',char(channel_list(i)));
            for j = 1:length(ROI_displaylist)-1
                fprintf(fileID,'\n Tests for reference ROI: %s =========\r\n',char(ROI_displaylist(j)));
                for k = j+1:length(ROI_displaylist)
                    fprintf(fileID,'\n Test with ROI: %s \r\n',char(ROI_displaylist(k)));
                    % group 1
                    g1 = squeeze(bi_corr_gr1(j,k,:,nb_chan));
                    fprintf(fileID,'\n Values for group 1 \r\n');
                    for ind = 1:length(g1)
                        fprintf(fileID,'%2.4f \n',g1(ind));
                    end
                    mean_g1 = mean(g1);
                    var_g1 = var(g1);
                    
                    % group 2
                    g2 = squeeze(bi_corr_gr2(j,k,:,nb_chan));
                    fprintf(fileID,'\n Values for group 2 \r\n');
                    for ind = 1:length(g2)
                        fprintf(fileID,'%2.4f \n',g2(ind));
                    end
                    mean_g2 = mean(g2);
                    var_g2 = var(g2);
                    % group 3
                    g3 = squeeze(bi_corr_gr3(j,k,:,nb_chan));
                    fprintf(fileID,'\n Values for group 3 \r\n');
                    for ind = 1:length(g3)
                        fprintf(fileID,'%2.4f \n',g3(ind));
                    end
                    mean_g3 = mean(g3);
                    var_g3 = var(g3);
                    % group 4
                    g4 = squeeze(bi_corr_gr4(j,k,:,nb_chan));
                    fprintf(fileID,'\n Values for group 4 \r\n');
                    for ind = 1:length(g4)
                        fprintf(fileID,'%2.4f \n',g4(ind));
                    end
                    mean_g4 = mean(g4);
                    var_g4 = var(g4);
                    
                    % Evaluate Anova and print it
                    matrix = [g1, g2, g3, g4]; % Important!! Groups must be the same size!!!! ====
                    [p,tbl] = anova2(matrix, 1,'off');
                    fprintf(fileID,'%s \t\t %s \t\t %s \t %s \t\t %s \t\t %s \n -------------------------------------------------------------------------\n',...
                        tbl{1,1},tbl{1,2},tbl{1,3},tbl{1,4},tbl{1,5},tbl{1,6});
                    fprintf(fileID,'%s \t %2.4f \t %d \t %2.4f \t %2.4f \t %2.4f \n',...
                        tbl{2,1},tbl{2,2},tbl{2,3},tbl{2,4},tbl{2,5},tbl{2,6});
                    fprintf(fileID,'%s \t\t %2.4f \t %d \t %2.4f \t %2.4f \t %2.4f \n',...
                        tbl{3,1},tbl{3,2},tbl{3,3},tbl{3,4},tbl{3,5},tbl{3,6});
                    fprintf(fileID,'%s \t\t %2.4f \t %d \t %2.4f \n',...
                        tbl{4,1},tbl{4,2},tbl{4,3},tbl{4,4});
                    fprintf(fileID,'%s \t\t %d \t %d \n',tbl{5,1},tbl{5,2},tbl{5,3});
                    % Making compare group figure
                    c = categorical({char(name_g1),char(name_g2),char(name_g3),char(name_g4)});
                    e1 = sqrt(var_g1/nb_path_g1);
                    e2 = sqrt(var_g2/nb_path_g2);
                    e3 = sqrt(var_g3/nb_path_g3);
                    e4 = sqrt(var_g4/nb_path_g4);
                    fig = figure('visible','off');
                    bar(c,[mean_g1 mean_g2 mean_g3 mean_g4]);
                    hold on;
                    eb = errorbar(c,[mean_g1 mean_g2 mean_g3 mean_g4],[-e1 -e2 -e3 -e4],[e1 e2 e3 e4]); 
                    eb.Color = [0 0 0];                            
                    eb.LineStyle = 'none'; 
                    xlabel('Groups');
                    ylabel('Average');
                    fig_save_path = char(strcat(save_path,filesep,answer,'_figures',filesep,channel_list(i),'_',ROI_displaylist(j),'_',ROI_displaylist(k)));
                    print(fig,fig_save_path,'-djpeg');
                    delete(fig);
                end
            end
        end
    end
    fclose(fileID);
    disp('Processing done!')
end