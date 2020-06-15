% T-test analysis script for bilateral correlation
% Author: Brent Fortin-Boily
% Updated : 10-11-2019

clear all;

% enter group 1 datas paths
channel_list = {'Flow','Green','HbO','HbR','HbT','Red','Yellow'};
name_g1 = inputdlg('Add first group dataset name');
fprintf(char(strcat('Add data path for :',name_g1,'\n')));
main_path_g1 = uigetdir();
datas_paths_g1 = dir(main_path_g1);
nb_path_g1 = 0;

for indi = 3:size(datas_paths_g1,1)
    nb_path_g1 = nb_path_g1 + 1;
    paths_g1(nb_path_g1) = strcat(cellstr(main_path_g1),filesep,datas_paths_g1(indi).name);
    for i = 1:length(channel_list)
        pos = strfind(char(paths_g1(nb_path_g1)),filesep);
        temp_path = char(paths_g1(nb_path_g1));
        complete_path = char(strcat(paths_g1(nb_path_g1),filesep,'bc_',temp_path(pos(end)+1:end)...
            ,filesep,'bilateral_correlation_',channel_list(i),'_',temp_path(pos(end)+1:end),'.mat'));
        load(complete_path);
        bi_corr_gr1(:,:,nb_path_g1,i)= bi_corr;
        % z_gr1_datas(:,:,nb_path_g1,i) = homologue_fisher(bi_corr);
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
    for i = 1:length(channel_list)
        pos = strfind(char(paths_g2(nb_path_g2)),filesep);
        temp_path = char(paths_g2(nb_path_g2));
        load(char(strcat(paths_g2(nb_path_g2),filesep,'bc_',temp_path(pos(end)+1:end)...
            ,filesep,'bilateral_correlation_',channel_list(i),'_',temp_path(pos(end)+1:end),'.mat')));
        bi_corr_gr2(:,:,nb_path_g2,i)= bi_corr;
%         z_gr2_datas(:,:,nb_path_g2,i) = homologue_fisher(bi_corr);
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
    ROI_displaylist ={'left frontal','left motor','left cingulate','left somato','left retrospin','left visual'...
            ,'right frontal','right motor','right cingulate','right somato','right retrospin','right visual'};
    step = 1:1:12;
for i = 1:length(channel_list)
    temp = squeeze(bi_corr_gr1(:,:,:,i));
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

for i = 1:length(channel_list)
    temp = squeeze(bi_corr_gr2(:,:,:,i));
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
% making .txt file for saving results
fileID = fopen(complete_path,'w');

% saving results in .txt file
fprintf(fileID,'Group 1 "%s" datas paths \r\n',char(name_g1));
for i = 1:nb_path_g1
    fprintf(fileID,'%s \r\n',char(paths_g1(i)));
end
fprintf(fileID,'\nGroup 2 "%s" datas paths \n\r',char(name_g2));
for i = 1:nb_path_g2
    fprintf(fileID,'%s \r\n',char(paths_g2(i)));
end

% bi correlation, figure plots and .txt file save.
fprintf(fileID,'\n=============== Two-sample t-tests for bicorrelation ============== \r\n');
for i = 1:length(channel_list)
    fprintf(fileID,'\n-------------- Tests for %s channel -------------- \r\n',char(channel_list(i)));
    for j = 1:length(roi_list)-1
        fprintf(fileID,'\n Tests for reference ROI: %s \r\n',char(roi_list(j)));
        for k = j+1:length(roi_list)
            fprintf(fileID,'\n Test with ROI: %s \r\n',char(roi_list(k)));
            g1 = squeeze(bi_corr_gr1(j,k,:,i));
            mean_g1 = mean(g1);
            var_g1 = var(g1);
            g2 = squeeze(bi_corr_gr2(j,k,:,i));
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
            fig_save_path = char(strcat(save_path,filesep,answer,'_figures',filesep,channel_list(i),'_',roi_list(j),'_',roi_list(k)));
            print(fig,fig_save_path,'-djpeg');
            delete(fig);
        end
    end
end
fclose(fileID);
