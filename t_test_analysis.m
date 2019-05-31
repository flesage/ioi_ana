
% enter group 1 datas paths
go = true;

channel_list = {'Flow','Green','HbO','HbR','HbT','Red','Yellow'};
name_g1 = inputdlg('Add first group dataset name');
fprintf(char(strcat('Add data paths for :',name_g1,'\n')));
nb_path_g1 = 0;
while(go)
    nb_path_g1 = nb_path_g1 + 1;
    paths_g1(nb_path_g1) = cellstr(uigetdir());
    for i = 1:length(channel_list)
        load(char(strcat(paths_g1(nb_path_g1),filesep,'bilateral_correlation_',channel_list(i),'.mat')));
        z_gr1_datas(:,i,nb_path_g1) = homologue_fisher(correlation);
    end
    fprintf('Path added! \n');
    answer = questdlg('Do you want to add another ROI?',...
        'Before continuing...',...
        'Yes','No','No');
    if(strcmp(answer,'No'))
        go = false;
    end
end

% enter group 2 datas paths

go = true;
name_g2 = inputdlg('Add second group dataset name');
fprintf(char(strcat('Add data paths for :',name_g2,'\n')));
nb_path_g2 = 0;
while(go)
    nb_path_g2 = nb_path_g2 + 1;
    paths_g2(nb_path_g2) = cellstr(uigetdir());
    for i = 1:length(channel_list)
        load(char(strcat(paths_g2(nb_path_g2),filesep,'bilateral_correlation_',channel_list(i),'.mat')));
        z_gr2_datas(:,i,nb_path_g2) = homologue_fisher(correlation);
    end
    fprintf('Path added! \n');
    answer = questdlg('Do you want to add another ROI?',...
        'Before continuing...',...
        'Yes','No','No');
    if(strcmp(answer,'No'))
        go = false;
    end
end

% t-test between groups
%%
roi_list = {'frontal','motor','cingulate','somato','retrospin','visual'};
answer = inputdlg('Add a file name');
save_path = cellstr(uigetdir());
complete_path = char(strcat(save_path,filesep,answer,'.txt'));
folder_save_path = char(strcat(save_path,filesep,answer,'_figures'));
if(~exist(folder_save_path, 'dir'))
    mkdir(folder_save_path);
end
fileID = fopen(complete_path,'w');
fprintf(fileID,'Group 1 "%s" datas paths \r\n',char(name_g1));
for i = 1:nb_path_g1
    fprintf(fileID,'%s \r\n',char(paths_g1(i)));
end
fprintf(fileID,'\nGroup 2 "%s" datas paths \n\r',char(name_g2));
for i = 1:nb_path_g2
    fprintf(fileID,'%s \r\n',char(paths_g2(i)));
end

fprintf(fileID,'\n=============== Two-sample t-tests for Fischer correlation ============== \r\n');
for i = 1:length(channel_list)
    fprintf(fileID,'\n-------------- Tests for %s channel -------------- \r\n',char(channel_list(i)));
    for j = 1:length(roi_list)
        fprintf(fileID,'\nTest for %s ROI\r\n',char(roi_list(j)));
        g1 = squeeze(z_gr1_datas(j,i,:));
        mean_g1 = mean(g1);
        var_g1 = var(g1);
        g2 = squeeze(z_gr2_datas(j,i,:));
        mean_g2 = mean(g2);
        var_g2 = var(g2);
        [h,p] = ttest2(g1,g2,'Vartype','unequal','Alpha',0.05);
        fprintf(fileID,'group:        %s        %s \r\n',char(name_g1),char(name_g2));
        fprintf(fileID,'mean:        %2.4f     %2.4f \r\n',mean_g1,mean_g2);
        fprintf(fileID,'variance:    %2.4f     %2.4f \r\n',var_g1,var_g2);
        fprintf(fileID,'number of data: %d          %d \r\n',nb_path_g1,nb_path_g2);
        fprintf(fileID,'P value: %2.4f \r\n',p);
        if(h)
            fprintf(fileID,'Decision: significant difference \r\n');
        else
            fprintf(fileID,'Decision: non-significant difference \r\n');
        end
        c = categorical({char(name_g1),char(name_g2)});
        e1 = sqrt(var_g1/nb_path_g1);
        e2 = sqrt(var_g2/nb_path_g2);
        fig = figure();
        bar(c,[mean_g1 mean_g2]);
        hold on;
        eb = errorbar(c,[mean_g1 mean_g2],[-e1 -e2],[e1 e2]); 
        eb.Color = [0 0 0];                            
        eb.LineStyle = 'none'; 
        xlabel('Groups');
        ylabel('Average');
        fig_save_path = char(strcat(save_path,filesep,answer,'_figures',filesep,channel_list(i),'_',roi_list(j)));
        print(fig,fig_save_path,'-djpeg');
        delete(fig);
    end
end

