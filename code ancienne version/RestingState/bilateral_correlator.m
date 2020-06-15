clear all
warning off
%% Description
% This program allows to check correlation between ROIs in a dataset
% Update: 10-05-2019
% Author: Brent Fortin-Boily

% define ROIs and channels
ROI_list ={'left_frontal','left_motor','left_cingulate','left_somato','left_retrospin','left_visual'...
    ,'right_frontal','right_motor','right_cingulate','right_somato','right_retrospin','right_visual'};
channel_list = {'Flow','Green','HbO','HbR','HbT','Red','Yellow'};

add = true;
nbPath = 0;
% open loop until the user want to stop analyze datas
while(add)
    nbPath = nbPath + 1;
    fprintf('Add a data path \n');
    path(nbPath) = cellstr(uigetdir()); % getting data path
    
    % checking channel and ROIs are present in folder
    for ind_channel = 1:length(channel_list)
        fprintf('Checking file: %s \n',char(path(nbPath)));
        allGood = true;
        for i = 1:length(ROI_list)
            folder_path = char(strcat(path(nbPath),filesep,'image_correlation',filesep,ROI_list(i)));
            roi_file_path = char(strcat(path(nbPath),filesep,ROI_list(i),'.mat'));
            corr_file_path = char(strcat(folder_path,filesep,channel_list(ind_channel),'.mat'));
            if(isfolder(folder_path)&& isfile(roi_file_path) && isfile(corr_file_path))
                roi_exist(nbPath,i) = 1;
            else
                roi_exist(nbPath,i) = 0;
                warn = char(strcat('There is files missing for that ROI: ',ROI_list(i),', Channel: ', channel_list(ind_channel)));
                fprintf('%s \n',warn);
                allGood = false;
            end    
        end
        if(allGood)
            fprintf('No file missing! \n');
        end
        
        % bilateral correlation
        for i = 1:length(ROI_list)
            if(roi_exist(nbPath,i))
                folder_path = char(strcat(path(nbPath),filesep,'image_correlation',filesep,ROI_list(i)));
                roi_file_path = char(strcat(path(nbPath),filesep,ROI_list(i),'.mat'));
                corr_file_path = char(strcat(folder_path,filesep,channel_list(ind_channel),'.mat'));
                roi_selected = load(roi_file_path);
                roi_selected = roi_selected.ROIs{1,1}.mask;
                load(corr_file_path);
                for j = 1:length(ROI_list)
                    if(roi_exist(nbPath,j))
                        roicmp_file_path = char(strcat(path(nbPath),filesep,ROI_list(j),'.mat'));
                        roi_cmp = load(roicmp_file_path);
                        roi_cmp = roi_cmp.ROIs{1,1}.mask;
                        corr_cmp = corr_map(roi_cmp);
                        correlation(i,j) = mean(corr_cmp);
                    else
                        correlation(i,j) = nan;
                    end
                end
            else
                correlation(i,:) = nan;
            end
        end
        
        % saving in .mat and .xlsx format
        filename = char(strcat(path(nbPath),filesep,'bilateral_correlation','_',channel_list(ind_channel)));
        mat_filename = char(strcat(filename,'.mat'));
        xlsx_filename = char(strcat(filename,'.xlsx'));
        
        save(mat_filename,'correlation'); % saving in .mat
        
        % building table for .xlsx format
        reference=ROI_list';
        left_frontal = correlation(:,1);
        left_motor = correlation(:,2);
        left_cingulate = correlation(:,3);
        left_somato = correlation(:,4);
        left_retrospin = correlation(:,5);
        left_visual = correlation(:,6);
        right_frontal = correlation(:,7);
        right_motor = correlation(:,8);
        right_cingulate = correlation(:,9);
        right_somato = correlation(:,10);
        right_retrospin = correlation(:,11);
        right_visual = correlation(:,12);
        T = table(reference,left_frontal,left_motor,left_cingulate,left_somato,left_retrospin,left_visual,right_frontal,right_motor,right_cingulate,...
            right_somato,right_retrospin,right_visual);
        writetable(T,xlsx_filename); % saving in .xlsx
% Display bilateral correlation figure
ROI_displaylist ={'left frontal','left motor','left cingulate','left somato','left retrospin','left visual'...
    ,'right frontal','right motor','right cingulate','right somato','right retrospin','right visual'};
        step = 1:1:12;
        fig=figure('Units','normalized','position',[0 0 0.90 0.90]);
        imagesc(correlation,[-1,1]);
        set(gca,'XTick',step);
        set(gca,'YTick',step);
        xticklabels(ROI_displaylist);
        yticklabels(ROI_displaylist);
        xlabel('compared roi');
        ylabel('selected map');
        str_title = char(strcat(path(nbPath),'_',channel_list(ind_channel)));
        colormap jet
        colorbar
        print(fig,filename,'-djpeg');
        delete(fig);
    end

    % asking user to analyze a new dataset
    answer = questdlg('Do you want to add another ROI?',...
        'Before continuing...',...
        'Yes','No','No');
    if(strcmp(answer,'No'))
        add = false;
    end
end




