function out = ReadInfoFile(FolderPath)

if( ~strcmp(FolderPath(end), filesep) )
    FolderPath = strcat(FolderPath, filesep);
end

fid = fopen([FolderPath 'info.txt']);
bCamIndex = 0;
fgetl(fid);
indS = 1;
while ~feof(fid)
    tline = fgetl(fid);
%     disp(tline)
    Pos = strfind(tline, ':');
    % Look fot "tabs" if there is "Events" info table in "info.txt" file:
    if isempty(Pos) && ~isempty(tline)
        Pos = regexp(tline, '\t');
    end
    tline(1:Pos) = regexprep(tline(1:Pos), ' ', '_');
    if( length(Pos) == 1 )
        Param = tline(1:(Pos-1));
        Value = (tline((Pos+2):end));
        % Skip empty values:
        if isempty(Value) 
            continue
        end
            
        if( startsWith(Param, 'Illumination') )
            if( endsWith(Param, 'CameraIdx') )
                bCamIndex = 1;
                eval(['out.' Param(1:13) '.CamIdx = ' Value ';']);
            elseif( endsWith(Param, 'FrameIdx') )
                eval(['out.' Param(1:13) '.FrameIdx = ' Value ';']);
                bCamIndex = 1;
            elseif( regexp(Param, '[0-9]') )
                eval(['out.' Param(1:13) '.ID  = ' Param(13:end) ';']);
                eval(['out.' Param(1:13) '.Color  =  '''  Value  ''';']);
            end
        elseif( all(isstrprop(erase(Value, ' '), 'digit')) )% Detect arrays of integers.
            eval(['out.' Param ' = [' Value '];']);            
        elseif( isnan(str2double(Value)) ) 
            eval(['out.' Param ' = ''' Value ''';']);
        else
            eval(['out.' Param ' = ' Value ';']);
        end
    elseif( length(Pos) == 3 )  %Digital stim
        %PosEnd = strfind(tline, ',');
        if startsWith(tline, 'id', 'IgnoreCase', true) % Skip table header line.
            continue
        end
        str = regexp(tline,'\t', 'split');               
        eval(['out.Stim' int2str(indS) '.ID = ' str{1} ';']);
        eval(['out.Stim' int2str(indS) '.name = ''' str{2} ''';']);
        eval(['out.Stim' int2str(indS) '.code = ' str{3} ';']);
        eval(['out.Stim' int2str(indS) '.Duration = ' str{4} ';']);
        indS = indS + 1;
    elseif( length(Pos) == 4 )  %Digital stim with optogen
        PosEnd = strfind(tline, ',');
        StimName = tline((Pos(1)+2):(PosEnd(1)-1));
        StimCode = tline((Pos(2)+2):(PosEnd(2)-1));
        StimDuration = tline((Pos(3)+2):(PosEnd(3)-1));
        eval(['out.Stim' int2str(indS) '.name = ''' StimName ''';']);
        eval(['out.Stim' int2str(indS) '.code = ' StimCode ';']);
        eval(['out.Stim' int2str(indS) '.Duration = ' StimDuration ';']);
        indS = indS + 1;        
    end
    
    if( bCamIndex )
        out.MultiCam = 1;
    else
        out.MultiCam = 0;
    end
end