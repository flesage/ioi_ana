function out = ReadInfoFile(FolderPath)

if( ~strcmp(FolderPath(end), filesep) )
    FolderPath = strcat(FolderPath, filesep);
end

fid = fopen([FolderPath 'info.txt']);
fgetl(fid);
indS = 1;
while ~feof(fid)
    tline = fgetl(fid);
%    disp(tline)
    Pos = strfind(tline, ':');
    tline(1:Pos) = regexprep(tline(1:Pos), ' ', '_');
    if( length(Pos) == 1 )
        Param = tline(1:(Pos-1));
        Value = (tline((Pos+2):end));
        if( contains(Param, 'Illumination')||contains(Param,'Name') )
            eval(['out.' Param ' = ''' Value ''';']);
        else
            eval(['out.' Param ' = ' Value ';']);
        end
    elseif( length(Pos) == 3 )
        PosEnd = strfind(tline, ',');
        StimName = tline((Pos(1)+2):(PosEnd(1)-1));
        StimCode = tline((Pos(2)+2):(PosEnd(2)-1));
        StimDuration = tline((Pos(3)+2):end);
        eval(['out.Stim' int2str(indS) '.name = ''' StimName ''';']);
        eval(['out.Stim' int2str(indS) '.code = ' StimCode ';']);
        eval(['out.Stim' int2str(indS) '.Duration = ' StimDuration ';']);
        indS = indS + 1;
    end
end