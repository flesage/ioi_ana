function ConvertToMat(folderPath)

if( ~strcmp(folderPath(end), filesep) )
    folderPath = strcat(folderPath, filesep);
end

disp(['Starting .mat conversion for: ' folderPath]);
ChanList = dir([folderPath '*Chan.dat']);
InfList = dir([folderPath 'Data_*.mat']);
for indC = 1:size(ChanList,1)
    ChanName = ChanList(indC).name;
    fid = fopen([folderPath ChanName]);
    dat = fread(fid,inf,'single');
    idx = arrayfun(@(X) contains(InfList(X).name,['_' ChanName(1)],'IgnoreCase',true), 1:size(InfList,1)); 
    Infos = matfile([folderPath InfList(find(idx)).name]);
    dat = reshape(dat, Infos.datSize(1,1), Infos.datSize(1,2), []);
    Stim = Infos.Stim;
    
    outFName = [folderPath 'mat_'];
    switch(ChanName(1))
        case 'g'
            disp('Saving green channel.');
            outFName = strcat(outFName, 'Green.mat');
        case 'r'
            disp('Saving red channel.');
            outFName = strcat(outFName, 'Red.mat');
        case 'y'
            disp('Saving amber channel.');
            outFName = strcat(outFName, 'Yellow.mat');
        case 'f'
            disp('Saving fluo channel.');
            outFName = strcat(outFName, 'Fluo.mat');
    end
    F = matfile(outFName, 'Writable', true);
    
    idCond = nonzeros(unique(Stim));
    nbCond = length(idCond);
    for indS = 1:nbCond
        iStim = idCond(indS) == Stim;
        iS = find(diff(iStim)>0) + 1;
        iE = find(diff(iStim)<0) + 1;
        if( Stim(1) == idCond(indS) )
            iS = cat(1, 1, iS);
        end
        iL = min(iE - iS);
        Cond = zeros(size(dat,1),size(dat,2),iL,length(iS));
        for indR = 1:length(iS)
           Cond(:,:,:,indR) = dat(:,:,iS(indR) + (0:(iL-1))); 
        end
        eval(['F.Cond' int2str(idCond(indS)) ' = Cond;'])
    end 
end

disp('Mat conversion is done.')

end
