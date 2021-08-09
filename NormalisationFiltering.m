function OutDat = NormalisationFiltering(FolderData, FileData, lowFreq, highFreq, bDivide)

if( ~strcmp(FolderData(end),filesep) )
    FolderData = strcat(FolderData, filesep);
end

fprintf('Opening: %s \n', [FolderData FileData '.dat']);
Infos = matfile([FolderData FileData '.mat']);
fid = fopen([FolderData FileData '.dat']);
OutDat = fread(fid, inf, '*single');
OutDat = reshape(OutDat, Infos.datSize(1,1), Infos.datSize(1,2), []);


% Temporal filtering butterworth
f = fdesign.lowpass('N,F3dB', 4, lowFreq, Infos.Freq); %Fluo lower Freq
lpass = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, highFreq, Infos.Freq);   %Fluo Higher Freq
hpass = design(f,'butter');

dims = size(OutDat);
Hd = zeros(dims,'single');
PrcLims = round(linspace(1, dims(1), 11));
fprintf('Progress: ');
if( bDivide )
    for ind = 1:dims(1)
        Hd(ind,:,:) = reshape(single(filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(squeeze(OutDat(ind,:,:))'))'),dims(2),dims(3));
        Hd(ind,:,:) = squeeze(Hd(ind,:,:))./...
            reshape(single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(squeeze(OutDat(ind,:,:))'))'),dims(2),dims(3));
        if( any(ind == PrcLims) )
            idx = find(ind == PrcLims);
            fprintf('%d%%..', 10*(idx-1));
        end
    end
else
    for ind = 1:dims(1)
        Hd(ind,:,:) = reshape(single(filtfilt(hpass.sosMatrix, hpass.ScaleValues, double(squeeze(OutDat(ind,:,:))'))'),dims(2),dims(3));
        Hd(ind,:,:) = squeeze(Hd(ind,:,:))-...
            reshape(single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(squeeze(OutDat(ind,:,:))'))'),dims(2),dims(3));
        if( any(ind == PrcLims) )
            idx = find(ind == PrcLims);
            fprintf('%d%%..', 10*(idx-1));
        end
    end
end
fprintf('\n');    
OutDat = Hd;
end