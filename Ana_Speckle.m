function out = Ana_Speckle(FolderName)

AcqInfoStream = readtable([FolderName filesep 'info.txt'],...
    'Delimiter',':','ReadVariableNames',false, 'ReadRowNames',true);

tExposure = AcqInfoStream{'ExposureMsec',1};
if( iscell(tExposure) )
    tExposure = str2double(cell2mat(tExposure));
end
tExposure = tExposure/1000.;

fprintf('Opening files.\n');
% Parameters
speckle_window_size = 5;
speckle_int_time = tExposure;

FileList = dir([FolderName filesep 'Data_speckle.mat']);
if( isempty(FileList) )
    disp(['No speckle data files found in ' FolderName ' Folder.']);
    disp('Speckle Analysis will not run');
    return;
else
    Iptr = matfile([FolderName filesep 'Data_speckle.mat']);
    nx = Iptr.datSize(1,1);
    ny = Iptr.datSize(1,2);
    nt = Iptr.datLength;
    tFreq = Iptr.Freq;

    Dptr = memmapfile(Iptr.datFile, 'Format', 'single');
end

% Need to convert to contrast
fprintf('Flow Conversion:\n');
speckle_window = ones(speckle_window_size);
OPTIONS.GPU = 0;
OPTIONS.Power2Flag = 0;
OPTIONS.Brep = 0;
dat = zeros(nt - 1, nx, ny, 'single');
prcflg = linspace(1, nt-1, 11); indP = 2;
for i3 = 1:nt-1
    if( i3 >= prcflg(indP) )
         fprintf('%d%%...', 10*(indP-1));
         indP = indP+1;
     end
    tmp_laser = Dptr.Data(((i3-1)*nx*ny + 1):(i3*nx*ny));
    tmp_laser = reshape(tmp_laser, nx, []);  
    std_laser=stdfilt(tmp_laser,speckle_window);
    mean_laser = convnfft(tmp_laser,speckle_window,'same',1:2,OPTIONS)/speckle_window_size^2;
    contrast=std_laser./mean_laser;
    dat(i3, :, :) = single(private_flow_from_contrast(contrast,speckle_int_time));
end
clear tmp_laser std_laser contrast mean_laser;

fprintf('\nFiltering:\n');
f = fdesign.lowpass('N,F3dB', 4, 1, tFreq);
ph = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, 1/120, tFreq);
pb = design(f,'butter');
nBlocs = round((nt-1)*ny*nx/(2^25));
marks = round(linspace(0, double(ny*nx), nBlocs));
if( length(marks) > 12 )
    prcflg = linspace(1, length(marks)-1, 11);
else
    prcflg = 1:length(marks);
end
indP = 1;
for ind = 1:(length(marks) - 1)
     if( ind >= prcflg(indP) )
         fprintf('%d%%...', round(100*prcflg(indP)/length(marks)));
         indP = indP+1;
     end
     pF = dat(:, (marks(ind)+1):marks(ind+1) );
 
     S_ioi= filtfilt(ph.sosMatrix, ph.ScaleValues, double(pF));
     %S_base = filtfilt(pb.sosMatrix,pb.ScaleValues, double(pF));
     S_base = repmat(mean(pF(1:12,:),1),size(pF,1),1);
     S_base(S_base(:) < 1e3) = 1e3;
     Snorm = single(S_ioi./S_base);
     dat(:, (marks(ind)+1):marks(ind+1) )= Snorm;
end
fprintf('100%%.');
fprintf('\nSaving...\n');
dat = permute(dat, [2 3 1]);
fFlow = fopen([FolderName filesep 'Flow.dat'], 'w');
fwrite(fFlow, dat, 'single');
fclose(fFlow);
fptr = matfile([FolderName filesep 'Flow_infos.mat'], 'Writable', true);
fptr.Stim = Iptr.Stim;
fptr.datLength = Iptr.datLength-1;
fptr.datSize = Iptr.datSize;
fptr.Freq = Iptr.Freq;
fptr.datFile = [FolderName filesep 'Flow.dat'];
fprintf('Done!\n');
end

function speed = private_flow_from_contrast(contrast,T)
% trouve la vitesse � partir des images de contraste
%[nx ny] = size(contrast);

% Correct for points that cannot be...
contrast(isnan(contrast)|contrast<0)=0;

% Not sure about this, verify
contrast2=contrast(3:end-2,3:end-2);
mmean=mean(contrast2(1:end));
sstd=std(contrast2(1:end));

% Build non-linear curve between contrast and correlation time (tau)
tau=(logspace(-11,-2,30).^.5); % Correlation time
K  = ((tau/(2*T)).*(1-exp(-2*T*ones(size(tau))./tau))).^(1/2);
% Find values for which the mean contrast is in the middle
[dummy index1]=find(K>(mmean-3*sstd),1);
[dummy index2]=find(K>(mmean+3*sstd),1);
if isempty(index1), index1=1; end
if isempty(index2)||index2==index1, index2=30; end

% For these values, build a log-linear vector on which contrast is computed
Tau2=(logspace(log10(tau(index1)),log10(tau(index2)),40));
K  = ((Tau2/(2*T)).*(1-exp(-2*T*ones(size(Tau2))./Tau2))).^(1/2);

% Add limit points for interpolation
Tau2=[Tau2(1) Tau2 Tau2(end)];
K= [ 0 K 1e30];
% Interpolate contrast image on these values to obtain correlation time
Tau3=interp1(K,Tau2,contrast); %This is SLOW
% Get speed from correlation time
speed=1./Tau3;
end