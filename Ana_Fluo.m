function out = Ana_Fluo(FolderName)

fprintf('Opening files.\n');
FileList = dir([FolderName filesep 'Data_Fluo.mat']);
if( isempty(FileList) )
    disp(['No Fluorescence data files found in ' FolderName ' Folder.']);
    disp('Fluorescence Analysis will not run');
    return;
end
Iptr = matfile([FolderName filesep 'Data_Fluo.mat']);

FileList = dir([FolderName filesep 'Data_HBs.mat']);
if( ~isempty(FileList) )
%      fprintf('*** HB correction ***.\n');
%      HBInfo = matfile([FolderName filesep 'Data_Hbs.mat']);
%      nx = HBInfo.datSize(1,1);
%      ny = HBInfo.datSize(1,2);
%      HB_nbFrames = HBInfo.datLength;
%      
%      HBOptr = memmapfile(HBInfo.datFileHbO, 'Format', 'single');
%      HBRptr = memmapfile(HBInfo.datFileHbR, 'Format', 'single');
%      extHbO = (HBOptr.Data(:) + 60)*1e-6;
%      extHbR = (HBRptr.Data(:) + 40)*1e-6;
%      clear HBOptr HBRptr
%      [e_hbo, e_hbr] = ioi_get_extinctions(450, 700, 1000);
%      load('SysSpecs.mat');
%      
%      HbCorr = zeros(nx*ny*HB_nbFrames,1, 'single');
%      marks = round(linspace(1,double(nx*ny*HB_nbFrames),nx));
%      Prctage = linspace(1,nx,20);
%      iP = 1;
%      Illumination = RoyalBlue.*FF02_472_30.*GCaMP6_abs;
%      Detection = GCaMP6_em.*FF01_496_LP.*Camera;
%      for ind = 2:length(marks)
%          disp(ind)
%          if( ind > Prctage(iP) )
%             fprintf('%d%%...',100*round(Prctage(iP)/nx));
%             iP = iP+1;
%          end
%          eHbO = bsxfun(@times, extHbO(marks(ind-1):marks(ind)), e_hbo);
%          eHbR = bsxfun(@times, extHbR(marks(ind-1):marks(ind)), e_hbr);
%          Tmp = bsxfun(@times, Illumination, exp(-0.02*(eHbO + eHbR)));
%          Tmp = bsxfun(@times, sum(Tmp,2), Detection);
%          HbCorr(marks(ind-1):marks(ind)) = sum(Tmp,2);
%      end
%      
%      HBRptr = memmapfile(HBInfo.datFileHbR, 'Format', 'single');
%      HbR = reshape([HBRptr.Data(:)], nx, ny, []);
%      
end
nt = Iptr.datLength;
nx = Iptr.datSize(1,1);
ny = Iptr.datSize(1,2);
Fs = Iptr.Freq;
Dptr = memmapfile(Iptr.datFile, 'Format', 'single');
Fluo = reshape([Dptr.Data(:)], nx, ny, []);

Fluo = imfilter(Fluo, fspecial('gaussian',5,3.3));

f = fdesign.lowpass('N,F3dB', 4, 1/120, Fs);
pb = design(f,'butter');
Fnorm = permute(Fluo,[3 1 2]);
L1 = round(size(Fnorm,1)/5);
L2 = round(4*size(Fnorm,1)/5);
Fnorm = cat(1, Fnorm(L1:-1:1, :, :), Fnorm, Fnorm(end:-1:L2, :, :));
Fnorm = permute(filtfilt(pb.sosMatrix, pb.ScaleValues, double(Fnorm)),[2 3 1]);
Fnorm = Fluo./Fnorm;

end

