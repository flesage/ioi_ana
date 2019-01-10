function HBLookupTableGen()
[e_hbo, e_hbr] = ioi_get_extinctions(450, 700, 50);
e_hbo = single(e_hbo);
e_hbr = single(e_hbr);
pathlength = ioi_path_length_factor(450, 700, 50, 0.100, 'Dunn');
pathlength = single(pathlength);
SP = load('SysSpecs.mat');
Illumination = single(SP.RoyalBlue(1:20:end).*SP.FF02_472_30(1:20:end));
Detection =  single(SP.Camera(1:20:end).*SP.FF01_496_LP(1:20:end));

[HbO, HbR] = meshgrid(single(-50:0.025:50), single(-50:0.025:50));

Iillz = Illumination.*SP.GCaMP6_abs(1:20:end);
Iillz = sum(Iillz(:));
Idetz = Detection.*SP.GCaMP6_em(1:20:end);
Idetz = sum(Idetz(:));
IllumFact = Illumination.*SP.GCaMP6_abs(1:20:end)./Iillz;
DetectFact = Detection.*SP.GCaMP6_em(1:20:end)./Idetz;

extinct = bsxfun(@times, HbO, permute(e_hbr,[3 1 2])) + ...
    bsxfun(@times, HbR, permute(e_hbr,[3 1 2]));
extinct = extinct*1e-6;
extinct = bsxfun(@times, extinct, permute(pathlength,[3 1 2]));
extinct = 10.^(-extinct);

tmp = bsxfun(@times, extinct, permute(IllumFact,[3 1 2]));
CorrIllum = sum(tmp,3);

tmp = bsxfun(@times, extinct, permute(DetectFact,[3 1 2]));
CorrDetect = sum(tmp,3);

Corr = 1./(CorrIllum.*CorrDetect);

Folder = mfilename('fullpath');
Folder = Folder(1:strfind(Folder,'\HBLookup'));
save([Folder 'FactCorrHB.mat'], 'Corr');
end