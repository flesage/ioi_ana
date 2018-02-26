function out = Ana_Fluo(FolderName, OStream)
%%%%%%%%%%%%%%%%%%%%%%
% Validation & opening
%%%%%%%%%%%%%%%%%%%%%%
if( isempty(OStream) )
    fprintf('Opening files.\n');
else
    OStream.String = sprintf('%s\r%s',...
        'Opening files',...
        OStream.String);
    drawnow;
end

FileList = dir([FolderName filesep 'Data_Fluo.mat']);
if( isempty(FileList) )
    disp(['No Fluorescence data files found in ' FolderName ' Folder.']);
    disp('Fluorescence Analysis will not run');
    return;
end
Iptr = matfile([FolderName filesep 'Data_Fluo.mat'],'Writable',true);

FileList = dir([FolderName filesep 'Data_HBs.mat']);
HBInfo = matfile([FolderName filesep 'Data_Hbs.mat']);
nf = HBInfo.datLength;
%%%%%%%%%%%%%%%%%%%%%%%
% Dimension permutation
%%%%%%%%%%%%%%%%%%%%%%%
Dptr = memmapfile(Iptr.datFile, 'Format', 'single', 'Writable', true);
nx = Iptr.datSize(1,1);
ny = Iptr.datSize(1,2);
nfF = Iptr.datLength;
Fluo = reshape([Dptr.Data(:)], nx, ny, []);
if( nf < nfF )
    Fluo = Fluo(:,:,1:nf);
    Iptr.datLength = nf;
end
Fluo = permute(Fluo,[3 1 2]);
clear Dptr;
fid = fopen(Iptr.datFile, 'w');
fwrite(fid, Fluo,'single');
clear Fluo Dptr fid;

%%%%%%%%%%%%%%%
% Memory Test (limited to 5GB)
%%%%%%%%%%%%%%%
SpP = nfF*4*256;
NbP = floor(5e8/SpP);
if( NbP > double(nx*ny) )
    NbP = double(nx*ny);
end

Ffluo = round(0.5*(Iptr.Freq));
Fbase = round(60*(Iptr.Freq));

%%%%%%%%%%%%%%%
% Hb correction
%%%%%%%%%%%%%%%
if( ~isempty(FileList) )
    
    if( isempty(OStream) )
        fprintf('*** HB correction ***.\n');
    else
        OStream.String = sprintf('%s\r%s',...
            '*** HB correction ***',...
            OStream.String);
        drawnow;
    end
       
    % Constants
    %%%%%%%%%%%
    [e_hbo, e_hbr] = ioi_get_extinctions(450, 700, 1000);
    pathlength = ioi_path_length_factor(450, 700, 1000, 0.100, 'Dunn');
    SP = load('SysSpecs.mat');
    Illumination = SP.RoyalBlue.*SP.FF02_472_30;
    Detection =  SP.Camera.*SP.FF01_496_LP;
 
    % For loop
    %%%%%%%%%%
    Tags = round(linspace(double(0), double(nx*ny), double(20)));
    indT = 1;
    StaticStr = OStream.String;
    for indP = 1:NbP:nx*ny
        
        if( (nx*ny - indP) < NbP )
            NbP = (nx*ny - indP);
        end
        HbOptr = memmapfile(HBInfo.datFileHbO,'Offset', 4*(indP - 1)*nf,...
            'Format', 'single','Repeat', NbP*nf);
        HbRptr = memmapfile(HBInfo.datFileHbR,'Offset', 4*(indP - 1)*nf,...
            'Format', 'single','Repeat', NbP*nf);
        Fptr = memmapfile(Iptr.datFile,'Offset', 4*(indP - 1)*nf,...
            'Format', 'single', 'Writable', true, 'Repeat', NbP*nf);
        
        extHbO = bsxfun(@times, e_hbo, (HbOptr.Data(:) + 60)*1e-6);
        extHbR = bsxfun(@times, e_hbr, (HbRptr.Data(:) + 40)*1e-6);
        
        Att = exp(-bsxfun(@times, pathlength,(extHbO + extHbR)));
        aI = bsxfun(@times, Illumination, Att);
        
        Ci = sum(bsxfun(@times, Illumination, SP.GCaMP6_abs))./...
            sum(bsxfun(@times, aI, SP.GCaMP6_abs),2);
        aI = bsxfun(@times, Detection, Att);
        Cd = sum(bsxfun(@times, Detection, SP.GCaMP6_em))./...
            sum(bsxfun(@times, aI, SP.GCaMP6_em),2);
        
        Corr = Ci.*Cd;
        data = Fptr.Data(:).*Corr;
        data = reshape(data, nf, []);
        Sfluo = medfilt1(data,Ffluo,[],1,'truncate');
        Sbase = medfilt1(data,Fbase,[],1,'truncate');
        Fnorm = Sfluo./Sbase; 
        Fptr.Data(:) = Fnorm(:);
        
        if( indP >= Tags(indT) )
            P = round((100*Tags(indT))/nx*ny);
            if( isempty(OStream) )
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            else
                OStream.String = sprintf('%s\r%s',...
                    ['Completion: ' int2str(P) '%'],...
                    StaticStr);
                drawnow;
            end
            indT = indT + 1;
        end
    end
end
OStream.String = StaticStr;
OStream.String = sprintf('%s\r%s',...
    'Done.',...
    OStream.String);
drawnow;
clear Corr Ci Cd aI Att Detection Illumination e_* HbR* HbO* ind nfF pathlength Fptr

end

