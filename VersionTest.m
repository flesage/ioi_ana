function VersionTest(FolderName, Binning)


tOld = dir( [FolderName filesep 'IOI_scan.seq']);
tNew = dir( [FolderName filesep 'img_*.bin']);
if( numel(tOld) > 0 )
    disp('System version detected: IOI 1.0');
    OpenIOI_OldSyst(FolderName, Binning);
elseif( tNew > 0 )
    header = memmapfile([FolderName filesep 'img_00000.bin'], ...
    'Offset', 0, 'Format', 'int32', 'repeat', 1);
    %FLAG for header version here!
    headerVer = header.Data;
    if( headerVer == 2 )
        disp('System version detected: IOI 2.1');
        OpenIOI_NewSyst(FolderName, Binning, 2);
    elseif( headerVer == 1 )
        disp('System version detected: IOI 2.0');
        OpenIOI_NewSyst(FolderName, Binning, 1);
    end
else
      disp('No Valid data set found.');
end

end