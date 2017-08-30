function V = VersionTest(FolderName)


tOld = dir( [FolderName filesep 'IOI_scan.seq']);
tNew = dir( [FolderName filesep 'img_*.bin']);
if( numel(tOld) > 0 )
    disp('System version detected: IOI 1.0');
    V = '1.0';
    %OpenIOI_OldSyst(FolderName, Binning);
elseif( numel(tNew) > 0 )
    header = memmapfile([FolderName filesep 'img_00000.bin'], ...
    'Offset', 0, 'Format', 'int32', 'repeat', 1);
    %FLAG for header version here!
    headerVer = header.Data;
    if( headerVer == 3 )
        disp('System version detected: IOI 2.2');
        V = '2.2';
    elseif( headerVer == 2 )
        disp('System version detected: IOI 2.1');
        V = '2.1';
    elseif( headerVer == 1 )
        V = '2.0';
        disp('System version detected: IOI 2.0');
    end
else
    V = '0.0';
      disp('No Valid data set found.');
end

end