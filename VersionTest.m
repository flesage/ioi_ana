function VersionTest(FolderName, Binning)

tOld = dir( [FolderName filesep 'IOI_scan.seq']);
tNew = dir( [FolderName filesep 'img_*.bin']);
if( numel(tOld) > 0 )
    OpenIOI_OldSyst(FolderName, Binning);
elseif( tNew > 0 )
     OpenIOI_NewSyst(FolderName, Binning);
else
      disp('No Valid data set found.');
end

end