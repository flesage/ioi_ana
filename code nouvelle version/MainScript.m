FolderPath = ""; %Folder containing the dataset to be converted
BinningRatio = 4; %Binning value: 1 = no binning (1x1); 2 = 2x2 binning, etc.

ImagesClassification(FolderPath, BinningRatio)
ConvertToTiff(FolderPath);