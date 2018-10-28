% Export to R importer
load('/main/calen/occluding/images/ImgStats/ImgStatsDS1.mat')
ImgStatsExport = ImgStatsDS1;

filePath = ['~/Dropbox/Calen/Dropbox/'];

L = ImgStatsExport.L(:);
L = single(L);
save([filePath, '/L.mat'], 'L', '-v7.3');

C = ImgStatsExport.C(:);
C = single(C);
save([filePath, '/C.mat'], 'C', '-v7.3');

patchX = repmat(ImgStatsExport.sampleCoords(:,1), 1491,1);
save([filePath, '/patchX.mat'], 'patchX', '-v7.3');

patchY = repmat(ImgStatsExport.sampleCoords(:,2), 1491,1);
save([filePath, '/patchY.mat'], 'patchY', '-v7.3');

imgNr = repelem(1:1491, 90297)';
save([filePath, '/imgNr.mat'], 'patchY', '-v7.3');

%imgFilePath = repmat(ImgStatsExport.imgDir(:), 90297, 1);
%save([filePath, '/imgFilePath.mat'], 'imgFilePath');

%bNatural = repmat(ImgStatsExport.bNatural(:), 90297, 1);
%save([filePath, '/bNatural.mat'], 'bNatural');

ImgStatsExportStruct.SaVertical   = ImgStatsExport.Sa(:,:,1); ImgStatsExportStruct.SaVertical   = ImgStatsExportStruct.SaVertical(:);
ImgStatsExportStruct.SaHorizontal = ImgStatsExport.Sa(:,:,2); ImgStatsExportStruct.SaHorizontal = ImgStatsExportStruct.SaHorizontal(:);
ImgStatsExportStruct.SaBowtie     = ImgStatsExport.Sa(:,:,3); ImgStatsExportStruct.SaBowtie     = ImgStatsExportStruct.SaBowtie(:);
ImgStatsExportStruct.SaSpot       = ImgStatsExport.Sa(:,:,4); ImgStatsExportStruct.SaSpot       = ImgStatsExportStruct.SaSpot(:);

verticalBin = zeros(3, 1491 * 90297);

vertical = ImgStatsExport{1};
for i = 1:10
    for j = 1:10
        for k = 1:10
            verticalIdx = vertical{i,j,k};
            verticalBin = [verticalBin;repmat([i,j,k], length(verticalIdx),1)];
        end
    end
end


save('~/Dropbox/Calen/Work/ImgStatsExportStruct.mat', 'ImgStatsExportStruct')
