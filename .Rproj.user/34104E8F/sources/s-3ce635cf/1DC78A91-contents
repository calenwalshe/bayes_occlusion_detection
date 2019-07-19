nX = numel(unique(ImgStatsDS1.sampleCoords(:,2)));
nY = numel(unique(ImgStatsDS1.sampleCoords(:,1)));

minIdx = 1;
maxIdx = size(ImgStatsDS1.L, 1);
imsquare = reshape(ImgStatsDS1.L(:,1), [nX, nY]);
startIdx = (2 * nY + 2);
endIdx = maxIdx - (2*nY + 2);

coord_surround = [0 -474  -476    -2   472   474   476     2  -472];

L = ImgStatsDS1.L;
C = ImgStatsDS1.C;
S = ImgStatsDS1.Sa;

for targ = 1:1
    S = S(:,:,targ);      
    for lBin = 1:10
        for cBin = 1:10            
            for sBin = 1:10
                display(lBin)
                display(cBin)            
                display(sBin)                
                patchIdx    = ImgStatsDS1.patchIndex{targ}{lBin, cBin, sBin};
                [I,J]       = ind2sub(size(L),patchIdx);
                validIdx    = find(I > startIdx & I < endIdx);
                
                patchIdx = patchIdx(validIdx);
                L_statBin = zeros(size(patchIdx, 1), 10);
                C_statBin = zeros(size(patchIdx, 1), 10);
                S_statBin = zeros(size(patchIdx, 1), 10);
             
                for i = 1:size(validIdx,1)
                    centIdx        = patchIdx(i);
                    surroundIdx    = patchIdx(i) + coord_surround;
                    L_statBin(i,:) = [L(surroundIdx), centIdx];
                    C_statBin(i,:) = [C(surroundIdx), centIdx];                                        
                                
                    S_statBin(i,:) = [S(surroundIdx), centIdx];
                end
                
                dlmwrite(['/main/calen/occluding/surround_prediction/target_', num2str(targ), '-stat_', 'L_', 'Bin_L',...
                    num2str(lBin), '_C', num2str(cBin), '_S', num2str(sBin), '.txt'], L_statBin);                
                dlmwrite(['/main/calen/occluding/surround_prediction/target_', num2str(targ), '-stat_', 'C_', 'Bin_L',...
                    num2str(lBin), '_C', num2str(cBin), '_S', num2str(sBin), '.txt'], C_statBin);
                dlmwrite(['/main/calen/occluding/surround_prediction/target_', num2str(targ), '-stat_', 'S_', 'Bin_L',...
                    num2str(lBin), '_C', num2str(cBin), '_S', num2str(sBin), '.txt'], S_statBin);                
            end
        end
    end
end
