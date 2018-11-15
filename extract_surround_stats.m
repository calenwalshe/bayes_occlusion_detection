nX = numel(unique(ImgStatsDS1.sampleCoords(:,1)));
nY = numel(unique(ImgStatsDS1.sampleCoords(:,2)));

minIdx = 1;
maxIdx = size(ImgStatsDS1.L,1);

startIdx = nX + 1; % boundaries do not have surround
endIdx   = maxIdx - nX;
coord_surround = [0, -nX, -nX - 1, -1, (nX - 1), nX, (nX + 1), 1, (1 - nX)];

imL = ImgStatsDS1.L(:,1);

L = ImgStatsDS1.L;

statSurround = cell(10,10,10,4);
for targ = 1:4
    for lBin = 1:10
        for cBin = 1:10
            for sBin = 1:10
                
                %[imX, imY]  = ind2sub([nX, nY], I);

                %imX = setdiff(imX, minmax(ImgStatsDS1.sampleCoords(:,1)));
                %imY = setdiff(imY, minmax(ImgStatsDS1.sampleCoords(:,2)));
                
                patchIdx    = ImgStatsDS1.patchIndex{targ}{lBin, cBin, sBin};
                [I,J]       = ind2sub(size(L),patchIdx);            
                
                validIdx    = find(I > startIdx & I < endIdx);
                
                patchIdx    = patchIdx(validIdx);
                
                statBin = zeros(size(patchIdx,1), size(coord_surround,2) + 5);                
                for i = 1:size(patchIdx,1)                
                    centIdx     = patchIdx(i);
                    surroundIdx = patchIdx(i) + coord_surround;

                    statBin(i,:) = [L(surroundIdx), centIdx, lBin, cBin, sBin, targ];                
                end
                
                statSurround{lBin, cBin, sBin, targ} = statBin;
            end
        end
    end
end