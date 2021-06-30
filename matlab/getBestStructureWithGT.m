function [bestS,bestShape] = getBestStructureWithGT(Xs,SFactor,S_gt)
bestS = [];
bestSErr = inf;
if isempty(Xs)
    bestS = [];
    bestShape = [];
    return;
end
for u = 1:size(Xs,3)
    X = Xs(:,:,u);
    Shat = inv(X)*SFactor;
    Shat(3,:) =0;
    [medSErr,meanSErr,pointErrs] = structureError(Shat,S_gt,false);
    
    if meanSErr < bestSErr
        bestShape = meanSErr; 
        bestS = Shat;
    end
end
