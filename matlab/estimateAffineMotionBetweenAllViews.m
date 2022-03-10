function AsMat = estimateAffineMotionBetweenAllViews(qs)
%computes affine motion between all pairs of views. AsMat is a cell array
%where AsMat(j,i) contains the 2x2 affine from view j to view i 
M = length(qs); %number of views
AsMat = cell(M,M);
for i = 1:M
    for j = 1:M
        [A] = ptCorrespond2Affine_2D(qs{i},qs{j});
        AsMat{j,i} = A(1:2,1:2);
    end
end


