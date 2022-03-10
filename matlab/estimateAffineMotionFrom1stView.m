function As1M = estimateAffineMotionFrom1stView(qs)
%function that computes the affine motion matrices from the first views to
%all other views.
M = length(qs);
for i = 1:M
    [A] = ptCorrespond2Affine_2D(qs{1},qs{i});
    As1M(:,:,i) = A;
end


