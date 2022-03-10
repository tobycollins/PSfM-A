function [normalSolutions, rotationSolutions, alphaSolutions, Xs] = estimateNormalsFromAffineMotion(As, forceSingleUpgradeMatrix)
%estimates the surface normals from the affine motion factor (As). As can
%either be a 3x3xM metrix (where M is the number of views) or a 2x2xM
%matrix. In the first case, each As(:,:,i) is supposed to be a homogeneous matrix.
%
%if forceSingleUpgradeMatrix, then only a single upgrade matrix will be computed, which is the 
% one with the lowest reconstruction error. WARNING!! for some scenes (especially in the minimal case of 3 views, 
% there will be more than 1 valid upgrade matrix, so if forceSingleUpgradeMatrix = true, 
% you may not obtain the correct upgrade matrix. It is therefore safer to
% use forceSingleUpgradeMatrix = false

%Function return normalSolutions: a cell array of size UxM, where U is the
%number of found upgrade matrices. normalSolutions{u,i} gives the two
%surface normal solutions for view M, using the uth upgrade matrix
%
%rotationSolutions: a cell array of size UxM, where U is the
%number of found upgrade matrices. rotationSolutions{u,i} gives the two
%rotation solutions for view M, using the uth upgrade matrix. Rotation
%solutions are the rotations that transform the surface plane, defined in
%world coordinates on the plane z=0, to camera coordinates.
%
%alphaSolutions: a matrix of size UxM, where alphaSolutions(u,i) gives the
%alpha factor (scaling factor) for view i using the uth upgrade matrix. 
%
%Xs: The found upgrade matrices of size 2x2xU

%step 1: compute upgrade matrices:
[Xs] = oUpgradeLS(As(1:2,1:2,:));
Xs = uniqueUpgrades(Xs);
if forceSingleUpgradeMatrix
    [Xs,cr] = findBestFittingUpgradeMatrix(Xs, As);
end
%step 2: compute scen plane - to - image affine motions, using upgrade
%matrix
numberOfUpgrades = size(Xs,3);
rotationSolutions = cell(numberOfUpgrades, size(As,3));
normalSolutions = cell(numberOfUpgrades, size(As,3));
alphaSolutions = zeros(numberOfUpgrades, size(As,3));
for u = 1:numberOfUpgrades
    X = Xs(:,:,u);
    for i = 1:size(As,3)
        A = As(1:2,1:2,i);
        AUpgraded = A*X;
        [R22, R1, R2, alpha] = ss22Decomposition(AUpgraded);
        rotationSolutions{u,i} = cat(3,R1, R2);
        normalSolutions{u,i} = [R1(:,3), R2(:,3)];
        alphaSolutions(u,i) = alpha;

    end
end