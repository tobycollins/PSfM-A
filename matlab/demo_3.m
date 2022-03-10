function demo_3(motionEstimationMode)
%a simple demo that shows how to use PSfM-O to reconstruct surface normals
%using inter-image affine motion. In this demo, the motion is computed with respect to a reference image
%(without loss of generality, this is the first image).

%first we need to create a simulatated scene:
sceneOpts.M = 3; %number of views
sceneOpts.N = 3; %number of points (used to compute affine motion, with a minimum of 3 points)
sceneOpts.maxTilt = 80; % maximum tilt angle of structure in each view (degrees)
sceneOpts.colinearThresh = 0.05; % prevents structure that are very colinear. Implemented with e2/e1 > colinearThresh where (e1, e2) are the PCA eigenvalues of the structure points
sceneOpts.sigma = 0; % no noise

%generate the scene as a set of camera poses (Rgt, Tgt), a structure matrix
%(Sgt) and a set of observed points in each image (qs):
[Rgt,Tgt,Sgt,qs] = generateRandomScene(sceneOpts);


%estimate affine motion:
switch motionEstimationMode
    case 'oneReferenceView'
        %estimates the affine motion from 1st view. This computes a 3x3xN matrix
        %(AFactor) where each AFactor(:,:,i) gives the affine motion from view 1 to view
        %i. Affine motion matrices are computed in homogeneous coordinates.
        AFactor = estimateAffineMotionFrom1stView(qs);

    case 'interViewMotion'
        %estimate the affine motion between all views, then factorizes it
        %with a rank-2 decomposition
        AsCell = estimateAffineMotionBetweenAllViews(qs);

        %compute rank-2 factorization of AsMat:
        AsMat = cell2mat(AsCell);
        [UL,S,VR] = svd(AsMat);

        %extract left rank-2 factor, and reshape it into a 2x2xM matrix (required
        %by estimateNormalsFromAffineMotion)
        cnt = 1;
        AFactor = zeros(2,2,size(AsCell,1));
        for i=1:size(AsCell,1)
            AFactor(:,:,i) = UL(cnt:cnt+2-1,1:2);
            cnt = cnt + 2;
        end
    otherwise
        disp('demo_3 requires one argument (specifying how the affine motion factor is computed. Use either oneReferenceView or interViewMotion')
end


%Now that we have the affine motion from 1st view, we can reconstruct the
%normals. Note that there will be two normal solutions for each view (flip ambiguity).
%In total, if there are u upgrade matrices, there will be u x 2^M normal
%solutions, where M is the number of views.

forceSingleUpgradeMatrix = false; %if forceSingleUpgradeMatrix, then only a single upgrade matrix will be computed, which is the one with the lowest reconstruction error. WARNING!! for some scenes (especially in the minimal case of 3 views, there will be more than 1 valid upgrade matrix, so if forceSingleUpgradeMatrix = true, you may not obtain the correct upgrade matrix.

[normalSolutions, rotationSolutions, alphaSolutions, Xs]  = estimateNormalsFromAffineMotion(AFactor, forceSingleUpgradeMatrix);

%check that solutions are correct using ground truth:
normalErrors = ones(1,sceneOpts.M)*inf;
for i=1:sceneOpts.M
    for j = 1:size(normalSolutions,1)
        n1 = normalSolutions{j,i}(:,1);
        n2 = normalSolutions{j,i}(:,2);

        nGT = Rgt(1:3,end,i);
        df1 = norm(n1-nGT);
        df2 = norm(n2-nGT);
        normalErrors(i) = min(normalErrors(i), min(df1, df2));
    end
end

assert(mean(normalErrors)<1e-10)
disp(['demo_3 completed with normal error of ', num2str(mean(normalErrors))])