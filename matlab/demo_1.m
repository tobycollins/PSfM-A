function demo_1()
%a simple demo that 
%1) creates a simulated scene with co-planar points, camera poses and point
%correspondences
%2) Solves plane-bases SfM with the ortographic camera according to the
%PAMI paper
%3) evaluates structure and pose errors.


%define scene simulation parameters:
sceneOpts.M = 10; %number of views
sceneOpts.N = 7; %number of points
sceneOpts.sigma = 0; % IID gaussian noise added to image points
sceneOpts.maxTilt = 80; % maximum tilt angle of structure in each view (degrees)
sceneOpts.colinearThresh = 0.05; % prevents structure that are very colinear. Implemented with e2/e1 > colinearThresh where (e1, e2) are the PCA eigenvalues of the structure points

%generate the scene as a set of camera poses (Rgt, Tgt), a structure matrix
%(Sgt) and a set of observed points in each image (qs):
[Rgt,Tgt,Sgt,qs] = generateRandomScene(sceneOpts);

%perform 2D affine factorization of point observatioin matrix:
[Qmat,AFactor,SFactor,As] = affineFactorize2D(qs);

%find orthographic camera least-squares upgrade solutions:
[Xs] = oUpgradeLS(As);

%make sure the upgrade solutions are unique to within machine precision
Xs = uniqueUpgrades(Xs);

%now evaluate accuracy of upgrade:
%Find the upgraded structure matrix that best aligns with Sgt
[Shat,structErr] = getBestStructureWithGT(Xs,SFactor,Sgt);

%align upgraded structure matrix with Sgt:
[~,ShatAligned,~]=absor(Shat(1:2,:),Sgt(1:2,:),'doScale',false);
ShatAligned(3,:) = 0;

%find the camera poses
[Rhat,that] = resectCamerasPMAR(ShatAligned,qs);

%evaluate camera pose errors
[RNorms,tNorms] = sceneErrorL2(ShatAligned,Rhat,that,Sgt,Rgt,Tgt);

%if all is good, these assertions should pass
assert(structErr < 1e-6);
assert(norm(RNorms) < 1e-6);
assert(norm(tNorms) < 1e-6);









