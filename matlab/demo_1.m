function demo_1()

sceneOpts.M = 10; %number of views
sceneOpts.N = 7; %number of points
sceneOpts.sigma = 0; % IID gaussian noise added to image points
sceneOpts.maxTilt = 80; % maximum tilt angle of structure in each view (degrees)
sceneOpts.colinearThresh = 0.05; % prevents structure that are very colinear. Implemented with e2/e1 > colinearThresh where (e1, e2) are the PCA eigenvalues of the structure points

[Rgt,Tgt,Sgt,qs] = generateRandomScene(sceneOpts);

[Qmat,AFactor,SFactor,As] = affineFactorize2D(qs);
[Xs] = closedFormUpgrade(As);
Xs = uniqueUpgrades(Xs);

[Shat,structErr] = getBestStructureWithGT(Xs,SFactor,Sgt);
[Rhat,that] = resectCamerasPMAR(Shat,qs);

[RNorms,tNorms] = sceneErrorL2(Shat,Rhat,that,Sgt,Rgt,Tgt);
''







