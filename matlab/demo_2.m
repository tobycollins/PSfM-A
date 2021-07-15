function demo_2()
%a simple demo that shows reconstruction accuracy w.r.t. image noise. 
%For each image noise level, a random 3D secne is created, which is then
%reconstructed with plane-bases SfM with the ortographic camera according to the
%PAMI paper


%define scene simulation parameters:
sceneOpts.M = 10; %number of views
sceneOpts.N = 7; %number of points
sceneOpts.maxTilt = 80; % maximum tilt angle of structure in each view (degrees)
sceneOpts.colinearThresh = 0.05; % prevents structure that are very colinear. Implemented with e2/e1 > colinearThresh where (e1, e2) are the PCA eigenvalues of the structure points

sigmas = linspace(0,0.1,10);
numberOfTrials = 50; % for each sigma, multiple trials are run, with each trial corresponding to a different random scene
SErrAll = [];
RErrAll = [];
tErrAll = [];
for s=sigmas
    sceneOpts.sigma = s;
    
    
    Serr = [];
    Rerr = [];
    tErr = [];
    
    for j=1:numberOfTrials
        
        %generate the scene as a set of camera poses (Rgt, Tgt), a structure matrix
        %(Sgt) and a set of observed points in each image (qs):
        [Rgt,Tgt,Sgt,qs] = generateRandomScene(sceneOpts);
        
        %perform 2D affine factorization of point observation matrix:
        [Qmat,AFactor,SFactor,As] = affineFactorize2D(qs);
        
        %find orthographic camera least-squares upgrade solutions:
        [Xs] = oUpgradeLS(As);
        
        %make sure the upgrade solutions are unique to within machine precision
        Xs = uniqueUpgrades(Xs);
        
        %now evaluate accuracy of upgrade:
        %Find the upgraded structure matrix that best aligns with Sgt
        [Shat,structErr] = getBestStructureWithGT(Xs,SFactor,Sgt);
        
        %align upgraded structure matreix with Sgt:
        [~,ShatAligned,~]=absor(Shat(1:2,:),Sgt(1:2,:),'doScale',false);
        ShatAligned(3,:) = 0;
        
        %find the camera poses
        [Rhat,that] = resectCamerasPMAR(ShatAligned,qs);
        
        %evaluate camera pose errors
        [RNorms,tNorms] = sceneErrorL2(ShatAligned,Rhat,that,Sgt,Rgt,Tgt);
        
        
        Serr = [Serr,structErr];
        Rerr = [Rerr,RNorms];
        tErr = [tErr,tNorms];
        
    end
    SErrAll = [SErrAll;Serr];
    RErrAll = [RErrAll;Rerr];
    tErrAll = [tErrAll;tErr];
end

figure(1);
clf;
subplot(2,2,1);
errorbar(sigmas,mean(SErrAll,2)',std(SErrAll,[],2)');
xlabel('Image noise (standard deviation)');
ylabel('Structure error');

subplot(2,2,2);
errorbar(sigmas,mean(RErrAll,2)',std(RErrAll,[],2)');
xlabel('Image noise (standard deviation)');
ylabel('Rotation error');

subplot(2,2,3);
errorbar(sigmas,mean(tErrAll,2)',std(tErrAll,[],2)');
xlabel('Image noise (standard deviation)');
ylabel('Translation error');






