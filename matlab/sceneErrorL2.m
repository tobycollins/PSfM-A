function [RNorms,tNorms] = sceneErrorL2(Shat,Rhat,that,Sgt,Rgt,Tgt)

M = length(Rgt);
RNorms = zeros(1,M);
tNorms = zeros(1,M);


for i=1:M
    RsEst = Rhat{i};
    tsEst = that{i};
    rn = zeros(2,1);
    for j=1:2
       R = RsEst(:,:,j) -  Rgt(:,:,i);
       rn(j) = norm(R);
    end
    [a,b] = min(rn);
    tNorms(i) = norm(tsEst(:,b) - Tgt(:,i));
    RNorms(i) = a;
    
end

