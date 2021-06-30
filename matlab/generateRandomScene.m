function [Rs,Ts,S,qs] = generateRandomScene(sceneOpts)
%generate structure
trm = false;
while ~trm
   S = rand(2,sceneOpts.N);
   S = S - repmat(mean(S,2),1, sceneOpts.N);
   [u,s,v] = svd(S);
   e1 = sqrt(s(1,1));
   e2 = sqrt(s(2,2));
   
   if e2/e1 > sceneOpts.colinearThresh
       
       trm = true;
   end
end
S(3,:) = 0;

Rs = zeros(3,3,sceneOpts.M);
Ts = zeros(2,sceneOpts.M);

for i=1:sceneOpts.M
    trm = false;
    while ~trm
       [R,~,~] = svd(rand(3,3));
       R(:,3) = cross(R(:,1),R(:,2));
       ang = 360*acos(R(3,3))/(2*pi);
       if abs(ang)<sceneOpts.maxTilt
           trm = true;
       end
    end
    Rs(:,:,i) = R;
    Ts(:,i) = (rand(2,1) - 0.5); 
end
for i=1:sceneOpts.M
    Q = Rs(:,:,i) * S;
    q = Q(1:2,:) + repmat(Ts(:,i),1,sceneOpts.N);
    ns = randn(2,sceneOpts.N)*sceneOpts.sigma;
    qs{i} = q + ns;
end