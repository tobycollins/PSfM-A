function [Qmat,AFactor,SFactor,As] = affineFactorize2D(qs_n)
%factorisation:
numViews = length(qs_n);
for i=1:numViews
    qs_ = qs_n{i};
    qBar = mean(qs_,2);
    qsMeaned{i} = qs_;
    qsMeaned{i}(1,:) = qsMeaned{i}(1,:) - qBar(1);
    qsMeaned{i}(2,:) = qsMeaned{i}(2,:) - qBar(2);
end
M = qsMeaned';
M = cell2mat(M);
Qmat = M;
[UU,SS,VV] = svd(M);
S22 = sqrt(SS(1:2,1:2));

AFactor = UU(:,1:2)*S22;
SFactor = S22*VV(:,1:2)';


mm1 = mean(SFactor(1,:));
mm2 = mean(SFactor(2,:));
S_ = SFactor;
S_(1,:) = SFactor(1,:)-mm1;
S_(2,:) = SFactor(2,:)-mm2;


[EE,VV] = eigs(S_*S_');
G1 = inv(EE);
G2 = [1/sqrt(VV(1,1)),0;0,1/sqrt(VV(2,2))];
GG = G2*G1;
if det(GG)<0
    GG = [-1,0;0,1]*GG;
end

GG = GG*100;
AFactor = AFactor*inv(GG);
SFactor = GG*SFactor;

cc=1;
for i=1:numViews
    As(:,:,i) = AFactor(cc:cc+1,:);
    cc=cc+2;
end

for i=1:numViews
    dts(i) = sign(det(As(:,:,i)));
end
if median(dts)<0
    SFactor = [1,0;0,-1]*SFactor;
    AFactor = AFactor*[1,0;0,-1];
end


cc=1;
for i=1:numViews
    As(:,:,i) = AFactor(cc:cc+1,:);
    cc=cc+2;
end
