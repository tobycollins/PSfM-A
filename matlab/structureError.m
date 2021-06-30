function [errMed,errMean,pointErrs] = structureError(Shat,Strue,withScaleCorrect)
Strue(3,:) = 0;
Shat(3,:) = 0;

dd = distance_vec(Strue,Strue);
mm = max(dd(:));
Strue = 100*Strue/mm;
Shat = 100*Shat/mm;

[errMed,errMean,s,R,T,pointErrs]  = compError(Shat,Strue,withScaleCorrect);


function [errMed,errMean,s,R,T,pointErrs] = compError(S_rec,ps,withScaleCorrect)

if norm(mean(ps,2))>1e-10
    error('ground truth model points should be zero centred');
end


Sbar = mean(S_rec,2);
S_rec(1,:) = S_rec(1,:)-Sbar(1);
S_rec(2,:) = S_rec(2,:)-Sbar(2);


AA = ptCorrespond2Affine_2D(S_rec(1:2,:)',ps(1:2,:)');
if det(AA(1:2,1:2))<0
    S_rec(2,:) = -S_rec(2,:);
end


if nargin<3
    withScaleCorrect =0;
end
[regParams,Bfit,ErrorStats]=absor(S_rec,ps,'doScale',withScaleCorrect);

s=[];
R = [];
T = [];

resids = Bfit(1:2,:)-ps(1:2,:);
df = distance_vec(resids,[0;0]);

errMed = median(df);
errMean = mean(df);

pointErrs = df(:);