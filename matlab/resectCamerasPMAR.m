function [Rhat,that] = resectCamerasPMAR(Shat,qs)
for i=1:length(qs)
    [Rp,Rm,tp,tm,cost] = PMAR(Shat,qs{i},'OR',[]);
    Rhat{i} = cat(3,Rp,Rm);
    that{i} = cat(2,tp,tm);
end