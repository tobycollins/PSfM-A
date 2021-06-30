function S = upgradeStructure(Xs,S_aff)
S = [];
cc =0;
if isempty(Xs)
    return;
else
    for i=1:size(Xs,3)
        X =  Xs(:,:,i);
        
        if sum(isnan(X(:)))>0|rank(X)<2
            
        else
            cc = cc+1;
            
            S(:,:,cc) = [inv(X)*S_aff;zeros(1,size(S_aff,2))];
            
        end
        
        
    end
    
end
