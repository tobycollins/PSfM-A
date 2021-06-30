function Xs_ = uniqueUpgrades(Xs)
Xs_= [];
for i=1:size(Xs,3)
   X = Xs(:,:,i);
   if isIn(Xs_,X)
      
   else
       Xs_ = cat(3,Xs_,X);
   end
end


function t = isIn(A,B)
t = false;
if isempty(A)
    t = false;
    
else
   for i=1:size(A,3)
      A_ = A(:,:,i);
      if norm(A_(:)-B(:)) < eps
         t = true;
         return;
      end
   end
end
    
