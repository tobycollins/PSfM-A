function [X, cr] = findBestFittingUpgradeMatrix(Xs, As)
%finds the best fitting upgrade matrix, given a set of upgrade matrices
%(Xs) and the corresponding affine motion matrix As.
%returns best upgrade matrix (X) and the cost ratio between best and second
%best matrix. 

M = size(As,3); %number of views
U = size(Xs,3); %number of upgrade matrices

costs = zeros(U,1);
for u = 1:U
    X = Xs(:,:,u);
    alphas =  zeros(M,1);
    for i = 1:M
        B = As(1:2,1:2,i)*X;
        [R22, R1, R2, alpha] = ss22Decomposition(B);
        alphas(i) = alpha;
    end
    meanAlpha = mean(alphas);

    for i = 1:M
        B = As(1:2,1:2,i)*X / meanAlpha;
        [R22, R1, R2, alpha] = ss22Decomposition(B);
        costs(u) = costs(u) + norm(R22 - B);
    end
end

[m, indx] = min(costs);
X = Xs(:,:,indx);
csts = sort(costs)/M;
if length(csts)>1
    cr = (csts(2) / csts(1));
else
    cr = 0.0;
end
