function [R22, R1, R2, alpha] = ss22Decomposition(A)
%decomposes a 2x2 matrix A int A = alpha*R22, where R22 is a 2x2
%sub-stiefel matric. The two matrics R1 and R2 are the full matrices for
%which R22 is the top-left 2x2 submatrix
[u,s,v] = svd(A);
alpha = s(1,1);
R22 = A./alpha;
a = sqrt(1 - norm(R22(1,:))^2);
b = sqrt(1 - norm(R22(2,:))^2);

R1 = [R22, [a;b]];
d1 = R1(1,:)*R1(2,:)';
d2 = R1(1,:)*[R1(2,1), R1(2,2), -R1(2,3)]';

if abs(d1) > abs(d2)
    R1 = [R22, [a;-b]];
end
R2 = [R1(:,1:2),-R1(:,end)];

R1 = [R1; cross(R1(1,:), R1(2,:))];
R2 = [R2; cross(R2(1,:), R2(2,:))];


