function [A,MM,M_mat] = ptCorrespond2Affine_2D(ps1,ps2)
numPts = size(ps1,1);

if numPts == 3;
     %exact solution
    px1=ps1(1,1);
    py1=ps1(1,2);
    px2=ps1(2,1);
    py2=ps1(2,2);
    px3=ps1(3,1);
    py3=ps1(3,2);

   
    M_mat=[[          (py2-py3)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),        -(-py3+py1)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),          (py1-py2)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),                                                                     0,                                                                     0,                                                                     0]
        [                                                                     0,                                                                     0,                                                                     0,          (py2-py3)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),        -(-py3+py1)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),          (py1-py2)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3)]
        [        -(-px3+px2)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),          (px1-px3)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),         -(px1-px2)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),                                                                     0,                                                                     0,                                                                     0]
        [                                                                     0,                                                                     0,                                                                     0,        -(-px3+px2)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),          (px1-px3)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),         -(px1-px2)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3)]
        [  (px2*py3-px3*py2)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3), -(px1*py3-px3*py1)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),  (px1*py2-px2*py1)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),                                                                     0,                                                                     0,                                                                     0]
        [                                                                     0,                                                                     0,                                                                     0,  (px2*py3-px3*py2)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3), -(px1*py3-px3*py1)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3),  (px1*py2-px2*py1)/(-px2*py1+px1*py2+px3*py1-px3*py2-px1*py3+px2*py3)]];


    aParams = M_mat*ps2(:);
    A=[aParams(1:2),aParams(3:4),aParams(5:6)];
    A(3,3)=1;
    MM=[];

else

    %M = [];
    %B = [];
    % for i=1:numPts
    %    x_p = ps1(i,1);
    %    y_p = ps1(i,2);
    %    M = [M; [ x_p,   0, y_p,   0,   1,   0]];
    %    M = [M; [   0, x_p,   0, y_p,   0,   1]];
    %    B = [B;ps2(i,1);ps2(i,2)];
    % end
    z = zeros(numPts,1);
    M1 = [ps1(:,1),z,ps1(:,2),z,[z+1],z];
    M2 = [z,ps1(:,1),z,ps1(:,2),z,[z+1]];
    MM = [M1;M2];
    BB = [ps2(:,1);ps2(:,2)];
    %As = MM\BB;

    %As = inv(MM'*MM)*MM'*BB;
    M_mat = inv(MM'*MM)*MM';
    As=M_mat*BB;
    %As = M\B;
    A = [As(1:2),As(3:4),As(5:6)];
    A = [A;[0,0,1]];

end




