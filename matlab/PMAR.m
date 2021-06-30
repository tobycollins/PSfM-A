function [Rp,Rm,tp,tm,cost,g,singlesolution] = PMAR(X,Y,camera,d,verb)
vlds = ~isnan(sum(X,1)) & ~isnan(sum(Y,1));
X = X(:,vlds);
Y = Y(:,vlds);



% [Rp,Rm,tp,tm,cost,g,singlesolution] = PMAR(X,Y,camera,d,verb)
%
% Solves L2-optimal correspondence-based planar resection of the metric
% affine cameras.
%
% Inputs:
%
%  - X [(3 x m) or (2 x m) matrix]
%   The coplanar model points. If the number of rows is two, the last
%   coordinate is taken as 0.
%
%  - Y [(2 x m) matrix]
%   The image points.
%
%  - (opt) camera [string or integer] -- default: 'wp'
%   The camera model:
%   * 1, 'or', 'OR': the orthographic camera
%   * 2, 'wp', 'WP': the weak-perspective camera
%   * 3, 'pp', 'PP': the paraperspective camera, which requires d
%   Use '' or [] for default.
%
%  - (opt) d [(2 x 1) vector]
%   The projection direction for the paraperspective camera.
%   If unspecified for the paraperspective camera, the camera is switched
%   to weak-perspective. Ignored for the orthographic and weak-perspective
%   cameras.
%
%  - (opt) verb [integer] -- default: 0
%   Verbosity level:
%   * 0: silent
%   * 1: slightly verbose
%   * 2: very verbose.
%
% Outputs:
%
%  - Rp, Rm [(3 x 3) orthonormal matrices]
%   The rotation matrix for the two solutions.
%
%  - tp, tm [(3 x 1) vectors]
%   The translation vector for the two solutions.
%
%  - cost [scalar]
%   The cost of the estimate, as the sum of squared distances for all input
%   points.
%
%  - g [scalar]
%   The scale for the WP and PP cameras.
%
%  - singlesolution [binary]
%   Flag to indicate that there's a single solution.
%   This activates in the O2 case for the paraperspective camera.
%   When false, it does not imply that the two solutions do not match.

%%%
% setup default and verify parameters
%%%
% the default camera is weak-perspective
if ~exist('camera','var'), camera = 'wp'; end
if isempty(camera), camera = 'wp'; end
% the default verbosity level is silent
if ~exist('verb','var'), verb = 0; end
% standardize the camera variable and verify the parameters
switch camera
    case 0
    case 'or'
    case 'OR'
        camera = 'or';
    case 1
    case 'wp'
    case 'WP'
        camera = 'wp';
    case 2
    case 'pp'
    case 'PP'
        camera = 'pp';
        if ~exist('d','var')
            warning(['(PMAR) paraperspective camera chosen but no ' ...
                'projection direction specified, switching to ' ...
                'weak-perspective']);
            camera = 'wp';
        end
    otherwise
        error(['(PMAR) unknown camera model ' camera]);
end
if verb>1, fprintf('(PMAR) chosen camera: %s\n',camera); end

%%%
% prepare the data
%%%
% number of correspondences
m = size(X,2);
% compute the centroids of the inputs points
x = mean(X,2);
y = mean(Y,2);
% centre the input points
Xp = X-repmat(x,1,m);
Yp = Y-repmat(y,1,m);
% reduce the problem
[U,Sigma,V] = svd(Xp);
detU = det(U);
Z = detU*Yp*V(:,1:2);
Zp = detU*Yp*V(:,3:end);
O = sum(sum(Zp.^2));
s1 = Sigma(1,1);
s2 = Sigma(2,2);
if verb>1, fprintf('(PMAR) base cost: %d\n',O); end

%%%
% solve
%%%
switch camera
    case 'or'
        %%%
        % start with the special case 2xx
        %%%
        if verb>1, fprintf('(PMAR) computing the solution in O2\n'); end
        singlesolution = true;
        W = diag([s1,s2]);
        s = sign(det(Z));
        a = s*s1*Z(1,1)+s2*Z(2,2);
        b = s*s1*Z(2,1)-s2*Z(1,2);
        n = sqrt(a^2+b^2);
        a = a/n;
        b = b/n;
        B = [s*a -b ; s*b a];
        hO = sum(sum((B*W-Z).^2));
        if verb>1, fprintf('(PMAR) additional cost: %d\n',hO); end
        %%%
        % try if better solution in cases 122 and 022
        %%%
        if verb>1, fprintf('(PMAR) computing the solutions in SS2\n'); end
        % form and apply the formulation transformation
        S = [0 -1 ; 1 0];
        e2 = Z(1,2)^2+Z(2,2)^2;
        e = sqrt(e2);
        N = [ Z(:,2)'*S ; Z(:,2)' ];
        f = s1^2;
        Z = [Z(:,2)'*S*Z(:,1) 0 ; Z(:,2)'*Z(:,1) e2]/(e*s1);
        s2 = s2/s1;
        s1 = 1;
        W = diag([1,s2]);
        invW = diag([1,1/s2]);
        % form and solve the sextic
        p = sextic(Z(1,1),Z(2,1),Z(2,2),s2);
        try
            r = roots(p);
        catch
            disp('pmar error: ')
            err.X = X;
            err.Y = Y;
            save err err;
            return;
            
        end
        % go through each root
        for j = 1:length(r)
            beta = real(r(j));
            if verb>1, fprintf('(PMAR)\troot %02d\n\t\treal(beta) = %15d\tsextic''s value = %15d\timag(root) = %15d\n',j,beta,polyval(p,beta),imag(r(j))); end
            % recover B from beta, if possible
            if(abs(1+beta)<eps || abs(s2^2+beta)<eps)
                if verb>1, fprintf('\tB cannot be recovered, as beta is too close to sigma1 or sigma2\n'); end
                continue;
            end
            invG = inv(W^2+beta*eye(2));
            Eb = Z*W*invG^2*W*Z';
            d = trace(Eb)-det(Eb)-1;
            if d<0, d = 0; end
            sqrtd = sqrt(d);
            % recovers the two solutions for q
            % these are the same if d = 0, but this is rare
            % do qp
            qp1 = [ Eb(1,2) - sqrtd ; 1-Eb(1,1) ];
            qp2 = [ 1-Eb(2,2) ; Eb(1,2) + sqrtd ];
            if norm(qp1)<eps
                %added by Toby feb 2021
                continue;
            end
            
            if norm(qp1)>norm(qp2), qp = qp1; else qp = qp2; end
            qp = qp / norm(qp);
            % do qm
            qm1 = [ Eb(1,2) + sqrtd ; 1-Eb(1,1) ];
            qm2 = [ 1-Eb(2,2) ; Eb(1,2) - sqrtd ];
            if norm(qm1)>norm(qm2), qm = qm1; else qm = qm2; end
            qm = qm / norm(qm);
            % form the two solutions for B
            Bp = [qp S*qp]*[qp'*Z*W*invG ; qp'*S'*Z*invW];
            Bm = [qm S*qm]*[qm'*Z*W*invG ; qm'*S'*Z*invW];
            % detect which solutions are in SS22
            
            svp = svd(Bp);
            svm = svd(Bm);
            
            if verb>2
                fprintf('\t\tsvp(1) = %d, svp(2) = %d, svm(1) = %d, svm(2) = %d\n',svp(1),svp(2),svm(1),svm(2));
            end
            Bp = Bp/svp(1);
            Bm = Bm/svm(1);
            
            Op = f*sum(sum((Bp*W-Z).^2));
            if Op<hO
                singlesolution = false;
                hO = Op;
                B = Bp;
            end
            if verb>1, fprintf('\t\tBp is in SS22, with sv1 = %d, sv2 = %d and additional cost %d\n',svp(1),svp(2),Op); end
            
            Om = f*sum(sum((Bm*W-Z).^2));
            if Om<hO
                singlesolution = false;
                hO = Om;
                B = Bm;
            end
            if verb>1, fprintf('\t\tBm is in SS22, with sv1 = %d, sv2 = %d and additional cost %d\n',svm(1),svm(2),Om); end
            
        end
        % finalize
        if singlesolution
            Rp = detU*[ B [0;0] ; [0 0 s] ]*U';
            Rm = Rp;
            tp = y-[eye(2) [0;0]]*Rp*x;
            tm = tp;
        else
            % to undo the formulation transformation
            B = N'*B/e;
            %%% solve the rank-1 factorization
            % form the entries of B*B'
            a = B(1,1)^2+B(1,2)^2;
            b = B(1,1)*B(2,1)+B(1,2)*B(2,2);
            c = B(2,1)^2+B(2,2)^2;
            % compute u by rank-1 decomposition
            u = [ sqrt(1-a) ; sign(-b)*sqrt(1-c) ];
            % assemble the rotation part
            Qb = [ B u ];
            Qp = [ Qb ; cross(Qb(1,:),Qb(2,:)) ];
            Qm = Qp.*[1 1 -1 ; 1 1 -1 ; -1 -1 1];
            Rp = detU*Qp*U';
            Rm = detU*Qm*U';
            % retrieve the translation
            tp = y-Rp(1:2,:)*x;
            tm = y-Rm(1:2,:)*x;
        end
        cost = O+hO;
        g = 1;
        if verb>0, fprintf('(PMAR) final cost: %d\n',cost); end
    case 'wp'
        % no hard choice
        singlesolution = false;
        % the matrix to factor in gamma and u
        B = Z*diag([1/s1,1/s2]);
        %%% rank-1 type-4 problem
        % form the entries of B*B'
        a = B(1,1)^2+B(1,2)^2;
        b = B(1,1)*B(2,1)+B(1,2)*B(2,2);
        c = B(2,1)^2+B(2,2)^2;
        % the leading eigenvalue gives gamma-square directly
        g2 = .5*(a+c+sqrt((a-c)^2+4*b^2));
        g = sqrt(g2);
        % compute u by rank-1 decomposition
        u = [ sqrt(1-a/g2) ; sign(-b/g2)*sqrt(1-c/g2) ];
        % assemble the rotation part
        Qb = [B/g u];
        Qp = [ Qb ; cross(Qb(1,:),Qb(2,:)) ];
        Qm = Qp.*[1 1 -1 ; 1 1 -1 ; -1 -1 1];
        Rp = detU*Qp*U';
        Rm = detU*Qm*U';
        % retrieve the translation
        tp = y-detU*[B g*u]*U'*x;
        tm = y-detU*[B -g*u]*U'*x;
        % assign the cost
        cost = O;
        if verb>0, fprintf('(PMAR) final cost: %d\n',cost); end
    case 'pp'
        % no hard choice
        singlesolution = false;
        % the matrix to factor in gamma and u
        B = Z*diag([1/s1,1/s2]);
        % Cholesky and rotation to Z axis
        a = 1+d(1)^2; b = a+d(2)^2;
        ap = sqrt(a); bp = sqrt(b);
        app = 1/ap; bpp = 1/bp;
        c = app*d(1); dd = bpp*d(2);
        H = [ap 0 ; c*d(2) app*bp ];
        invH = [app 0 ; -c*dd ap*bpp ];
        Rd = [ invH' [-bpp*d(1);-dd] ; c app*dd bpp ];
        Bt = invH*B;
        %%% rank-1 type-4 problem
        % form the entries of B*B'
        a = Bt(1,1)^2+Bt(1,2)^2;
        b = Bt(1,1)*Bt(2,1)+Bt(1,2)*Bt(2,2);
        c = Bt(2,1)^2+Bt(2,2)^2;
        % the leading eigenvalue gives gamma square directly
        g2 = .5*(a+c+sqrt((a-c)^2+4*b^2));
        g = sqrt(g2);
        % compute ut by rank-1 decomposition, then u
        ut = [ sqrt(1-a/g2) ; sign(-b/g2)*sqrt(1-c/g2) ];
        u = H*ut;
        % assemble the rotation part
        Qb = [Bt/g ut];
        Qp = [ Qb ; cross(Qb(1,:),Qb(2,:)) ];
        Qm = Qp.*[1 1 -1 ; 1 1 -1 ; -1 -1 1];
        Rp = detU*Rd*Qp*U';
        Rm = detU*Rd*Qm*U';
        tp = y-detU*[B g*u]*U'*x;
        tm = y-detU*[B -g*u]*U'*x;
        cost = O;
        if verb>0, fprintf('(PMAR) final cost: %d\n',cost); end
end

function c = sextic(Z11,Z21,Z22,s2)
% construction of the sextic
t2 = Z11.^2;
t3 = Z21.^2;
t4 = Z22.^2;
t5 = t2.^2;
t6 = t3.^2;
t7 = t4.^2;
t8 = s2.^2;
t9 = t2.*t3.*2.0;
t10 = t3.*t4.*6.0;
t11 = t8.^2;
t12 = t7.*4.0;
t13 = t4.*t5.*2.0;
t14 = t3.*t4.*2.0;
t15 = t5.*t8.*4.0;
t16 = t6.*t8.*4.0;
t17 = t2.*t3.*t8.*8.0;
t18 = t3.*t4.*t8.*1.8e1;
t19 = t5.*t11.*6.0;
t20 = t7.*t8.*8.0;
t21 = t6.*t11.*6.0;
t22 = t11.^2;
t23 = t2.*t3.*t11.*1.2e1;
t24 = t2.*t7.*t8.*2.0;
t25 = t3.*t4.*t8.*6.0;
t26 = t3.*t4.*t11.*1.8e1;
t27 = t5.*t8.*t11.*4.0;
t28 = t7.*t8.*2.0;
t29 = t7.*t11.*4.0;
t30 = t6.*t8.*t11.*4.0;
t31 = t2.*t3.*t8.*t11.*8.0;
t32 = t3.*t4.*t11.*6.0;
t33 = t3.*t4.*t8.*t11.*6.0;
t34 = t5.*t22;
t35 = t7.*t11;
t36 = t6.*t22;
t37 = t2.*t4.*t7.*t8;
t38 = t2.*t3.*t22.*2.0;
t39 = t3.*t4.*t8.*t11.*2.0;
t40 = t2.*t3.*t7.*t11.*2.0;
c = [t5+t6+t7+t9+t14-t2.*t4.*2.0;t5.*2.0+t6.*2.0+t10+t12+t15+t16+t17+t25+t28+t2.*t3.*4.0-t2.*t4.*6.0-t2.*t4.*t8.*6.0;t5+t6+t7.*6.0+t9+t10+t13+t18+t19+t20+t21+t23+t24+t32+t35-t2.*t4.*6.0-t2.*t5-t2.*t6.*3.0-t3.*t5.*3.0-t2.*t7-t3.*t6-t3.*t7-t4.*t6.*2.0+t5.*t8.*8.0+t6.*t8.*8.0+t2.*t3.*t8.*1.6e1-t2.*t4.*t8.*1.8e1-t2.*t4.*t11.*6.0-t4.*t5.*t8-t3.*t7.*t8.*2.0-t4.*t6.*t8-t4.*t7.*t8-t2.*t3.*t4.*t8.*2.0;t12+t13+t14+t15+t16+t17+t18+t26+t27+t29+t30+t31+t39-t2.*t4.*2.0-t2.*t7.*2.0-t3.*t7.*2.0-t4.*t6.*2.0+t7.*t8.*1.2e1+t5.*t11.*1.2e1+t6.*t11.*1.2e1-t2.*t4.*t8.*1.8e1-t2.*t5.*t8.*4.0+t2.*t3.*t11.*2.4e1-t2.*t6.*t8.*1.2e1-t3.*t5.*t8.*1.2e1-t2.*t4.*t11.*1.8e1+t2.*t7.*t8.*4.0-t3.*t6.*t8.*4.0+t4.*t5.*t8.*4.0-t3.*t7.*t8.*8.0-t4.*t6.*t8.*8.0-t4.*t7.*t8.*4.0+t2.*t7.*t11.*2.0-t4.*t5.*t11.*2.0-t3.*t7.*t11.*2.0-t4.*t6.*t11.*2.0-t2.*t3.*t4.*t8.*4.0-t2.*t3.*t4.*t11.*4.0-t2.*t4.*t8.*t11.*2.0;t7+t19+t20+t21+t23+t24+t25+t26+t33+t34+t36+t37+t38-t2.*t7-t3.*t7+t7.*t11.*6.0+t2.*t3.*t7-t2.*t4.*t8.*6.0-t2.*t4.*t11.*1.8e1+t4.*t5.*t8.*5.0-t2.*t5.*t11.*6.0-t3.*t7.*t8.*1.0e1-t4.*t6.*t8.*7.0-t2.*t6.*t11.*1.8e1-t3.*t5.*t11.*1.8e1-t4.*t7.*t8.*6.0+t2.*t7.*t11.*5.0-t3.*t6.*t11.*6.0+t4.*t5.*t11.*2.0-t5.*t7.*t8.*2.0-t3.*t7.*t11.*7.0-t4.*t6.*t11.*1.0e1+t5.*t8.*t11.*8.0+t6.*t8.*t11.*8.0-t2.*t3.*t4.*t8.*2.0+t2.*t4.*t5.*t8-t2.*t3.*t4.*t11.*8.0+t2.*t4.*t6.*t8+t3.*t4.*t5.*t8.*2.0+t2.*t3.*t7.*t11+t2.*t3.*t8.*t11.*1.6e1-t2.*t4.*t8.*t11.*6.0-t4.*t5.*t8.*t11-t4.*t6.*t8.*t11-t2.*t3.*t4.*t8.*t11.*2.0;t27+t28+t29+t30+t31+t32+t33+t40+t5.*t22.*2.0+t6.*t22.*2.0-t2.*t4.*t11.*6.0-t3.*t7.*t8.*4.0-t4.*t7.*t8.*4.0+t2.*t7.*t11.*4.0+t4.*t5.*t11.*4.0-t5.*t7.*t8.*2.0-t3.*t7.*t11.*8.0-t4.*t6.*t11.*8.0-t5.*t7.*t11.*2.0+t2.*t3.*t22.*4.0-t2.*t3.*t4.*t11.*4.0+t2.*t3.*t7.*t8.*2.0+t2.*t4.*t7.*t8.*2.0+t2.*t4.*t5.*t11.*2.0+t2.*t4.*t6.*t11.*2.0+t3.*t4.*t5.*t11.*4.0-t2.*t4.*t8.*t11.*6.0-t2.*t5.*t8.*t11.*4.0-t2.*t6.*t8.*t11.*1.2e1-t3.*t5.*t8.*t11.*1.2e1-t3.*t6.*t8.*t11.*4.0-t4.*t6.*t8.*t11.*4.0-t2.*t3.*t4.*t8.*t11.*4.0;t34+t35+t36+t37+t38+t39+t40-t4.*t7.*t8+t2.*t7.*t11-t3.*t7.*t11.*3.0-t5.*t7.*t11.*2.0-t2.*t5.*t22-t2.*t6.*t22.*3.0-t3.*t5.*t22.*3.0-t3.*t6.*t22-t2.*t4.*t8.*t11.*2.0+t4.*t5.*t8.*t11-t4.*t6.*t8.*t11.*3.0-t2.*t3.*t4.*t8.*t11.*2.0+t2.*t4.*t5.*t8.*t11+t2.*t4.*t6.*t8.*t11+t3.*t4.*t5.*t8.*t11.*2.0];