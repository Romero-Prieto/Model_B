function [y,lambda,B,BB,DD] = BSplineEM(x,n,q,p,lambda,Y)
% Eilers and Marx (1996). Flexible Smoothing with B-splines and Penalties, Statistical Science 11(2):89-121

x      = diag(diag(x));
dx     = (max(x) - min(x))/n;
T      = min(x) + (-q:n - 1)*dx;
P      = (x - T)/dx;
B      = (x >= T) & (x < T + dx);
B(end) = 1;

r      = [2:numel(T) 1];
for k = 1:q
    B = (P.*B + (k + 1 - P).*B(:,r))/k;
end
BB     = B'*B;

D      = diff(eye(size(B,2)),p);
DD     = D'*D;

if isequal(numel(lambda),0)
    lambda = ones(size(Y,2),1)*10^-2;
    sEL    = find(~isnan(lambda));
    I      = NaN(size(lambda));
    r      = 0;
    E      = NaN(numel(lambda),200);
    while numel(sEL) > 0 & r < 200
        [I(sEL) E(sEL,r + 1)] = dGCV(B,BB,DD,lambda(sEL),Y(:,sEL));
        lambda(sEL)           = lambda(sEL).*exp(I(sEL));
        sEL                   = find(abs(I) > 10^-10);
        r                     = r + 1;
    end
    r
end
[y,e]  = GCV(B,BB,DD,lambda,Y);
end

function [y,e] = GCV(B,BB,DD,lambda,Y)
if isequal(lambda,ones(size(lambda))*lambda(1))
    Q = pinv(BB + lambda(1)*DD);
    y = B*Q*B'*Y;
    s = sum((Y - y).^2,1);
    t = sum(diag(Q*BB));
    e = s/(size(B,1) - t)^2;
else
    for j = 1:numel(lambda)
        Q{j}   = pinv(BB + lambda(j)*DD);
        y(:,j) = B*Q{j}*B'*Y(:,j);
        s      = sum((Y(:,j) - y(:,j)).^2);
        t      = sum(diag(Q{j}*BB));
        e(j)   = s/(size(B,1) - t)^2;
    end
end
end

function [I,E] = dGCV(B,BB,DD,lambda,Y)
delta = 0.05;
for j = 1:3
    [~,e{j}] = GCV(B,BB,DD,lambda.*exp(delta*(j - 2)),Y);
end
de    = (e{3} - e{1})'/(2*delta);
d2e   = (e{3} - 2*e{2} + e{1})'/(delta^2);
I     = -de./(max(abs(d2e),10^-5).*sign(d2e));
E     = e{2}';
end