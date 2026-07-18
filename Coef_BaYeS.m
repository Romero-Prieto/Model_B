function beta = Coef_BaYeS(q,x,ReG,N,Uf)

ages          = size(x{1},1) - 1;
if isequal(char(table2array(ReG(1,4))),'q')
    F   = table(cellstr(char(kron(ones(ages,1),'q'))),cellstr(char(kron(ones(ages,1),char(x{2}(1))))),cellstr(char(x{2}(2:end))),ones(ages,1)/ages);
    LAB = "q(0, " + x{2}(2:end) + ")";
elseif isequal(char(table2array(ReG(1,4))),'m')
    F = table(cellstr(char(kron(ones(ages,1),'m'))),cellstr(char(x{2}(1:end - 1))),cellstr(char(x{2}(2:end))),ones(ages,1)/ages);
    LAB = "m(" + x{2}(1:end - 1) + ", " + x{2}(2:end) + ")";
end

CO            = string(ReG{1,1}) + "(" + string(ReG{1,2}) + ", " + string(ReG{1,3}) + ")";
sEL           = ~ismember(LAB,CO);
Z             = table2array(ReG(1,5));
G             = ReG(1,1:3);
inputs        = Set(q,x,{F,G},[]);
y             = log(inputs{1});
sET           = find(max(isinf(y)')' == 0);
y             = y(sET,sEL)';
X             = log(inputs{2}(sET,:))';
X             = X.^((0:Z)');
if isequal(numel(N),0)
    f = 1;
    N = size(y,2);
else
    f = N/size(y,2);
end

clc;
[~,s,~]       = svd(y*y');

R             = max(sum(cumsum(diag(s))/sum(diag(s)) < (1 - 10^-8)),Z + 1);
[U,s,v]       = svds(y,R);
Y             = (v*s)';             

ages          = size(Y,1);
YY            = Y*Y'*f;
XX            = X*X'*f;
XY            = X*Y'*f;
theta         = XX\XY;
if isequal(numel(Uf),0)
    EE       = YY - XY'*theta - theta'*XY + theta'*XX*theta;
    [Uf,~,~] = svds(U*EE*U',R);
else
    Uf       = Uf(:,1:R);
end

I             = sparse(eye(ages));
priorB        = {zeros(1 + Z,ages),ones(1 + Z,ages)*10^-4};
priorPsi      = {ones(ages,1)*10^-2,ones(ages,1)*10^-2};
priorT        = {triu(NaN(ages)) + tril(zeros(ages),-1),triu(NaN(ages)) + tril(ones(ages),-1)*10^-2};
G             = sum(~isnan(priorT{1}));
G             = G(G > 0);
lisT          = find(~isnan(priorT{2}));

C             = sparse(zeros(ages^2,numel(lisT)));
for i = 1:numel(lisT)
    C(lisT(i),i) = 1;
end
B             = theta;
Psi           = sparse(diag(gamrnd(priorPsi{1},priorPsi{2})));
T             = eye(ages);

rng(0);
IterAT        = 5000;
BurnIn        = 500;
for i = 1:BurnIn + IterAT
    EE        = YY - XY'*B - B'*XY + B'*XX*B;
    for j = 1:ages
        Psi(j,j)   = gamrnd(priorPsi{1}(j) + (N)/2,1/(priorPsi{2}(j) + 1/2*T(:,j)'*EE*T(:,j)));      
    end        
    
    L         = mat2cell(C'*kron(Psi,EE)*C + sparse(diag(priorT{2}(lisT))),G,G);
    V         = mat2cell(zeros(size(C,2)),G,G);
    for j = 1:numel(G)
        V{j,j} = sparse(pinv(full(L{j,j})));
        L{j,j} = chol(L{j,j})';        
    end

    L         = cell2mat(L);
    V         = cell2mat(V);
    Mu        = -C'*kron(Psi,EE)*reshape(I,[],1) + sparse(diag(priorT{2}(lisT))*priorT{1}(lisT));
    T(lisT)   = V*Mu + L'\normrnd(0,1,numel(lisT),1);
    Phi       = T*Psi*T';
    clear L V Mu

    V         = kron(Phi,XX) + diag(reshape(priorB{2},[],1));
    L         = chol(sparse(V))';
    V         = pinv(V);
    Mu        = kron(Phi,XX)*reshape(theta,[],1) + diag(reshape(priorB{2},[],1))*reshape(priorB{1},[],1);
    B         = V*Mu + L'\normrnd(0,1,numel(B),1);
    B         = reshape(B',Z + 1,ages);
    sB        = U*B';
    b         = NaN(numel(sEL),Z + 1);
    b(sEL,:)  = sB;
    b(~sEL,:) = ((0:Z) == 1); 
    clear L V Mu sB

    [u,~,~]   = svds(U*pinv(full(Phi))*U',R);
    su        = u.*sign(sum(Uf.*u));
    u         = NaN(numel(sEL),ages);
    u(sEL,:)  = su;
    u(~sEL,:) = 0;

    u         = [b u];
    %k         = (y - u(:,1))'*u(:,2:end);

    if i > BurnIn
        for j = 1:size(u,2)
            beta{j}(:,i - BurnIn) = u(:,j);
        end
        %K = [K;k];
    else
        %K = [];
    end
end