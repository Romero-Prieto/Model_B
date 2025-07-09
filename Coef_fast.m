function beta = Coef_fast(q,x,ReG)

ages    = size(x{1},1) - 1;
if isequal(char(table2array(ReG(1,4))),'q')
    F = table(cellstr(char(kron(ones(ages,1),'q'))),cellstr(char(kron(ones(ages,1),char(x{2}(1))))),cellstr(char(x{2}(2:end))),ones(ages,1)/ages);
elseif isequal(char(table2array(ReG(1,4))),'m')
    F = table(cellstr(char(kron(ones(ages,1),'m'))),cellstr(char(x{2}(1:end - 1))),cellstr(char(x{2}(2:end))),ones(ages,1)/ages);
end
Z       = table2array(ReG(1,5));
G       = ReG(1,1:3);
inputs  = Set(q,x,{F,G},[]);
y       = log(inputs{1});
sET     = find(max(isinf(y)')' == 0);
y       = y(sET,:);
X       = log(inputs{2}(sET,:));
X       = X.^(0:Z);
B       = pinv(X'*X)*X'*y;
E       = y - X*B;
[u,s,~] = svd(E'*E);
s       = diag(s);
beta    = [B',u];
s(1:5)/sum(s)