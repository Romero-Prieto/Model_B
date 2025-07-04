clear
pATh                    = "/Users/lshjr3/Documents/DHS/OuTPuT/";
load(char(pATh + "ReSuLTs/BoOTStrAp_DHS.mat"),'LT','iNFo','xE');
fileattrib(pATh + "ReSuLTs")
RESolUTioN              = 300;

sEL                     = (iNFo.R == '1. Model B' | iNFo.R == '2. Rest of the World');
LT                      = LT{3}(sEL,:);
iNFo                    = iNFo(sEL,:);
sEL                     = (iNFo.R == '1. Model B');

Z                       = 2;
ReG                     = table({'q'},{'0'},{'60m'},{'q'},Z);
[B,Q,~,x]               = Coef(ReG);
B                       = B{3};
sET                     = find(ismember(xE{2},x{2}));
I                       = zeros(sET(end - 1),numel(sET) - 1);
for i = 2:numel(sET)
    I(sET(i - 1):sET(i) - 1,i - 1) = 1;
end

F                       = table([],[],[],ones(0,1));
G                       = table({'q';'q'},{'0';'12m'},{'28d';'60m'});
sET                     = {F,ReG(1,1:3)};
Qm{1}                   = Q{3}';
Beta{1}                 = B;

p                       = [[4 10]/365.25;[2 4]/12;[5 7]/12;[8 10]/12;[11 24]/12];
P                       = {find(ismember(xE{2},["0";"1d";"2d";"3d";"7d";"14d";"28d";"3m";"6m";"12m";"24m";"60m"])),find(ismember(xE{2},["0";"6m";"12m";"18m";"24m";"36m";"48m";"60m"]))};
sR                      = find(ismember(xE{2},["0";"28d";"12m";"24m";"60m"]));

nMx                     = cell2mat(LT(:,1))';
nMx                     = cell2mat(mat2cell(prctile(nMx(2:end,:),50)',ones(numel(sEL),1)*size(LT{1,1},1))');
q                       = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(xE{1},1).*nMx,1))];
qR                      = q;
[qS,nMxS]               = smoothingLT(q,xE,sR,100,10,0.25,p,P);

ndx                     = q(2:end,:) - q(1:end - 1,:);
nLx                     = min(ndx./nMx,diff(xE{1},1).*(1 - q(1:end - 1,:)));
nMx                     = (I'*ndx)./(I'*nLx);
qI                      = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(x{1},1).*nMx,1))];
Qm{2}                   = qI(:,sEL);

ndx                     = qS(2:end,:) - qS(1:end - 1,:);
nLx                     = min(ndx./nMxS,diff(xE{1},1).*(1 - qS(1:end - 1,:)));
nMxS                    = (I'*ndx)./(I'*nLx);
qS                      = 1 - [ones(1,size(nMxS,2));exp(-cumsum(diff(x{1},1).*nMxS,1))];
B{1}                    = Coef_fast(qS(:,sEL),x,ReG);
Beta{2}                 = B;

iN                      = iNFo(sEL,:); 
J                       = [10 23 50 240];
for j = 1:numel(J)
    smoothingLT(q(:,J(j)),xE,sR,100,10,0.25,p,P);
    exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/smoothing example " + iN.country(J(j)) + " " + iN.fILe(J(j)) + ".png"),'Resolution',RESolUTioN);
end
clear iN J

for i = 1:numel(Beta)
    q          = Qm{2};
    X          = Set(q,x,{F,ReG(1,1:3)});
    X          = log(X{2}).^(0:Z);
    qm         = [zeros(1,size(q,2));exp(X*Beta{i}{1}(:,1:Z + 1)')'];
    
    q          = Set(q,x,{F,G});
    qm         = Set(qm,x,{F,G});
    sCaTTeR{i} = q{2}./qm{2};
    clear q qm X
end

H                       = [25 50 100]';
K                       = (-1.0:0.5:1.0)';
X                       = [kron(log(H/1000),ones(size(K))),kron(ones(size(H)),K)];
X                       = [X(:,1).^(0:Z),X(:,2)];
for i = 1:numel(Beta)
    Xs    = X;
    for j = 1:2
        q         = [zeros(1,size(Xs,1));exp(Beta{i}{1}(:,1:Z + 2)*Xs')];
        q         = Set(q,x,{F,G});
        qs{j,i}   = q{2};
        Xs(:,end) = 0;
    end
    clear q Xs
end

for i = 1:numel(Beta)
    OuT{1,i} = mat2cell(qs{1,i}./qs{2,1},ones(size(H))*numel(K),2);
    OuT{2,i} = mat2cell(qs{1,i}./qs{2,2},ones(size(H))*numel(K),2);
end


mPIX                     = 538756;
pix                      = 1/37.7952755906;
z                        = min(sqrt(mPIX/(10*18/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*[0 0 2*9 10]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};

for i = 1:2
    nexttile(i)
    ax{i}                            = gca;
    ax{i}.FontName                   = 'Times New Roman';
    ax{i}.FontSize                   = 9*z;
    ax{i}.XAxis.TickLabelFormat      = '%.1f';
    ax{i}.YAxis.TickLabelFormat      = '%.1f';
    ax{i}.YAxis.Exponent             = 0;

    ax{i}.XAxis.TickValues           = 0.00:0.50:5.00;
    ax{i}.XAxis.MinorTickValues      = 0.00:0.25:5.00;
    ax{i}.YAxis.TickValues           = 0.00:0.50:5.00;
    ax{i}.YAxis.MinorTickValues      = 0.00:0.25:5.00;
    ax{i}.XAxis.MinorTick            = 'on';
    ax{i}.YAxis.MinorTick            = 'on';
    ax{i}.XTickLabelRotation         = 0;
    ax{i}.YTickLabelRotation         = 0;
    ax{i}.LabelFontSizeMultiplier    = 1;

    xlim([0 2])
    ylim([0 2.5])
    
    if ~isequal(mod(i,3),2)
        xlabel('$\mathit{q}\mathrm{(}\mathrm{28}\mathit{d}\mathrm{)}/\mathit{q^A}\mathrm{(}\mathrm{28}\mathit{d;k=0}\mathrm{)}$','Interpreter','latex','FontSize',9*z);
        ylabel('$\mathit{q}\mathrm{(}\mathrm{12}\mathit{m,}\mathrm{5}\mathit{y}\mathrm{)}/\mathit{q^A}\mathrm{(}\mathrm{12}\mathit{m,}\mathrm{5}\mathit{y;k=0}\mathrm{)}$','Interpreter','latex','FontSize',9*z);
    else
        xlabel('$\mathit{q}\mathrm{(}\mathrm{28}\mathit{d}\mathrm{)}/\mathit{q^B}\mathrm{(}\mathrm{28}\mathit{d;k=0}\mathrm{)}$','Interpreter','latex','FontSize',9*z);
        ylabel('$\mathit{q}\mathrm{(}\mathrm{12}\mathit{m,}\mathrm{5}\mathit{y}\mathrm{)}/\mathit{q^B}\mathrm{(}\mathrm{12}\mathit{m,}\mathrm{5}\mathit{y;k=0}\mathrm{)}$','Interpreter','latex','FontSize',9*z);
    end 
    
    grid on;
    grid minor;
    box on;
    hold on;
end

nexttile(1)
for i = 1:numel(H)
    plot(OuT{1,1}{i}(:,1),OuT{1,1}{i}(:,2),'LineWidth',1.5,'color',coloR{i});
end
text(OuT{1,2}{3}(1,1),OuT{1,2}{3}(1,2),"$\mathit{k=}\mathrm{" + string(K(1)) + "}$",'Interpreter','latex','HorizontalAlignment','left','FontSize',10);
text(OuT{1,2}{3}(end,1),OuT{1,2}{3}(end,2),"$\mathit{k=}\mathrm{" + string(K(end)) + "}$",'Interpreter','latex','HorizontalAlignment','right','FontSize',10);
plot([-1 3],[1 1],'LineWidth',1.0,'color','k');
plot([1 1],[-1 3],'LineWidth',1.0,'color','k');
scatter(sCaTTeR{1}(:,1),sCaTTeR{1}(:,2),7.5,'filled','MarkerFaceColor',coloR{4},'MarkerFaceAlpha',.5)
legend("$\textrm{U5M = " + string(H) + "}$",'Interpreter','latex','FontSize',9*z,'FontAngle','oblique','Location','southoutside','NumColumns',4,'Box','off');
title('$\textrm{Model A (Log-quadratic)}$','Interpreter','latex','FontSize',11);

nexttile(2)
for i = 1:numel(H)
    plot(OuT{2,2}{i}(:,1),OuT{2,2}{i}(:,2),'LineWidth',1.5,'color',coloR{i});
end
text(OuT{2,2}{3}(1,1),OuT{2,2}{3}(1,2),"$\mathit{k=}\mathrm{" + string(K(1)) + "}$",'Interpreter','latex','HorizontalAlignment','left','FontSize',10*z);
text(OuT{2,2}{3}(end,1),OuT{2,2}{3}(end,2),"$\mathit{k=}\mathrm{" + string(K(end)) + "}$",'Interpreter','latex','HorizontalAlignment','right','FontSize',10*z);
plot([-1 3],[1 1],'LineWidth',1.0,'color','k');
plot([1 1],[-1 3],'LineWidth',1.0,'color','k');
scatter(sCaTTeR{2}(:,1),sCaTTeR{2}(:,2),7.5,'filled','MarkerFaceColor',coloR{4},'MarkerFaceAlpha',.5)
title('$\textrm{Model B}$','Interpreter','latex','FontSize',11*z);
exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/scatters.png"),'Resolution',RESolUTioN);



QM{1}                   = qI;
QM{2}                   = qI;
QM{3}                   = qI;
Beta                    = Beta(1:2);

iNFo.R(sEL)             = "Model B";
iNFo.R(~sEL)            = "Rest of the World";
Region                  = tabulate(iNFo.R);
Region                  = sortrows(string(Region(:,1)));
for i = 1:numel(Region)
    group{1}{i}  = find(iNFo.R == Region(i));
    lEGeND{1}{i} = char("$\textrm{" + Region(i) +  " (" + numel(group{1}{i}) + ")}$");
end

iNFo.R(sEL)       = iNFo.SubRegion(sEL);
Region            = tabulate(iNFo.R(sEL));
Region            = sortrows(string(Region(:,1)));
for i = 1:numel(Region)
    group{2}{i}  = find(iNFo.R == Region(i));
    lEGeND{2}{i} = char("$\textrm{" + Region(i) +  " (" + numel(group{2}{i}) + ")}$");
end

Region            = tabulate(iNFo.SubRegion);
Region            = sortrows(string(Region(:,1)));
for i = 1:numel(Region)
    group{3}{i}  = find(iNFo.SubRegion == Region(i));
    lEGeND{3}{i} = char("$\textrm{" + Region(i) +  " (" + numel(group{3}{i}) + ")}$");
end
models            = {char("$\textrm{Log-quadratic model (N = " + size(Qm{1},2) + ")}$"),char("$\textrm{Model B (N = " + size(Qm{2},2) + ")}$"),char("$\textrm{Model B2 (N = " + size(Qm{2},2) + ")}$")};

ages              = size(x{1},1) - 1;
F                 = table(cellstr(char(kron(ones(ages,1),'q'))),cellstr(char(kron(ones(ages,1),char(x{2}(1))))),cellstr(char(x{2}(2:end))),diff(x{1})/x{1}(end));
G                 = ReG(1,1:3);

g                 = 1000;
H                 = log(.1./(2.^[6 -3]));
H                 = H(1):(H(2) - H(1))/(g - 1):H(2);
K                 = [-3 3];
K                 = K(1):(K(2) - K(1))/(g - 1):K(2);
xo                = [kron(H',ones(g,1)),kron(ones(g,1),K')];

for i = 1:numel(Beta)
    [~,qs]   = Model(Beta{i},xo,x,{F,[]},{NaN(g^2,size(F,1)),NaN(g^2,0)},zeros(g^2,0));
    rever{i} = reshape(any(qs(2:end,:) <= qs(1:end - 1,:))',g,g);    
end



z                        = min(sqrt(mPIX/(numel(Beta)*(numel(group) + 1)*9^2/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*9*[0 0 numel(Beta) (numel(group) + 1)]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(numel(group),numel(Beta),'Padding','compact','TileSpacing','compact');
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};
coloR                    = [coloR coloR coloR]
title(TL,'$\textrm{Undetermined Outcomes}$','Interpreter','latex','FontSize',12*z);

for i = 1:numel(group)*numel(Beta)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8*z;
    
    ax{i}.XScale                = 'log';
    ax{i}.XAxis.TickLabelFormat = '%.1g';
    ax{i}.XAxis.TickValues      = .1./(2.^(7:-1:-3))*1000;
    ax{i}.XAxis.MinorTickValues = 10:10:400;
    ax{i}.XAxis.MinorTick       = 'off';
    xlim([6.25 400])
    
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.YAxis.TickValues      = -3:1.00:3;
    ax{i}.YAxis.MinorTickValues = -3:0.50:3;
    ax{i}.YAxis.MinorTick       = 'off';        
    ylim([-3 3])
    
    if isequal(1 + mod(i - 1,numel(Beta)),1)
        ylabel('$\textrm{Best fitting \textit{k}}$','Interpreter','latex','FontSize',9*z);
    end
    if isequal(ceil(i/numel(Beta)),numel(group))
        xlabel('$\mathrm{_5}\mathit{q}\mathrm{_0}$ $\textit{deaths per 1,000 births}$ $\textrm{(log scale)}$','Interpreter','latex','FontSize',9*z);
    end
    if isequal(ceil(i/numel(Beta)),1)
        title(models{1 + mod(i - 1,numel(Beta))},'Interpreter','latex','FontSize',9*z);
    end
    
    grid on;
    grid minor;
    box on;
    hold on;
end

for i = 1:numel(group)*numel(Beta)
    nexttile(i)
    h          = 1 + mod(i - 1,numel(Beta));
    r          = ceil(i/numel(Beta));
    imagesc(exp(H)*1000,K,rever{h},'alphaData',.25);
    colormap([[1 1 1];[0 0 0]]);
    
    xo         = Set(QM{r},x,sET);
    xo         = [log(xo{2}),zeros(size(xo{2}))];
    match      = Set(QM{r},x,{F,G});
    xa         = Match(match,Beta{h},x,{F,G},xo,50);
    
    for j = 1:numel(group{r})
        scatter(exp(xa(group{r}{j},1))*1000,xa(group{r}{j},2),15*z,'filled','MarkerFaceColor',coloR{j},'MarkerFaceAlpha',0.5,'MarkerEdgeColor',coloR{j},'MarkerEdgeAlpha',0.0);
    end
    for j = 1:numel(group{r})
        plot([6.25 400],ones(1,2)*mean(xa(group{r}{j},2)),'color',coloR{j},'LineWidth',1.25*z);
    end
    
    if h == 1
        legend(lEGeND{r},'Interpreter','latex','FontSize',8*z,'FontAngle','oblique','Location','southoutside','NumColumns',2,'Box','off');
    end
end
exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/undetermined outcomes.png"),'Resolution',RESolUTioN);


group{4}          = {(1:size(Qm{1},2))'};
lEGeND{4}         = {char("$\textrm{Log-quadratic Model: U5MD (" + numel(group{4}{1}) + ")}$")};
QM{4}             = Qm{1};
B                 = {Beta{2},Beta{1}};

W                 = diff(x{1},1)/x{1}(end);
Fs                = table([],[],[],ones(0,1));
G                 = table({'q';'q';'q';'q';'q';'q';'q';'k2'},{'0';'0';'0';'0';'0';'0';'28d';'0'},{'60m';'7d';'28d';'3m';'6m';'12m';'60m';'60m'});

mOD{1}            = {Fs,G([1 8],:)};
mOD{2}            = {Fs,G([1 2],:)};
mOD{3}            = {Fs,G([1 3],:)};
mOD{4}            = {Fs,G([1 4],:)};
mOD{5}            = {Fs,G([1 5],:)};
mOD{6}            = {Fs,G([1 6],:)};
mOD{7}            = {F,G(1,:)};
mOD{8}            = {Fs,G([7 8],:)};

pOIntS{2,1}       = {'$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)\,and\,}\mathit{q}\mathrm{(7}\mathit{d}\mathrm{)}$'};
pOIntS{3,1}       = {'$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)\,and\,}\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$'};
pOIntS{4,1}       = {'$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)\,and\,}\mathit{q}\mathrm{(3}\mathit{m}\mathrm{)}$'};
pOIntS{5,1}       = {'$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)\,and\,}\mathit{q}\mathrm{(6}\mathit{m}\mathrm{)}$'};
pOIntS{6,1}       = {'$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)\,and\,}\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$'};
pOIntS{7,1}       = {'$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)\,and\,all\,}\mathit{q}\mathrm{(}\mathit{x}\mathrm{)}$'};

lABs              = {{1} {2 3 4 5 6 7} {8}};
vARs              = {'$\mathrm{All\,}\mathit{q}\mathrm{(}\mathit{x}\mathrm{)}$','$\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','$\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)}$'};
vARs              = {vARs vARs vARs};
foRMaT            = {'%0.4f','%0.4f','%0.4f','%0.4f','%+0.4f','%+0.4f','%+0.4f','%+0.4f','%0.2f','%0.2f','%0.2f','%0.2f'};
sEt               = {'$\textrm{RMSE}$','$\mathit{\mu}\textrm{ (bias)}$','$\mathit{\varphi}\textrm{ (precision)}$'};

q                 = QM{1};
init              = Set(q,x,{F,ReG});
init              = [log(init{2}),zeros(size(init{2}))];
match             = Set(q,x,{F,G(1,:)});
xa                = Match(match,B{1},x,{F,G(1,:)},init,250);
s                 = group{1}{1};
k(s,1)            = xa(s,2);
xa                = Match(match,B{2},x,{F,G(1,:)},init,250);
s                 = group{1}{2};
k(s,1)            = xa(s,2);
k(:,2)            = NaN;

Region            = tabulate(iNFo.SubRegion);
Region            = sortrows(string(Region(:,1)));
h                 = NaN(numel(Region),1);
taBLeK            = table(Region,h,h,h,'VariableNames',{'Region','p50','LB','UB'});
for h = 1:numel(Region)
    s                 = (iNFo.SubRegion == Region(h));
    k(s,2)            = round(median(k(s,1)),4);
    k(s,3)            = round(mean(k(s,1)),4);
    taBLeK.p50(h)     = round(prctile(k(s,1),50),4);
    taBLeK.LB(h)      = round(prctile(k(s,1),2.5),4);
    taBLeK.UB(h)      = round(prctile(k(s,1),97.5),4);
end


s                 = group{1}{2};
k(s,2:3)          = 0;
lISt              = ["smooth";"raw"];
for h = 1:numel(lISt)
    for i = 1:numel(group{1})

        if isequal(lISt(h),'smooth')
            q                = qS(:,group{1}{i});
        elseif isequal(lISt(h),'raw')
            q                = QM{1}(:,group{1}{i});
        end
            
        init              = Set(q,x,{F,ReG});
        init              = [log(init{2}),zeros(size(init{2}))];
        
        if ~isequal(std(k(group{1}{i},2)),0)
            pOIntS{1,1}      = {'$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)\,only,\,}\mathit{k=\,}\mathrm{regional\,median}$'};
            pOIntS{8,1}      = {'$\mathit{q}\mathrm{(28}\mathit{d}\mathrm{,5}\mathit{y}\mathrm{)\,only,\,}\mathit{k=\,}\mathrm{regional\,median}$'};
        else
            pOIntS{1,1}      = {char("$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)\,only,\,}\mathit{k=\,}\mathrm{" + sprintf('%0.4f',k(group{1}{i}(1),2)) + "}$")};
            pOIntS{8,1}      = {char("$\mathit{q}\mathrm{(28}\mathit{d}\mathrm{,5}\mathit{y}\mathrm{)\,only,\,}\mathit{k=\,}\mathrm{" + sprintf('%0.4f',k(group{1}{i}(1),2)) + "}$")};
        end
            
        for j = 1:numel(mOD)
            if ismember(j,[1 8])
                init(:,2) = k(group{1}{i},2);
            else
                init(:,2) = 0;
            end
    
            match            = Set(q,x,mOD{j},init);
            [~,qp]           = Match(match,B{i},x,mOD{j},init,250);
            SE               = ismember(F(:,1:end - 1),mOD{j}{2});
            match            = Set(q,x,{F,[]});
            mATcH            = Set(qp,x,{F,[]});
            y                = log(match{1}(:,~SE)./mATcH{1}(:,~SE));
            
            [N,R]            = size(y);
            priorPhi         = {ones(1,R)*10^-2,ones(1,R)*10^-2};
            X                = ones(N,1);
            XX               = X'*X;
            Xy               = X'*y;
            yy               = y'*y;
            priorBi          = {pinv(XX)*Xy;ones(1,R)*10^-4};
    
            Bi               = zeros(1,R);
            Phi              = priorPhi{1};
            warmup           = 2500;
            iterations       = 10000;
    
            sample.Bi        = NaN(iterations,numel(SE));
            sample.Bi(:,SE)  = 0;
            sample.Phi       = NaN(iterations,numel(SE));
            sample.Phi(:,SE) = Inf;
            sample.MSE       = NaN(iterations,numel(SE));
            sample.MSE(:,SE) = 0;
            for r = 1:iterations + warmup
                EE        = diag(yy - Xy'*Bi - Bi'*Xy + Bi'*XX*Bi)';
                V         = 1./(Phi*XX + priorBi{2});
                M         = V.*(Phi.*Xy + priorBi{2}.*priorBi{1});
                Bi        = mvnrnd(M,V);
                Phi       = gamrnd(priorPhi{1} + N/2,1./(priorPhi{2} + 1/2*EE));
                
                if r > warmup
                    sample.Bi(r - warmup,~SE)  = Bi;
                    sample.Phi(r - warmup,~SE) = Phi;
                    sample.MSE(r - warmup,~SE) = Bi.^2 + 1./Phi;
                end
            end
            tABle{i}.RMSE(j,:)      = [prctile(sqrt(sample.MSE*W),[50 2.5 97.5]) mat2cell(prctile(sqrt(sample.MSE(:,[4 15 22])),[50 2.5 97.5])',ones(3,1),3)'];
            tABle{i}.Bias(j,:)      = [prctile(sample.Bi*W,[50 2.5 97.5]) mat2cell(prctile(sample.Bi(:,[4 15 22]),[50 2.5 97.5])',ones(3,1),3)'];
            tABle{i}.Precision(j,:) = [prctile(1./(1./sample.Phi*W),[50 2.5 97.5]) mat2cell(prctile(sample.Phi(:,[4 15 22]),[50 2.5 97.5])',ones(3,1),3)'];
            clear sample 
        end
    
        tABle{i}.Precision     = cell2mat(tABle{i}.Precision);
        tABle{i}.Precision     = mat2cell(tABle{i}.Precision,ones(size(tABle{i}.Bias,1),1),ones(1,size(tABle{i}.Bias,2))*3);
        nOTe                   = {'$\textrm{Entry points}$',char("$\mathrm{p50/}\mathit{[p2.5,p97.5]} \textrm{ of the posterior predictive distribution. Gibbs sampler algorithm, iterated 10,000 times after 2,500 rounds of warm-up. Not significant at 5\%, in pink.}\textbf{ " + lEGeND{1}{i} + ".}$")};
        nOTe                   = {'$\textrm{Entry points}$',''};
        tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,pOIntS,cell2mat([tABle{i}.RMSE tABle{i}.Bias tABle{i}.Precision]),0.145,0.075,[]);
        exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/RMSE " + lISt(h) + " (Group " + i + ").png"),'Resolution',RESolUTioN);
        
        nOTe                   = {'$\textrm{Entry points}$',char("$\mathrm{p50/}\mathit{[p2.5,p97.5]} \textrm{ of the posterior predictive distribution.}\textbf{ " + lEGeND{1}{i} + ".}$")};
        nOTe                   = {'$\textrm{Entry points}$',''};
        tABleBAyEs(sEt,vARs(1),foRMaT,lABs,nOTe,pOIntS,cell2mat(tABle{i}.RMSE),0.145,0.075,[]);
        exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/RMSE " + lISt(h) + " (Group " + i + " short).png"),'Resolution',RESolUTioN);
    end
end


foRMaT            = {'%0.4f','%0.4f','%0.4f','%0.4f','%0.4f','%0.4f','%0.4f','%0.4f','%+0.3f','%+0.3f','%+0.3f','%+0.3f','%+0.3f','%+0.3f','%+0.3f','%+0.3f','%0.0f','%0.0f','%0.0f','%0.0f','%0.0f','%0.0f','%0.0f','%0.0f'};
vARs              = {'$\mathrm{All\,}\mathit{q}\mathrm{(}\mathit{x}\mathrm{)}$','$\mathit{q}\mathrm{(7}\mathit{d}\mathrm{)}$','$\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','$\mathit{q}\mathrm{(3}\mathit{m}\mathrm{)}$','$\mathit{q}\mathrm{(6}\mathit{m}\mathrm{)}$','$\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','$\mathit{q}\mathrm{(48}\mathit{m}\mathrm{)}$','$\mathit{q}\mathrm{(5}\mathit{y}\mathrm{)}$'};
vARs              = {vARs vARs vARs};

lambda            = [0.025 0.05 0.10 0.25 0.50 0.75];
knots             = [5 10 15 20];

for i = 1:numel(group{1})
    q                 = QM{1}(:,group{1}{i});
    s                 = 1;
    for h = 1:numel(knots)
        for j = 1:numel(lambda)
            [qS,nMxS]       = smoothingLT(qR(:,group{1}{i}),xE,sR,100,knots(h),lambda(j),p,P);
            ndx             = qS(2:end,:) - qS(1:end - 1,:);
            nLx             = min(ndx./nMxS,diff(xE{1},1).*(1 - qS(1:end - 1,:)));
            nMxS            = (I'*ndx)./(I'*nLx);
            qS              = 1 - [ones(1,size(nMxS,2));exp(-cumsum(diff(x{1},1).*nMxS,1))];
            Qs{s}           = qS;
            clear qS nMxS ndx nLx
            if isequal(i,1)
                pARAmeTeRS{s,1} = {char("$\mathit{\lambda=\,}\mathrm{" + sprintf('%0.3f',lambda(j)) + "\,,\," + knots(h) + "\,}\mathit{kn.}$")};
                Bs{s}           = Coef_fast(Qs{s},x,ReG);
            end
            lABsII{h}{j}(1) = s;
            s               = s + 1;
        end
    end
    
    if isequal(i,1)
        lABsII{h + 1}{1} = s;
        pARAmeTeRS{s,1}  = {'$\textrm{No smoothing}$'};
        Bs{s}            = Coef_fast(q,x,ReG);
    end
    
    for j = 1:numel(Qs)
        y                = log(q(2:end,:)'./Qs{j}(2:end,:)');
        SE               = (std(y) == 0);
        y                = y(:,~SE);
        
        [N,R]            = size(y);
        priorPhi         = {ones(1,R)*10^-2,ones(1,R)*10^-2};
        X                = ones(N,1);
        XX               = X'*X;
        Xy               = X'*y;
        yy               = y'*y;
        priorBi          = {pinv(XX)*Xy;ones(1,R)*10^-4};

        Bi               = zeros(1,R);
        Phi              = priorPhi{1};
        warmup           = 2500;
        iterations       = 10000;

        sample.Bi        = NaN(iterations,numel(SE));
        sample.Bi(:,SE)  = 0;
        sample.Phi       = NaN(iterations,numel(SE));
        sample.Phi(:,SE) = Inf;
        sample.MSE       = NaN(iterations,numel(SE));
        sample.MSE(:,SE) = 0;
        for r = 1:iterations + warmup
            EE        = diag(yy - Xy'*Bi - Bi'*Xy + Bi'*XX*Bi)';
            V         = 1./(Phi*XX + priorBi{2});
            M         = V.*(Phi.*Xy + priorBi{2}.*priorBi{1});
            Bi        = mvnrnd(M,V);
            Phi       = gamrnd(priorPhi{1} + N/2,1./(priorPhi{2} + 1/2*EE));
            
            if r > warmup
                sample.Bi(r - warmup,~SE)  = Bi;
                sample.Phi(r - warmup,~SE) = Phi;
                sample.MSE(r - warmup,~SE) = Bi.^2 + 1./Phi;
            end
        end
        tABleII{i}.RMSE(j,:)      = [prctile(sqrt(sample.MSE*W),[50 2.5 97.5]) mat2cell(prctile(sqrt(sample.MSE(:,[1 4 6 9 15 21 22])),[50 2.5 97.5])',ones(7,1),3)'];
        tABleII{i}.Bias(j,:)      = [prctile(sample.Bi*W,[50 2.5 97.5]) mat2cell(prctile(sample.Bi(:,[1 4 6 9 15 21 22]),[50 2.5 97.5])',ones(7,1),3)'];
        tABleII{i}.Precision(j,:) = [prctile(1./(1./sample.Phi*W),[50 2.5 97.5]) mat2cell(prctile(sample.Phi(:,[1 4 6 9 15 21 22]),[50 2.5 97.5])',ones(7,1),3)'];
        clear sample 
    end

    tABleII{i}.Precision     = cell2mat(tABleII{i}.Precision);
    tABleII{i}.Precision     = mat2cell(tABleII{i}.Precision,ones(size(tABleII{i}.Bias,1),1),ones(1,size(tABleII{i}.Bias,2))*3);
    nOTe                     = {'$\textrm{Parameters}$',char("$\mathrm{p50/}\mathit{[p2.5,p97.5]} \textrm{ of the posterior predictive distribution. Gibbs sampler algorithm, iterated 10,000 times after 2,500 rounds of warm-up. Not significant at 5\%, in pink.}\textbf{ " + lEGeND{1}{i} + ".}$")};
    nOTe                     = {'$\textrm{Parameters}$',''};
    tABleBAyEs(sEt,vARs,foRMaT,lABsII(1:end - 1),nOTe,pARAmeTeRS,cell2mat([tABleII{i}.RMSE tABleII{i}.Bias tABleII{i}.Precision]),0.125,0.075,[]);
    exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/smoothing (Group " + i + ").png"),'Resolution',RESolUTioN);
    
    nOTe                     = {'$\textrm{Parameters}$',char("$\mathrm{p50/}\mathit{[p2.5,p97.5]} \textrm{ of the posterior predictive distribution.}\textbf{ " + lEGeND{1}{i} + ".}$")};
    nOTe                     = {'$\textrm{Parameters}$',''};
    tABleBAyEs(sEt,vARs(1),foRMaT,lABsII(1:end - 1),nOTe,pARAmeTeRS,cell2mat(tABleII{i}.RMSE),0.125,0.075,[]);
    exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/smoothing (Group " + i + " short).png"),'Resolution',RESolUTioN);
end



z                        = min(sqrt(mPIX/(numel(knots)*(numel(lambda) - 1)*9^2/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*9*[0 0 numel(knots) (numel(lambda) - 1)]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(numel(knots),numel(lambda),'Padding','compact','TileSpacing','compact');
mODlAB{1}                = '$\textrm{Model B (smooth)}$';
mODlAB{2}                = '$\textrm{Model B (smooth), }\mathit{k\,}\mathrm{= -1}$';
mODlAB{3}                = '$\textrm{Model B (smooth), }\mathit{k\,}\mathrm{= +1}$';
mODlAB{4}                = '$\textrm{Model B (raw)}$';
mODlAB{5}                = '$\textrm{Model B (raw), }\mathit{k\,}\mathrm{= -1}$';
mODlAB{6}                = '$\textrm{Model B (raw), }\mathit{k\,}\mathrm{= +1}$';
mODlAB{7}                = '$\textrm{Log-quadratic Model}$';
mODlAB{8}                = '$\textrm{Log-quadratic Model, }\mathit{k\,}\mathrm{= -1}$';
mODlAB{9}                = '$\textrm{Log-quadratic Model, }\mathit{k\,}\mathrm{= +1}$';

for i = 1:numel(Bs) - 1
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 10*z;
    ax{i}.XAxis.TickLabelFormat = '%.1f';
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    ax{i}.YScale                = 'log';
    ax{i}.XAxis.MinorTick       = 'off';
    ax{i}.YAxis.MinorTick       = 'off';
    ax{i}.XAxis.TickValues      = 0.00:1.00:5.00;
    ax{i}.XAxis.MinorTickValues = 0.00:0.50:5.00;
 
    if isequal(mod(i,numel(lambda)),1)
        ylabel(char("$\textrm{" + sprintf('%0.0f',knots(ceil(i/numel(lambda)))) +"}\textit{ knots}$"),'Interpreter','latex','FontSize',11*z);
    else
        ax{i}.YAxis.TickLabels = [];
    end
    
    if isequal(ceil(i/numel(lambda)),numel(knots))
        xlabel('$\mathit{age\,(in\,years)}$','Interpreter','latex','FontSize',11*z);
    else
        ax{i}.XAxis.TickLabels = [];
    end

    if i <= numel(lambda)
        title(char("$\mathrm{\lambda =" + sprintf('%0.3f',lambda(i)) +"}$"),'Interpreter','latex','FontSize',12*z);
    end
       
    grid on;
    grid minor;
    box on;
    hold on;
end

U5MR                     = .1./(2.^(2:-1:-2))*1000;
for j = 1:numel(U5MR)
    P       = cell(0);
    yLim    = [10 -10];
    qRW     = [zeros(1,3);exp(Bs{end}(:,1:(ReG.Z + 1))*log(U5MR(j)/1000).^(0:ReG.Z)' + Bs{end}(:,ReG.Z + 2)*[0 -1 1])];
    nMxRW   = -log((1 - qRW(2:end,:))./(1 - qRW(1:end - 1,:)))./diff(x{1},1);
    for h = 1:size(qRW,2)
        if min(nMxRW(:,h)) <= 0
            nMxRW(:,h) = NaN;
        end
    end
    s       = find(nMxRW <= 0);
    nMxR(s) = NaN;

    qv         = [zeros(1,3);exp(Beta{1}{1}(:,1:(ReG.Z + 1))*log(U5MR(j)/1000).^(0:ReG.Z)' + Beta{1}{1}(:,ReG.Z + 2)*[0 -1 1])];
    nMxv       = -log((1 - qv(2:end,:))./(1 - qv(1:end - 1,:)))./diff(x{1},1);
    for h = 1:size(qv,2)
        if min(nMxv(:,h)) <= 0
            nMxv(:,h) = NaN;
        end
    end
    
    for i = 1:numel(Bs) - 1
        nexttile(i)
        xlim([-0.1 5])
        qs         = [zeros(1,3);exp(Bs{i}(:,1:(ReG.Z + 1))*log(U5MR(j)/1000).^(0:ReG.Z)' + Bs{i}(:,ReG.Z + 2)*[0 -1 1])];
        nMxs       = -log((1 - qs(2:end,:))./(1 - qs(1:end - 1,:)))./diff(x{1},1);
        for h = 1:size(qs,2)
            if min(nMxs(:,h)) <= 0
                nMxs(:,h) = NaN;
            end
        end
 
        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxs(:,1),'color',coloR{7},'LineWidth',1.00*z);        
        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxs(:,2),'color',coloR{2},'LineWidth',0.50*z);
        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxs(:,3),'color',coloR{5},'LineWidth',0.50*z);
        
        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxRW(:,1),'color',[coloR{7} 0.25],'LineWidth',1.00*z);
        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxRW(:,2),'color',[coloR{2} 0.25],'LineWidth',0.50*z);
        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxRW(:,3),'color',[coloR{5} 0.25],'LineWidth',0.50*z);

        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxv(:,1),'color',coloR{1},'LineWidth',1.00*z);
        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxv(:,2),'color',coloR{4},'LineWidth',0.50*z);
        P{end + 1} = plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxv(:,3),'color',coloR{6},'LineWidth',0.50*z);
        yLim(1)    = min(min(nMxs(:,1)*0.90),yLim(1));
        yLim(2)    = max(max(nMxs(:,1)*1.10),yLim(2));

        if isequal(i,numel(lambda)*(numel(knots) - 1) + ceil(numel(lambda)/2))
            P{end + 1} = legend(mODlAB,'Interpreter','latex','FontSize',9*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
    for i = 1:numel(Bs) - 1
        nexttile(i)
        ylim(yLim)
    end
    title(TL,char("$\mathrm{U5MR=" + U5MR(j) + "}$"),'Interpreter','latex','FontSize',12*z);
    exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/smoothing nMx, U5MR = " + U5MR(j) + ".png"),'Resolution',RESolUTioN);

    for i = 1:numel(Bs) - 1
        nexttile(i)
        xlim([-0.025 1.25])
    end    
    exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/smoothing nMx, U5MR = " + U5MR(j) + " close.png"),'Resolution',RESolUTioN);
    
    for i = 1:numel(P)
        delete(P{i});
    end
end




load(char(pATh + "ReSuLTs/BoOTStrAp_DHS.mat"),'LT','iNFo','xE');
LT                      = LT{3};
sEL                     = (iNFo.R == '1. Model B');
nMx                     = cell2mat(LT(:,1))';
nMx                     = cell2mat(mat2cell(prctile(nMx(2:end,:),50)',ones(numel(sEL),1)*size(LT{1,1},1))');
q                       = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(xE{1},1).*nMx,1))];

ndx                     = q(2:end,:) - q(1:end - 1,:);
nLx                     = min(ndx./nMx,diff(xE{1},1).*(1 - q(1:end - 1,:)));
nMx                     = (I'*ndx)./(I'*nLx);
q                       = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(x{1},1).*nMx,1))];

y                       = log(q(2:end,:)) - Beta{2}{1}(:,1:Z + 1)*log(q(end,:)).^((0:Z)');
[N,R]                   = size(y);
X                       = Beta{2}{1}(:,Z + 2);
W                       = diag(diff(x{1},1)/x{1}(end));
XX                      = X'*W*X;
Xy                      = X'*W*y;
yy                      = y'*W*y;
priorBi                 = {pinv(XX)*Xy;ones(1,R)*10^-4};
priorPhi                = {ones(1,R)*10^-2,ones(1,R)*10^-2};

Bi                      = ones(1,R)*0;
Phi                     = ones(1,R)*10^-2;
warmup                  = 2500;
iterations              = 10000;
for r = 1:iterations + warmup
    EE       = diag(yy - Xy'*Bi - Bi'*Xy + Bi'*XX*Bi)';
    V        = 1./(Phi*XX + priorBi{2});
    M        = V.*(Phi.*Xy + priorBi{2}.*priorBi{1});
    Bi       = mvnrnd(M,V);
    Phi      = gamrnd(priorPhi{1} + N/2,1./(priorPhi{2} + 1/2*EE));

    if r > warmup
        sample.Bi(r - warmup,:)  = Bi;
        sample.Phi(r - warmup,:) = Phi;
        sample.MSE(r - warmup,:) = Bi.^2 + 1./Phi;
    end
end


z                        = min(sqrt(mPIX/(2*2*9^2/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*9*[0 0 2 2]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.05 0.05 0.05],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65]};

for i = 1:4
    if i == 1
        nexttile(i)
        ax{i}                            = gca;
        ax{i}.FontName                   = 'Times New Roman';
        ax{i}.FontSize                   = 9*z;
        ax{i}.XAxis.TickLabelInterpreter = 'latex';
        ax{i}.XAxis.TickLabelFormat      = '%.1f';
        ax{i}.YAxis.TickLabelFormat      = '%.2f';

        ax{i}.XAxis.TickValues           = 0.00:1.00:5.00;
        ax{i}.XAxis.MinorTickValues      = 0.00:0.25:5.00;
        %ax{i}.YAxis.TickValues           = 0.50:0.25:1.50;
        %ax{i}.YAxis.MinorTickValues      = 0.50:0.05:2.00;
        ax{i}.XAxis.MinorTick            = 'off';
        ax{i}.YAxis.MinorTick            = 'off';
        xlim([-0.1 5])
        %ylim([0.5 1.75])

        xlabel('$\mathit{age\,(in\,years)}$','Interpreter','latex','FontSize',10*z);
        ylabel('$\mathit{q}\mathrm{(}\mathrm{x}\mathrm{)}$','Interpreter','latex','FontSize',10*z);
        grid on;
        grid minor;
        box on;
        hold on;
    else
        nexttile(i)
        ax{i}                            = gca;
        ax{i}.FontName                   = 'Times New Roman';
        ax{i}.FontSize                   = 9*z;
        ax{i}.XAxis.TickLabelInterpreter = 'latex';
        ax{i}.XAxis.TickLabelFormat      = '%.1f';
        ax{i}.YAxis.TickLabelFormat      = '%.2f';
        ax{i}.YScale                     = 'log';
        
        if isequal(i,2)
            ax{i}.XAxis.TickValues      = 0.00:1.00:5.00;
            ax{i}.XAxis.MinorTickValues = 0.00:0.50:5.00;
            xlim([-0.1 5])
            xlabel('$\mathit{age\,(in\,years)}$','Interpreter','latex','FontSize',10*z);
        elseif isequal(i,3)
            ax{i}.XAxis.TickValues      = x{1}([1 2 5 6 7 8]);
            ax{i}.XAxis.TickLabels      = x{3}([1 2 5 6 7 8]);
            ax{i}.XAxis.MinorTickValues = x{1};
            xlim([-0.005 1/3])
            xlabel('$\mathit{age\,(in\,days\,and\,months)}$','Interpreter','latex','FontSize',10*z);
        elseif isequal(i,4)
            ax{i}.XAxis.TickValues      = x{1}(1:5);
            ax{i}.XAxis.TickLabels      = x{3}(1:5);
            ax{i}.XAxis.MinorTickValues = 0:28;
            xlim([-0.001 x{1}(5)])
            xlabel('$\mathit{age\,(in\,days)}$','Interpreter','latex','FontSize',10*z);
        end
        
        ax{i}.XAxis.MinorTick       = 'off';
        ax{i}.YAxis.MinorTick       = 'off';
        ylabel('$\mathit{_nM_x}$','Interpreter','latex','FontSize',11*z);
        grid on;
        grid minor;
        box on;
        hold on;
    end
end

lEGeND            = {'$\textrm{Log-quadratic model,}$ $\mathit{k = \pm 1}$','$\textrm{Model B, best fitting}$','$\textrm{p50 DHS}$','$\textit{Bootstrapping DHS}$'};
for j = 1:size(q,2)
    qm{1}             = [zeros(1,3);exp(Beta{1}{1}(:,1:Z + 1)*log(q(end,j)).^((0:Z)') + Beta{1}{1}(:,Z + 2)*[0 -1 1])];
    qm{2}             = [zeros(1,3);exp(Beta{2}{1}(:,1:Z + 1)*log(q(end,j)).^((0:Z)') + Beta{2}{1}(:,Z + 2)*prctile(sample.Bi(:,j),[50 2.5 97.5]))];
    qm{3}             = q(:,j);
    
    ndx               = LT{j,2}(2:end,:) - LT{j,2}(1:end - 1,:);
    nLx               = min(ndx./LT{j,1},diff(xE{1},1).*(1 - LT{j,2}(1:end - 1,:)));
    nMx               = (I'*ndx)./(I'*nLx);
    qm{4}             = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(x{1},1).*nMx,1))];
    
    for i = 1:numel(qm)
        mm{i} = -log((1 - qm{i}(2:end,:))./(1 - qm{i}(1:end - 1,:)))./diff(x{1},1);
        if ismember(i,[1 2])
            mm{i}(:,2:3) = [min(mm{i}(:,2:3),[],2) max(mm{i}(:,2:3),[],2)];
        end
    end
    
    if isequal(iNFo.sTArt(j),iNFo.sTArt(j))
       dir   = iNFo.fILe(j) + " " + iNFo.country(j) + " " + iNFo.sTArt(j);
    else
       dir   = iNFo.fILe(j) + " " + iNFo.country(j) + " " + iNFo.sTArt(j) + "-" + iNFo.eNd(j); 
    end
    if ismember(j,sEL)
        name = dir + ". Region: " + iNFo.Region(j) + " (" + iNFo.SubRegion(j) + ")" + ", U5MR = " + sprintf('%0.2f',q(end,j)*1000);
    else
        name = dir + ". Region: " + iNFo.Region(j) + ", U5MR = " + sprintf('%0.2f',q(end,j)*1000);
    end

    nexttile(1)
    WX{1}        = title(TL,name,'Interpreter','latex','FontSize',11*z);
    for i = 1:3
        WX{end + 1}  = plot(x{1},qm{i}(:,1),'LineWidth',1.25*z,'color',coloR{i});
    end
    for i = 1:2
        WX{end + 1}  = fill([x{1};flip(x{1})],[qm{i}(:,2);flip(qm{i}(:,3))],coloR{i},'FaceAlpha',.15,'EdgeAlpha',0.00,'LineWidth',0.25*z,'LineStyle','-','EdgeColor',coloR{i});
    end

    for r = 2:4
        nexttile(r)
        for i = 1:3
            WX{end + 1}  = plot(x{1}(1:end - 1) + diff(x{1},1)/2,mm{i}(:,1),'LineWidth',1.25*z,'color',coloR{i});
        end
        WX{end + 1}  = plot(NaN,NaN,'LineWidth',0.25*z,'color',[0.00 0.00 0.00 0.5]);
        WX{end + 1}  = plot(x{1}(1:end - 1) + diff(x{1},1)/2,mm{4},'LineWidth',0.25*z,'color',[0.00 0.00 0.00 0.0125]);
        for i = 1:2
            WX{end + 1}  = fill([x{1}(1:end - 1) + diff(x{1},1)/2;flip(x{1}(1:end - 1) + diff(x{1},1)/2)],[mm{i}(:,2);flip(mm{i}(:,3))],coloR{i},'FaceAlpha',.15,'EdgeAlpha',0.00,'LineWidth',0.25*z,'LineStyle','-','EdgeColor',coloR{i});
        end        
    end
    nexttile(3)
    WX{end + 1} = legend(lEGeND,'Interpreter','latex','FontSize',9*z,'FontAngle','oblique','Location','southoutside','NumColumns',2,'Box','off');
    
    if ismember(j,find(sEL))
        A            = char("00" + find(find(sEL) == j));
        saveas(gcf,char(pATh + "ReSuLTs/figures/fig In " + A(end - 2:end) + " " + iNFo.country(j)),'png');
    else
        A            = char("00" + sum(~ismember((1:j),find(sEL))));
        saveas(gcf,char(pATh + "ReSuLTs/figures/fig Out " + A(end - 2:end) + " " + iNFo.country(j)),'png');
    end

    for r = 1:numel(WX)
        delete(WX{r})
    end
    clear WX
end