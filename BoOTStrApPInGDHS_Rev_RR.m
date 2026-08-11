clear
R                       = 0;
pATh                    = "/Users/lshjr3/Documents/DHS/OuTPuT/";
fileattrib(pATh + "ReSuLTs");
RESolUTioN              = 300;
if R > 0
    options = detectImportOptions(char(pATh + "lISt.csv"));
    for i = 1:numel(options.VariableTypes)
        if isequal(options.VariableTypes{i},'char')
            options.VariableTypes{i} = 'string';
        end
    end
    iNFo                = readtable(char(pATh + "lISt.csv"),options);
    iNFo.ReTRo          = max(iNFo.ReTRo,10);
    redundant           = extractBetween(iNFo.fILe,1,2) == "IA" & extractBetween(iNFo.fILe + "AA",5,6) == "42";
    redundant           = redundant | ismember(iNFo.fILe,["IA23";"IA42";"IA52";"IA7E";"SNG0"]);
    iNFo                = iNFo(~redundant,:);
   
    iNFo.R              = ones(size(iNFo,1),1)*R;
    worked_examples     = ["CD61","MW7A"];
    WE                  = find(ismember(iNFo.fILe,worked_examples));
    WE                  = iNFo(WE,:);
    WE.R(:)             = 5000;
    iNFo                = [iNFo;WE]
    iNFo.j              = (1:size(iNFo,1))';
    
    options = detectImportOptions(char(pATh + iNFo.fILe(7) + ".csv"));
    for i = 1:numel(options.VariableTypes)
        if isequal(options.VariableTypes{i},'char')
            options.VariableTypes{i} = 'string';
        end
        if isequal(options.VariableTypes{i},'datetime')
            options.VariableOptions(1,i).InputFormat = 'dd/MM/yyyy';
        end
    end

    xE{1}               = [[(0:1:13) (14:7:28)]/365.25,[(2:1:12),(15:3:24),(36:12:60)]/12]';
    xE{2}               = [string(0);string([[(1:1:13) (14:7:28)]';(2:1:11)';(12:3:24)';(36:12:60)']) + char([kron('d',ones(16,1));kron('m',ones(18,1))])];
    xE{3}               = char(kron('$\mathit{',ones(35,1))) + [string(0);string([[(1:1:13) (14:7:28)]';(2:1:11)';(12:3:24)'/12;(36:12:60)'/12]) + char([kron('d',ones(16,1));kron('m',ones(10,1));kron('y',ones(8,1))])] + char(kron('}$',ones(35,1)));
    xE{4}               = round(xE{1}*365.25*24);
    structure           = {'stratum','cluster'};

    for i = 1:size(iNFo,1)
        clc;
        string(i) + ". " + iNFo.fILe(i)
        data              = readtable(char(pATh + iNFo.fILe(i) + ".csv"),options);
        R                 = iNFo.R(i);
        if isequal(std(day(data.interview + 1)),0)
            rng(0);
            date           = data.interview(data.k == 1);
            date           = datetime(year(date),month(date),1);
            p              = rand(size(date,1),1);
            date           = date + ceil(datenum(datetime(year(date),month(date) + 1,1) - date).*p);
            data.interview = repelem(date,data.K(data.k == 1));
            clear date p
        end        
        mAx               = median(data.interview);
        mIn               = datetime(year(mAx) - max(iNFo.ReTRo(i),10),month(mAx),day(mAx));
        T                 = [mIn mAx];
        mIn               = datetime(year(T(1)) - ceil(xE{1}(end)),month(T(1)),day(T(1)));
        mIn               = min(mIn,min(data.B_min));
        T                 = hours(T - mIn);

        rng(0);
        W                 = ones(size(data.W));
        for j = 1:numel(structure) - 1
            S   = groupsummary(data(:,structure(j:j + 1)),structure(j:j + 1));
            N   = groupcounts(S{:,structure(j)});
            if std(S{:,structure(j)}) > 0
                S.k = (1:sum(N))' - repelem(cumsum(N) - N,N);
                S.K = repelem(N,N);
                F   = [(1:sum(N))',ceil(rand(sum(N),R).*S.K) + repelem(cumsum(N) - N,N)];
            else
                S.k = (1:sum(N))';
                S.K = ones(N,1)*N;
                F   = [(1:sum(N))',ceil(rand(sum(N),R).*S.K)];
            end

            for h = 1:R + 1
                w{j}(:,h) = repelem(groupcounts([(1:sum(N))';F(:,h)]) - 1,S.GroupCount);
            end
            W   = W.*w{j};
            clear F S N h
        end
        w                 = W.*data.W;
        
        sEL               = (data.mother == 1);
        p                 = rand(sum(sEL),R + 1);
        p                 = round(hours((data.B_max(sEL) - data.B_min(sEL)).*p));
        B                 = hours(data.B_min(sEL) - mIn) + p;

        p                 = rand(sum(sEL),R + 1);
        d                 = round((p.*data.D_min(sEL) + (1 - p).*(data.D_max(sEL) - eps))*24);
        D                 = B + d;
        data.date         = hours(data.interview - mIn);
        O                 = min(D,data.date(sEL));
        
        sEX               = [data.sex(sEL) == [2 1] ones(sum(sEL),1)];
        for k = 1:size(sEX,2)
            wS{k}             = w(sEL,:).*sEX(:,k);
        end

        for j = 1:numel(xE{4}) - 1
            a              = max(min(B + xE{4}(j),min(O,T(2))),T(1));
            o              = max(min(B + xE{4}(j + 1),min(O,T(2))),T(1));
            eX             = (o - a)/(24*365.25);
            eV             = (d >= xE{4}(j) & d < xE{4}(j + 1) & D >= T(1) & D < T(2));
            for k = 1:size(sEX,2)
                exposure{k}(j,:)  = sum(eX.*wS{k});
                events{k}(j,:)    = sum(eV.*wS{k});
            end
            j
            clear a o eX eV
        end
        
        for k = 1:size(sEX,2)
            nMx               = events{k}./exposure{k};
            q                 = 1 - [ones(1,R + 1);exp(-cumsum(diff(xE{1},1).*nMx,1))];
            LT{k}{i,1}        = nMx;
            LT{k}{i,2}        = q;
            LT{k}{i,3}        = events{k};
            LT{k}{i,4}        = exposure{k};
            clear nMx q
        end
        iNFo.small(i)     = (min(LT{1}{i,2}(xE{2} == "7d",1),LT{2}{i,2}(xE{2} == "7d",1)) == 0);
        clear data d O D B mIn mAx w W wS j k p dATe T sEX exposure events sEL
    end
    exclusion                   = ["AF71";"AL71";"AM72";"BJ61";"DR41";"ES01";"GM61";"IA74";"MX01";"SL51";"SZ51";"TL71";"OS01";"MZ31";"ZW52";"KM61";"SNG0"];
    sET                         = ismember(iNFo.fILe,exclusion);
    iNFo.exclusion(sET)         = 1;

    N                           = 5000;
    sMAll                       = ["IA7E-R12";"IA7E-R13";"IA7E-R15";"IA7E-R17"];
    sET                         = ismember(iNFo.fILe,sMAll) | iNFo.N < N;
    iNFo.small(sET)             = 1;

    sEL                         = ismember(iNFo.SubRegion,["South Asia";"Western Africa";"Middle Africa";"Eastern Africa"]) & iNFo.small == 0 & ~ismember(iNFo.country,["Maldives";"Pakistan"]);
    iNFo.R(sEL)                 = "1. Model B";
    iNFo.R(~sEL)                = "2. Rest of the World";
    iNFo.R(iNFo.small == 1)     = "3. Small sample";
    iNFo.R(iNFo.exclusion == 1) = "4. Excluded";
    clear ans sET sEL N sMAll i exclusion redundant options WE
    
    for i = 1:size(LT,2)
        LTwe{i}           = LT{i}(end - numel(worked_examples) + 1:end,:);
        LT{i}             = LT{i}(1:end - numel(worked_examples),:);
    end
    iNFo                        = iNFo(1:end - numel(worked_examples),:);
    R                           = min(iNFo.R);
    save(char(pATh + "ReSuLTs/BoOTStrAp_DHS.mat"));
else
    load(char(pATh + "ReSuLTs/BoOTStrAp_DHS.mat"),'iNFo','LT','LTwe','worked_examples','pATh','R','structure','xE');
end

Z                           = 2;
order                       = Z + 2;
ReG                         = table({'q'},{'0'},{'60m'},{'q'},Z);
[BolsA,q,~,x]               = Coef(ReG);
Uf                          = BolsA{3}{1}(1:end - 1,2 + Z:end);
for j = 1:numel(q)
    B                              = Coef_BaYeS(q{j}',x,ReG,[],Uf,14);
    betaA{j}                       = B(:,1:order);
    for i = 1:order
        BetaA{j}{i} = prctile(betaA{j}{i},[50 2.5 97.5],2);
    end
    BetaA{j}                       = cell2mat(BetaA{j});
end

sET                         = find(ismember(xE{2},x{2}));
I                           = zeros(sET(end - 1),numel(sET) - 1);
for i = 2:numel(sET)
    I(sET(i - 1):sET(i) - 1,i - 1) = 1;
end

for i = 2:numel(x{2})
    leGend{i - 1}                  = {char("$\mathit{" + x{2}(i) + "}$")};
end


for i = 1:order
    if i <= Z + 1
        vARs{1}{i}  = char("$\mathrm{" + string(char(96 + i)) + "}\mathit{_x}$");
    else
        vARs{1}{i}  = char("$\mathrm{" + string(char(117 + (i - Z - 1))) + "}\mathit{_x}$");
    end
    foRMaT{i}   = '%0.4f';
end


ages                        = numel(x{1}) - 1;
F                           = table(cellstr(char(kron(ones(ages,1),'q'))),cellstr(char(kron(ones(ages,1),char(x{2}(1))))),cellstr(char(x{2}(2:end))),diff(x{1},1)/x{1}(end));
G                           = ReG(1,1:3);

sEL                         = (iNFo.R == "1. Model B" & iNFo.exclusion == 0);
ages                        = numel(xE{1}) - 1;
tables                      = numel(sEL);
p                           = [[4 10]/365.25;[2 4]/12;[5 7]/12;[8 10]/12;[11 24]/12];
P                           = {find(ismember(xE{2},["0";"1d";"2d";"3d";"7d";"14d";"28d";"3m";"6m";"12m";"24m";"60m"])),find(ismember(xE{2},["0";"6m";"12m";"18m";"24m";"36m";"48m";"60m"]))};
sR                          = find(ismember(xE{2},["0";"28d";"12m";"24m";"60m"]));
SmoothExample               = [13 18 26 27 54 278];
for j = 1:numel(LT)
    nMx                         = cell2mat(LT{j}(:,1))';
    nMxP                        = cell2mat(mat2cell(prctile(nMx(2:end,:),50)',ones(tables,1)*ages)');
    nMx                         = cell2mat(LT{j}(sEL,1));
    nMx                         = cell2mat(mat2cell(nMx(:,2:R + 1),ones(sum(sEL),1)*ages,R)');
    nMx                         = [nMxP nMx];
    q                           = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(xE{1},1).*nMx))];
    sET                         = q(end,:) > 0;

    if isequal(j,3)
        for i = 1:numel(SmoothExample)
            smoothingLT(q(:,SmoothExample(i)),xE,sR,100,10,0.25,p,P);
            exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/smoothing example " + iNFo.country(SmoothExample(j)) + " " + iNFo.fILe(SmoothExample(j)) + ".png"),'Resolution',RESolUTioN);    
        end
    end
    
    [qS,nMxS]                   = smoothingLT(q(:,sET),xE,sR,100,10,0.25,p,P);
    ndx                         = diff(qS,1);
    nLx                         = min(ndx./nMxS,diff(xE{1},1).*(1 - qS(1:end - 1,:)));
    nMxS                        = (I'*ndx)./(I'*nLx);
    qS                          = NaN(size(nMxS,1) + 1,numel(sET));
    qS(:,sET)                   = 1 - [ones(1,size(nMxS,2));exp(-cumsum(diff(x{1},1).*nMxS))];
    Qdhs_S{j}{1}                = qS(:,1:tables);
    Qdhs_S{j}{2}                = qS(:,tables + 1:end);

    ndx                         = diff(q,1);
    nLx                         = min(ndx./nMx,diff(xE{1},1).*(1 - q(1:end - 1,:)));
    nMx                         = (I'*ndx)./(I'*nLx);
    q                           = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(x{1},1).*nMx))];
    Qdhs_R{j}{1}                = q(:,1:tables);
    Qdhs_R{j}{2}                = q(:,tables + 1:end);
    
    B                           = Coef_BaYeS(Qdhs_S{j}{2},x,ReG,sum(sEL),Uf,6);
    betaB{j}                    = B(:,1:order);
    BolsB{j}                    = {Coef_fast(Qdhs_S{j}{2},x,ReG),ReG};
    for i = 1:order
        BetaB{j}{i} = prctile(betaB{j}{i},[50 2.5 97.5],2);
    end
    BetaB{j}                    = cell2mat(BetaB{j});
    X                           = [log(12.5/1000).^(0:Z) -.7];
    nMx                         = -diff(log(1 - [zeros(1,2);exp([BetaB{j}(:,1:3:end)*X' BolsB{j}{1}(:,1:order)*X'])]),1)./diff(x{1},1);
    hold on
    plot(x{1}(1:end - 1) + diff(x{1})/2,log(nMx))
    hold off
    
    q                           = Qdhs_R{j}{1};
    init                        = Set(q,x,{F,ReG});
    init                        = [log(init{2}),zeros(size(init{2}))];
    match                       = Set(q,x,{F,G});
    xo                          = Match(match,{BetaB{j}(:,1:3:end),ReG},x,{F,G},init,250);
    OuT{j}                      = [exp(xo(:,1)),xo(:,2)];
    clear q nMx qS nMxS ndx nLx sET B
end
if isequal(Z,2)
    save(char(pATh + "ReSuLTs/BoOTStrAp_DHS.mat"),'iNFo','LT','LTwe','worked_examples','pATh','R','RESolUTioN','structure','xE','x','sR','p','P','I','ReG','BetaA','betaA','BetaB','betaB','Qdhs_S','Qdhs_R','OuT','vARs','foRMaT','-v7.3','-nocompression')
else
    BetaA_c                     = BetaA;
    BetaB_c                     = BetaB;
    betaA_c                     = betaA;
    betaB_c                     = betaB;
    ReG_c                       = ReG;
    vARs_c                      = vARs;
    foRMaT_c                    = foRMaT;
    save(char(pATh + "ReSuLTs/Model_B" + Z + ".mat"),'x','ReG_c','BetaA_c','betaA_c','BetaB_c','betaB_c','vARs_c','foRMaT_c','-v7.3','-nocompression')
end


Region                      = tabulate(iNFo.R);
Region                      = sortrows(string(Region(:,1)));
for i = 1:numel(Region)
    sET               = find(iNFo.R == Region(i));
    info{i}           = iNFo(sET,{'fILe','survey','country','sTArt','eNd','DHS','Region','SubRegion'});
    for j = 1:numel(LT)
        lIFeTaBLeS{i}{j} = LT{j}(sET,2);
    end
        
    sET               = find(info{i}.sTArt ~= info{i}.eNd);
    yeAR              = string(info{i}.sTArt);
    yeAR(sET)         = yeAR(sET) + "-" + string(info{i}.eNd(sET));
    country           = info{i}.country(1);
    years             = yeAR(1);
    surveys           = info{i}.survey(1);
    files             = info{i}.fILe(1);
    tables            = 1;
    tABlE             = table(country,years,tables,surveys,files);
    clear country years tables surveys
    k                 = 1;
    K                 = k;
    for j = 2:size(info{i},1)
        if isequal(info{i}.country(j),tABlE.country(K))
            if ~isequal(mod(tABlE.tables(K),7),0)
                tABlE.years(k)   = tABlE.years(k) + ", " + yeAR(j);
                tABlE.surveys(k) = tABlE.surveys(k) + ", " + info{i}.survey(j);
                tABlE.files(k)   = tABlE.files(k) + ", " + info{i}.fILe(j);
                tABlE.tables(K)  = tABlE.tables(K) + 1;
            else
                k                = k + 1;
                tABlE(k,:)       = table("",yeAR(j),NaN,info{i}.survey(j),info{i}.fILe(j));
                tABlE.tables(K)  = tABlE.tables(K) + 1;
            end
        else
            k                = k + 1;
            K                = k;
            tABlE(k,:)       = table(info{i}.country(j),yeAR(j),1,info{i}.survey(j),info{i}.fILe(j));
        end
    end
    tABlE(k + 1,:)    = table("Total","",sum(recode(tABlE.tables,NaN,0)),"","");
    lABs              = {num2cell(1:k) {k + 1}};
    tABleSuMm({},{{'Surveys','Life Tables'}},{'left','right';10,11},lABs,{'$\mathrm{Country}$','$\mathit{}$'},tABlE.country,[tABlE.files string(tABlE.tables)],[.145 .380 .0])
    exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/Table (" + Region(i) + ").png"),'Resolution',RESolUTioN);
    taBleTexT{i}      = table(tABlE.country,tABlE.files,tABlE.tables);
    sEt{i}            = char(Region(i));
end

lABs                        = {num2cell(1:size(x{2},1) - 1)};
nOTe                        = {'$\mathrm{ages/}\mathit{coeff.}$','$\textrm{MCMC estimation via Gibbs sampling with 100,000 iterations after a warm-up period of 1000 iterations.}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$ $\textrm{conditional posterior distribution.}$'};
tABleBAyEs({'$\textrm{1. Model A}$','$\textrm{1. Model B}$'},[vARs vARs],foRMaT,lABs,nOTe,leGend,[BetaA{3} BetaB{3}],0.045,0.075,[]);
exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/Table Coefficients A & B.png"),'Resolution',RESolUTioN);

tABleBAyEs({'$\textrm{1. Female}$','$\textrm{2. Male}$','$\textrm{3. Both sexes}$'},[vARs vARs vARs],foRMaT,lABs,nOTe,leGend,[BetaB{1} BetaB{2} BetaB{3}],0.045,0.075,[]);
exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/Table Coefficients B.png"),'Resolution',RESolUTioN);

load(char(pATh + "ReSuLTs/Model_B3.mat"),'BetaA_c','BetaB_c','vARs_c','foRMaT_c')
tABleBAyEs({'$\textrm{1. Model B (log-quadratic, both sexes)}$','$\textrm{1. Model B (log-cubic, both sexes)}$'},[vARs vARs_c],[foRMaT foRMaT_c],lABs,nOTe,leGend,[BetaB{3} BetaB_c{3}],0.045,0.075,[]);
exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/Table Coefficients quadratic vs cubic, B.png"),'Resolution',RESolUTioN);

tABleBAyEs({'$\textrm{1. Model A (log-quadratic, both sexes)}$','$\textrm{1. Model A (log-cubic, both sexes)}$'},[vARs vARs_c],[foRMaT foRMaT_c],lABs,nOTe,leGend,[BetaA{3} BetaA_c{3}],0.045,0.075,[]);
exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/Table Coefficients quadratic vs cubic, A.png"),'Resolution',RESolUTioN);





clearvars -except pATh
load(char(pATh + "ReSuLTs/BoOTStrAp_DHS.mat"),'iNFo','LTwe','worked_examples','RESolUTioN','x','xE','I','ReG','BetaB','betaB','Qdhs_R')

rng(101);
sex                     = 3;
Z                       = ReG.Z;
W                       = diag(diff(x{1},1)/x{1}(end));
N                       = numel(x{1}) - 1;
warmup                  = 5000;
iterations              = 100000;
for j = 1:numel(worked_examples)
    s                = find(ismember(iNFo.fILe,worked_examples(j)));
    q                = Qdhs_R{sex}{1}(:,s);

    coefficients     = BetaB{sex}(:,1:3:end);
    V                = coefficients(:,end);
    X                = coefficients(:,1:end - 1)*(log(q(end,:)).^(0:Z))';
    e                = log(q(2:end)) - X;
    
    VV               = V'*W*V;
    Ve               = V'*W*e;
    ee               = e'*W*e;
    theta            = VV\Ve;
    ex1{j}.OLS       = theta;
    ex1{j}.VAR       = ((V'*W*V)\(e'*W*e) - theta^2)*N/(N - 1);
    ex1{j}.sd        = sqrt(ex1{j}.VAR);
    ex1{j}.se        = sqrt(ex1{j}.VAR)/sqrt(N);
    ex1{j}.sE        = sqrt((VV\ee - theta^2)/(N - 1));

    priork           = {0;10^-4};
    priorPhi         = {10^-2,10^-2};
    k                = 0;
    Phi              = 10^-2;
    for r = 1:iterations + warmup
        EE       = ee - 2*Ve*k + k*VV*k;
        Phi      = gamrnd(priorPhi{1} + N/2,1/(priorPhi{2} + 1/2*EE));
        
        s2       = 1./(Phi*VV + priork{2});
        Mu       = s2*(Phi*VV*theta + priork{2}*priork{1});
        k        = normrnd(Mu,sqrt(s2));
        
        if r > warmup
            ex1{j}.k(r - warmup,:)    = k;
            ex1{j}.Phi(r - warmup,:)  = Phi;
            ex1{j}.MSE(r - warmup,:)  = k.^2 + 1./Phi;
        end
        clear EE s2 Mu
    end
    k                = ex1{1}.k';
    q                = [zeros(1,iterations);exp(X + V*k)];
    ex1{j}.qBa       = prctile(q,[50 2.5 97.5],2);

    k                = ex1{j}.OLS + 1.96*ex1{j}.se*[0 -1 1];
    ex1{j}.qLse      = [zeros(1,3);exp(X + V*k)];
    k                = ex1{j}.OLS + 1.96*ex1{j}.sd*[0 -1 1];
    ex1{j}.qLsd      = [zeros(1,3);exp(X + V*k)];
    clear s q V e VV Ve ee theta priork priorPhi k Phi coefficients
end


warmup                  = 1000;
iterations              = 5000;
samples                 = 1000;
sEL                     = unidrnd(size(betaB{3}{1},2),samples,1);
W                       = diff(x{1},1)/x{1}(end);
for j = 1:numel(worked_examples)
    s                = find(ismember(iNFo.fILe,worked_examples(j)));
    q                = Qdhs_R{3}{1}(:,s);

    for i = 1:numel(sEL)
        X(:,i)   = [betaB{sex}{1}(:,sEL(i)) betaB{sex}{2}(:,sEL(i)) betaB{sex}{3}(:,sEL(i))]*(log(q(end)).^(0:Z))';
    end
    e                = log(q(2:end)) - X;
    V                = betaB{sex}{4}(:,sEL);

    VV               = sum(V.*W.*V)';
    Ve               = sum(V.*W.*e)';
    ee               = sum(e.*W.*e)';
    theta            = Ve./VV;
    ex2{j}.OLS       = theta;
    ex2{j}.VAR       = (sum(e.*W.*e)'./sum(V.*W.*V)' - theta.^2)*N/(N - 1);
    ex2{j}.sd        = sqrt(ex2{j}.VAR);
    ex2{j}.se        = sqrt(ex2{j}.VAR)/sqrt(N);
    ex2{j}.sE        = sqrt((ee./VV - theta.^2)/(N - 1));

    priork           = {zeros(samples,1);ones(samples,1)*10^-4};
    priorPhi         = {ones(samples,1)*10^-2,ones(samples,1)*10^-2};
    k                = zeros(samples,1);
    Phi              = ones(samples,1)*10^-2;    
    for r = 1:iterations + warmup
        EE       = ee - 2*Ve.*k + k.*VV.*k;
        Phi      = gamrnd(priorPhi{1} + N/2,1./(priorPhi{2} + 1/2*EE));
        
        s2       = 1./(Phi.*VV + priork{2});
        Mu       = s2.*(Phi.*VV.*theta + priork{2}.*priork{1});
        k        = normrnd(Mu,sqrt(s2));
        
        if r > warmup
            ex2{j}.k(:,r - warmup)    = k;
            ex2{j}.Phi(:,r - warmup)  = Phi;
            ex2{j}.MSE(:,r - warmup)  = k.^2 + 1./Phi;
        end
        clear EE s2 Mu
    end

    k                = ex2{j}.k';
    q                = {};
    for i = 1:iterations
        q{i}     = [zeros(1,samples);exp(X + V.*k(i,:))];
    end
    ex2{j}.qBa       = prctile(cell2mat(q),[50 2.5 97.5],2);
    q                = [zeros(1,samples);exp(X + V.*theta')];
    ex2{j}.qLC       = prctile(q,[50 2.5 97.5],2);
    clear s q X V e VV Ve ee theta priork priorPhi k Phi
end



for j = 1:numel(worked_examples)
    nMx              = LTwe{3}{j,1};
    q                = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(xE{1},1).*nMx))];
    ndx              = diff(q,1);
    nLx              = min(ndx./nMx,diff(xE{1},1).*(1 - q(1:end - 1,:)));
    nMx              = (I'*ndx)./(I'*nLx);
    q                = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(x{1},1).*nMx))];
    WE{j}            = q(:,2:end);
    T                = table(q','VariableNames',{'age_index'});
    writetable(T, char(pATh + "ReSuLTs/BoOTStrAp_WE_" + worked_examples(j) + ".csv"));
    clear q nMx T
end

warmup                  = 1000;
iterations              = 5000;
coefficients            = BetaB{sex}(:,1:3:end);
for j = 1:numel(worked_examples)
    q                = WE{j};
    V                = coefficients(:,end);
    X                = coefficients(:,1:end - 1)*(log(q(end,:)').^(0:Z))';
    e                = log(q(2:end,:)) - X;

    VV               = sum(V.*W.*V)';
    Ve               = sum(V.*W.*e)';
    ee               = sum(e.*W.*e)';
    theta            = Ve./VV;
    ex3{j}.OLS       = theta;
    ex3{j}.VAR       = (sum(e.*W.*e)'./sum(V.*W.*V)' - theta.^2)*N/(N - 1);
    ex3{j}.sd        = sqrt(ex3{j}.VAR);
    ex3{j}.se        = sqrt(ex3{j}.VAR)/sqrt(N);
    ex3{j}.sE        = sqrt((ee./VV - theta.^2)/(N - 1));
    
    samples          = numel(theta);
    priork           = {zeros(samples,1);ones(samples,1)*10^-4};
    priorPhi         = {ones(samples,1)*10^-2,ones(samples,1)*10^-2};
    k                = zeros(samples,1);
    Phi              = ones(samples,1)*10^-2;    
    for r = 1:iterations + warmup
        EE       = ee - 2*Ve.*k + k.*VV.*k;
        Phi      = gamrnd(priorPhi{1} + N/2,1./(priorPhi{2} + 1/2*EE));
        
        s2       = 1./(Phi.*VV + priork{2});
        Mu       = s2.*(Phi.*VV.*theta + priork{2}.*priork{1});
        k        = normrnd(Mu,sqrt(s2));
        
        if r > warmup
            ex3{j}.k(:,r - warmup)    = k;
            ex3{j}.Phi(:,r - warmup)  = Phi;
            ex3{j}.MSE(:,r - warmup)  = k.^2 + 1./Phi;
        end
        clear EE s2 Mu
    end

    k                = ex3{j}.k';
    q                = {i};
    for i = 1:iterations
        q{i}     = [zeros(1,samples);exp(X + V.*k(i,:))];
    end
    ex3{j}.qBa       = prctile(cell2mat(q),[50 2.5 97.5],2);
    q                = [zeros(1,samples);exp(X + V.*theta')];
    ex3{j}.qLS       = prctile(q,[50 2.5 97.5],2);
    clear s q X V e VV Ve ee theta priork priorPhi k Phi samples
end


options                     = detectImportOptions(char(pATh + "ReSuLTs/MW7A_turnbull.csv"));
turnbull                    = readtable(char(pATh + "ReSuLTs/MW7A_turnbull.csv"),options);
turnbull                    = {reshape(turnbull.qx,23,[])};
for j = 1:numel(turnbull)
    q                = turnbull{j};
    V                = coefficients(:,end);
    X                = coefficients(:,1:end - 1)*(log(q(end,:)').^(0:Z))';
    e                = log(q(2:end,:)) - X;

    VV               = sum(V.*W.*V)';
    Ve               = sum(V.*W.*e)';
    ee               = sum(e.*W.*e)';
    theta            = Ve./VV;
    ex3T{j}.OLS      = theta;
    ex3T{j}.VAR      = (sum(e.*W.*e)'./sum(V.*W.*V)' - theta.^2)*N/(N - 1);
    ex3T{j}.sd       = sqrt(ex3{j}.VAR);
    ex3T{j}.se       = sqrt(ex3{j}.VAR)/sqrt(N);
    ex3T{j}.sE       = sqrt((ee./VV - theta.^2)/(N - 1));
    
    samples          = numel(theta);
    priork           = {zeros(samples,1);ones(samples,1)*10^-4};
    priorPhi         = {ones(samples,1)*10^-2,ones(samples,1)*10^-2};
    k                = zeros(samples,1);
    Phi              = ones(samples,1)*10^-2;    
    for r = 1:iterations + warmup
        EE       = ee - 2*Ve.*k + k.*VV.*k;
        Phi      = gamrnd(priorPhi{1} + N/2,1./(priorPhi{2} + 1/2*EE));
        
        s2       = 1./(Phi.*VV + priork{2});
        Mu       = s2.*(Phi.*VV.*theta + priork{2}.*priork{1});
        k        = normrnd(Mu,sqrt(s2));   %%%
        
        if r > warmup
            ex3T{j}.k(:,r - warmup)   = k;
            ex3T{j}.Phi(:,r - warmup) = Phi;
            ex3T{j}.MSE(:,r - warmup) = k.^2 + 1./Phi;
        end
        clear EE s2 Mu
    end
    
    k                = ex3T{j}.k';
    q                = {i};
    for i = 1:iterations
        q{i}     = [zeros(1,samples);exp(X + V.*k(i,:))];
    end
    ex3T{j}.qBa      = prctile(cell2mat(q),[50 2.5 97.5],2);
    q                = [zeros(1,samples);exp(X + V.*theta')];
    ex3T{j}.qLS      = prctile(q,[50 2.5 97.5],2);
    clear s q X V e VV Ve ee theta priork priorPhi k Phi samples
end


for j = 1:numel(worked_examples)
    T           = table(x{2},x{1},ex1{j}.qBa,ex1{j}.qLse,ex1{j}.qLsd,ex2{j}.qBa,ex2{j}.qLC,ex3{j}.qBa,ex3{j}.qLS,'VariableNames',{'age','x','Bayes Ex1: Fitting k','OLS(se)','OLS(sd)','Bayes Ex2: Fitting k + Coef.','OLS Coef. only','Bayes Ex3: Fitting k + BS survey','OLS BS survey only'});
    writetable(T, char(pATh + "ReSuLTs/Uncertainty_" + worked_examples(j) + ".csv"));
    SD(:,j)     = [std(ex1{j}.k);ex1{j}.se;ex1{j}.sd;std(ex2{j}.k(:));std(ex2{j}.OLS);std(ex3{j}.k(:));std(ex3{j}.OLS)];
end

T           = table(x{2},x{1},ex3T{1}.qBa,ex3T{1}.qLS,'VariableNames',{'age','x','Bayes Ex3: Fitting k + BS survey (Turnbull)','OLS BS survey only (Turnbull)'});
writetable(T, char(pATh + "ReSuLTs/Uncertainty_MW7A_Turnbull.csv"));
SDt         = [NaN(2,1) [std(ex3T{1}.k(:));std(ex3T{1}.OLS)]];
SD          = [SD;SDt];
examples    = {'Bayes Ex1: Fitting k','OLS(se)','OLS(sd)','Bayes Ex2: Fitting k + Coef.','OLS Coef. only','Bayes Ex3: Fitting k + BS survey','OLS BS survey only','Bayes Ex3: Fitting k + BS survey (Turnbull)','OLS BS survey only (Turnbull)'};
T           = table(examples',SD(:,1),SD(:,2),'VariableNames',{'example',char(worked_examples(1)),char(worked_examples(2))});
writetable(T, char(pATh + "ReSuLTs/Uncertainty_sd_of_k.csv"));