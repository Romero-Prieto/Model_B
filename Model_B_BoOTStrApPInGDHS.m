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
    for i = 1:size(iNFo,1)
        clc;
        string(i) + ". " + iNFo.fILe(i)
        data              = readtable(char(pATh + iNFo.fILe(i) + ".csv"),options);
        if isequal(std(day(data.interview)),0)
            rng(0);
            date           = data.interview(data.k == 1);
            date           = datetime(year(date),month(date),1);
            p              = rand(size(date,1),1);
            date           = date + round(datenum(datetime(year(date),month(date) + 1,1) - date).*p,5);
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
        s                 = data(data.k == 1,{'cluster','K','woman','Women','iNDeX'});
        w                 = rand(size(s,1),R);
        w                 = [(1:size(s,1))',ceil(s.Women.*w) + s.iNDeX];
        W                 = NaN(size(data,1),R + 1);
        for r = 1:R + 1
            S      = tabulate([w(:,1);w(:,r)]);
            W(:,r) = repelem(S(:,2) - 1,s.K);
        end
        clear r s S w
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
        B                 = min(B,data.date(sEL));
        B                 = min(B,data.date(sEL) - d);
        
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

    clear ans sET sEL N sMAll i exclusion redundant options
    save(char(pATh + "ReSuLTs/BoOTStrAp_DHS.mat"));
else
    load(char(pATh + "ReSuLTs/BoOTStrAp_DHS.mat"));
end


Z                           = 2;
order                       = Z + 2;
ReG                         = table({'q'},{'0'},{'60m'},{'q'},Z);
[~,~,~,x]                   = Coef(ReG);
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
        vARs{1}{i}  = char("$\mathrm{" + string(char(111 + i + Z + 1)) + "}\mathit{_x}$");
    end
    foRMaT{i}   = '%0.4f';
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
    p                 = [[4 10]/365.25;[2 4]/12;[5 7]/12;[8 10]/12;[11 24]/12];
    P                 = {find(ismember(xE{2},["0";"1d";"2d";"3d";"7d";"14d";"28d";"3m";"6m";"12m";"24m";"60m"])),find(ismember(xE{2},["0";"6m";"12m";"18m";"24m";"36m";"48m";"60m"]))};
    sR                = find(ismember(xE{2},["0";"28d";"12m";"24m";"60m"]));

    if i <= 2
        sEL               = (iNFo.R == Region(i) & iNFo.exclusion == 0);
        for j = 1:numel(LT)
            nMx         = cell2mat(LT{j}(sEL,1))';
            nMx         = cell2mat(mat2cell(prctile(nMx(2:end,:),50)',ones(sum(sEL),1)*size(LT{j}{1,1},1))');
            q           = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(xE{1},1).*nMx))];
            [qS,nMxS]   = smoothingLT(q,xE,sR,100,10,0.25,p,P);
            ndx         = qS(2:end,:) - qS(1:end - 1,:);
            nLx         = min(ndx./nMxS,diff(xE{1},1).*(1 - qS(1:end - 1,:)));
            nMxS        = (I'*ndx)./(I'*nLx);
            %qS         = q(ismember(xE{2},x{2},:);
            qS          = 1 - [ones(1,size(nMxS,2));exp(-cumsum(diff(x{1},1).*nMxS))];
            B           = Coef_fast(qS,x,ReG);
            betaS{i}{j} = B(:,1:order);

            ndx         = q(2:end,:) - q(1:end - 1,:);
            nLx         = min(ndx./nMx,diff(xE{1},1).*(1 - q(1:end - 1,:)));
            nMx         = (I'*ndx)./(I'*nLx);
            %q          = q(ismember(xE{2},x{2},:);
            q           = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(x{1},1).*nMx))];
            B           = Coef_fast(qS,x,ReG);
            betaR{i}{j} = B(:,1:order);
            clear q nMx qS nMxS ndx nLx B 
        end
        sEt{i} = char(Region(i));
        clear sEL
    end
end

lABs                        = {num2cell(1:size(x{2},1) - 1)};
nOTe                        = {'$\mathrm{ages/}\mathit{coeff.}$',char("$\mathrm{OLS-estimation,}$ Bootstrapping each smoothed survey " + string(R) + " times with replacement")};
tABleX(sEt,[vARs vARs],foRMaT,lABs,nOTe,leGend,[betaS{1}{3} betaS{2}{3}],0.045,0.045,[])
exportgraphics(gcf,char(pATh + "ReSuLTs/T&F/Table Coefficients.png"),'Resolution',RESolUTioN);

sEL                         = (iNFo.R ~= '');
ages                        = size(x{1},1) - 1;
F                           = table(cellstr(char(kron(ones(ages,1),'q'))),cellstr(char(kron(ones(ages,1),char(x{2}(1))))),cellstr(char(x{2}(2:end))),diff(x{1},1)/x{1}(end));
G                           = ReG(1,1:3);

for j = 1:numel(LT)
    nMx        = cell2mat(LT{j}(sEL,1))';
    nMx        = cell2mat(mat2cell(prctile(nMx(2:end,:),50)',ones(sum(sEL),1)*size(LT{j}{1,1},1))');
    q          = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(xE{1},1).*nMx))];
    ndx        = q(2:end,:) - q(1:end - 1,:);
    nLx        = min(ndx./nMx,diff(xE{1},1).*(1 - q(1:end - 1,:)));
    nMx        = (I'*ndx)./(I'*nLx);
    q          = 1 - [ones(1,size(nMx,2));exp(-cumsum(diff(x{1}).*nMx,1))];
    %q          = q(ismember(xE{2},x{2},:);
    init       = Set(q,x,{F,ReG});
    init       = [log(init{2}),zeros(size(init{2}))];
    match      = Set(q,x,{F,G});
    xo         = Match(match,{betaS{1}{j},ReG},x,{F,G},init,250);
    OuT{j}     = [exp(xo(:,1)),xo(:,2)];
    Qdhs{j}    = q;
end