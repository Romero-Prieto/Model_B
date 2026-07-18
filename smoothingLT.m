function [qs,nMxS] = smoothingLT(q,x,s,R,r,lambda,p,P)


tables    = size(q,2);
y         = -log(1 - q);

if isequal(numel(lambda),0)
    lambda = zeros(size(q,2),0);
elseif isequal(numel(lambda),1)
    lambda = ones(size(q,2),1)*lambda;
end

sET       = sum(y == 0);
mIn       = union(sET,[]);
for i = 1:numel(mIn)
    rOWs              = mIn(i):numel(x{1});
    cOLs              = sET == mIn(i);
    S                 = (rOWs(1):(numel(rOWs) - 1)/(R - 1):rOWs(end))';
    X                 = interMonotonic(union(max(s,mIn(i)),mIn(i)),x{1}(union(max(s,mIn(i)),mIn(i))),S,[]);
    
    sEL               = isnan(X);
    for j = 1:size(p,1)
        sEL = sEL | (X > p(j,1) & X <= p(j,2));
    end
    xo                = X(~sEL);
    yi                = Monotonic(x{1}(rOWs),y(rOWs,cOLs),xo);
    K                 = xo > 0 & yi(:,1) > 0;
    [sM,lAMbda(cOLs)] = BSplineEM(log(xo(K)),r,3,3,lambda(cOLs),log(yi(K,:)));
    fIX               = find(any(sM(2:end,:) < sM(1:end - 1,:),1))
    sM(:,fIX)         = log(cumsum([exp(sM(1,fIX));max(diff(exp(sM(:,fIX)),1),0)],1));

    ySm               = yi;
    ySm(K,:)          = exp(sM);
    ySm               = Monotonic(xo,ySm,x{1}(rOWs));
    ys(rOWs,cOLs)     = ySm;
end


qs        = 1 - exp(-ys);
nMx       = diff(y,1)./diff(x{1},1);
nMxS      = diff(ys,1)./diff(x{1},1);

if isequal(tables,1)
    mPIX                     = 538756;
    pix                      = 1/37.7952755906;
    z                        = min(sqrt(mPIX/(20*30/pix^2)),1);
    fi                       = figure('Color',[1 1 1],'Position',z*[0 0 20 30]/pix,'Theme','light');
    axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
    hold(axes1,'on');
    TL                       = tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
    coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.05 0.05 0.05],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65]};
    title(TL,char("$\mathrm{" + r + "\,}\mathit{knots},\,\mathit{\lambda=\,}\mathrm{" + sprintf('%0.3f',lambda) +"}$"),'Interpreter','latex','FontSize',12*z);

    for i = 1:6
        nexttile(i)
        ax{i}                            = gca;
        ax{i}.FontName                   = 'Times New Roman';
        ax{i}.FontSize                   = 9*z;

        ax{i}.XAxis.MinorTickValues      = x{1};
        ax{i}.XAxis.TickLabelInterpreter = 'latex';
        ax{i}.XAxis.MinorTick            = 'off';
        ax{i}.YAxis.MinorTick            = 'off';
        ax{i}.XAxis.TickLabelFormat      = '%.1f';
        ax{i}.YAxis.TickLabelFormat      = '%.4f';
        ax{i}.YScale                     = 'log';
        ax{i}.XTickLabelRotation         = 0;
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.LabelFontSizeMultiplier    = 1;

        grid on;
        grid minor;
        box  on;
        hold on;
        if isequal(ceil(i/2),1)
            mIn                              = y(2)*0.95;
            mAx                              = y(end)*1.05;
            ax{i}.YAxis.TickValues           = .1./(2.^(7:-1:-6));
            ax{i}.YAxis.MinorTickValues      = .1./(2.^(7:-1:-6));

            if isequal(i,1)
                xlim([x{1}(2)*0.5 x{1}(end)])
                ax{i}.XScale                     = 'log';
                sET                              = P{1};
                ylabel('$\int_{0}^{x} m(y) \,dy\textit{   (log scale)}$','Interpreter','latex','FontSize',11*z);
            else
                xlim([-0.1 x{1}(end)])
                sET                              = P{2};                
            end
            plot(x{1},y,'-o','MarkerSize',3.0,'color',coloR{2},'LineWidth',0.5,'MarkerFaceColor',coloR{2});
            scatter(X(~sEL),yi,10,'filled','MarkerFaceColor',coloR{7},'MarkerFaceAlpha',1.0);
            plot(x{1},ys,'LineWidth',1.15,'color',coloR{3});
        elseif isequal(ceil(i/2),2)
            mIn                              = min([nMx(nMx > 0);nMxS])*0.95;
            mAx                              = max([nMx;nMxS])*1.05;
            ax{i}.YAxis.TickValues           = .1./(2.^(7:-1:-6));
            ax{i}.YAxis.MinorTickValues      = .1./(2.^(7:-1:-6));

            if isequal(i,3)
                xlim([x{1}(2)*0.5 x{1}(end)])
                ax{i}.XScale                     = 'log';
                sET                              = P{1};
                ylabel('$\mathit{_nM_x} \textit{   (log scale)}$','Interpreter','latex','FontSize',11*z);
            else
                xlim([-0.1 x{1}(end)])
                sET                              = P{2};
            end
            plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMx,'-o','MarkerSize',3.0,'color',coloR{2},'LineWidth',0.5,'MarkerFaceColor',coloR{2});
            plot(x{1}(1:end - 1) + diff(x{1},1)/2,nMxS,'LineWidth',1.15,'color',coloR{3});
        elseif isequal(ceil(i/2),3)
            mIn                              = 0;
            mAx                              = max([q;qs])*1.05;
            ax{i}.YScale                     = 'linear';

            if isequal(i,5)
                xlim([x{1}(2)*0.5 x{1}(end)])
                ax{i}.XScale                     = 'log';
                sET                              = P{1};
                ylabel('$\mathit{q}\mathrm{(}\mathit{x}\mathrm{)}\textit{   (linear scale)}$','Interpreter','latex','FontSize',11*z);
                xlabel('$\textit{age   (log scale)}$','Interpreter','latex','FontSize',11*z);
            else
                xlim([-0.1 x{1}(end)])
                sET                              = P{2};
                xlabel('$\textit{age}$','Interpreter','latex','FontSize',11*z);
            end
            plot(x{1},q,'-o','MarkerSize',3.0,'color',coloR{2},'LineWidth',0.5,'MarkerFaceColor',coloR{2});
            plot(x{1},qs,'LineWidth',1.15,'color',coloR{3});
            ax{i}.YAxis.MinorTickValues      = ax{i}.YAxis.TickValues;
        end
        ax{i}.XAxis.TickValues           = x{1}(sET);
        ax{i}.XAxis.TickLabels           = x{3}(sET);
        ax{i}.XAxis.MinorTickValues      = x{1}(sET);

        ylim([mIn mAx])
        for j = 1:size(p,1)            
            fill([p(j,1);p(j,2);p(j,2);p(j,1)],[mIn;mIn;mAx;mAx],'k','FaceAlpha',.15,'EdgeAlpha',0.00,'LineWidth',0.25,'LineStyle','-','EdgeColor','k');
        end      
    end
end
