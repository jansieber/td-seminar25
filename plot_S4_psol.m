%% Plot periodic orbits of coupled-oscillator demo
ddebiftool_path([pwd(),'/ddebiftool_25_03_07']);
format compact
format short g
load('S4_demo_psol_results.mat');
%%
txt={'fontsize',16};
ltx={'interpreter','latex','fontsize',20};
lw={'linewidth',2};
lww={'linewidth',4};
clr=lines();
gettauc=@(b)arrayfun(@(p)p.parameter(ip.tau_c),b.point);
gettaus=@(b)arrayfun(@(p)p.parameter(ip.tau_s),b.point);
gettauci=@(b,rg)arrayfun(@(p)p.parameter(ip.tau_c),b.point(rg));
gettausi=@(b,rg)arrayfun(@(p)p.parameter(ip.tau_s),b.point(rg));
ym=@(p)[2,0,0,0]*dde_coll_int(p,[0,1/2]);
getym=@(b)arrayfun(ym,b.point);
%% Plot 2\int_0^{1/2} x_1(s) ds 
pn=fieldnames(psolref);
figure(1);clf;tl=tiledlayout(4,3,"TileSpacing","tight");ax1=nexttile([4,2]);hold(ax1,'on');set(ax1,txt{:},lw{:},'box','on');
xlabel(ax1,'$\tau_\mathrm{c}$',ltx{:});
yl=ylabel(ax1,'$\int_0^{1/2}x_1(s)\mathrm{d}\,s$',ltx{:});
yl.Position([1,2])=[-0.2,0.75];
plargs={'max_nunst',2,'ax',ax1,'funcs',funcs,lww{:}};
Plot2dBranch(eq_1234ref,'y',@(p)p.x(1),'DisplayName','eq$(1234)$',plargs{:},'color',clr(1,:));
for i=1:length(pn)
    Plot2dBranch(psolref.(pn{i}),'y',ym,'DisplayName',psymrealname.(pn{i}),plargs{:},'color',clr(i+1,:));
end
hopf1=eq_1234ref.point(eq_1234bifind);
plot(ax1,hopf1.parameter(ip.tau_c),hopf1.x(1),'ko',lw{:},'DisplayName','equiv. Hopf (3)','MarkerSize',8);
sb_psol=psolref.p1234_1_4.point(psolbifind.p1234_1_4(1));
plot(ax1,sb_psol.parameter(ip.tau_c),ym(sb_psol),'ks','DisplayName','SB$\to(13)_{1/2}(24)_{1/2}$',...
    'MarkerSize',9,'MarkerFaceColor','k');
lg=legend(ax1,'Location','WestOutside');
lg.String=strrep(strrep(lg.String,'#',''),'>=','$\geq$');
lg.Interpreter='latex';
lg.FontSize=20;
lg.Box='off';
ax1.XLim=[0,4.5];
clp=[1,2,3,4];
lwp=[3,3,2,2];
lsp={'-','--','-',':'};
locp=[1,4,2,3];
for i=length(pn):-1:1
    iplot(i)=find(gettauc(psolref.(pn{i}))<2,1,'last');
    ppt(i)=psolref.(pn{i}).point(iplot(i));
    plot(ax1,ppt(i).parameter(ip.tau_c),ym(ppt(i)),'k^','MarkerSize',10,...
        'MarkerFaceColor',clr(i+1,:),'HandleVisibility','off',lw{:});
    axp(i)=nexttile(tl,3*locp(i));hold(axp(i),'on');
    plot(axp(i),NaN,NaN,'s','Color',clr(i+1,:),'DisplayName',psymrealname.(pn{i}),'MarkerFaceColor',clr(i+1,:));
    for k=1:4
        plot(axp(i),ppt(i).mesh,ppt(i).profile(k,:),'LineWidth',lwp(k),'Color',...
            clr(clp(k),:),'LineStyle',lsp{k},'DisplayName',sprintf('$x_%d$',k));
    end
    lgp(i)=legend(axp(i));
    lgp(i).Location='EastOutside';
    lgp(i).Interpreter='latex';
    %lgp(i).Title.String=psymrealname.(pn{i});
    lgp(i).Box='off';
    set(axp(i),txt{:},'box','on');
    axp(i).XTick=[];
    axp(i).YTick=[];
    axp(i).YLim=3*[-1,1];
end
%%
exportgraphics(figure(1),'../pics/coupled_oscbif1.pdf');
%% Plot two-parameter partial bifurcation diagram in (tau_s,tau_c) plane
figure(2);clf;ax2=gca;hold(ax2,'on');set(ax2,txt{:},lw{:},'box','on');
xlabel(ax2,'$\tau_\mathrm{c}$',ltx{:});
yl=ylabel(ax2,'$\tau_\mathrm{s}$',ltx{:});
yl.Rotation=0;
%yl.Position([1,2])=[-0.2,0.75];
fill(ax2,gettauc(pfbrref),gettaus(pfbrref),(clr(3,:)+1)/2,'HandleVisibility','off');
plot(ax2,gettauc(ehopfref),gettaus(ehopfref),lw{:},'Color',clr(1,:),'DisplayName','equiv. Hopf (3)');
plot(ax2,gettauc(shopfref),gettaus(shopfref),lw{:},'Color',clr(6,:),'DisplayName','sym. Hopf');
plot(ax2,gettauc(pfbrref),gettaus(pfbrref),lw{:},'Color',clr(2,:),'DisplayName','SB (per. orb.)');
plot(ax2,gettauc(trbr),gettaus(trbr),lw{:},'Color',clr(4,:),'DisplayName','equiv. torus bif. (3)');
pt_attr={'Color','k',lw{:},'MarkerSize',10};
plot(ax2,gettauci(shopfref,shopfbifind),gettausi(shopfref,shopfbifind),'s',pt_attr{:},'DisplayName','"double" Hopf (equiv.)','MarkerFaceColor',clr(6,:))
plot(ax2,gettauci(ehopfref,ehopfbifind),gettausi(ehopfref,ehopfbifind),'s',pt_attr{:},'DisplayName','"double" Hopf (equiv.)','MarkerFaceColor',clr(6,:),'HandleVisibility','off');
[~,imx]=max(gettaus(trbr));
plot(ax2,gettauci(trbr,imx),gettausi(trbr,imx),'o',pt_attr{:},'DisplayName','1:2 resonance (equiv.)','MarkerFaceColor',clr(3,:));
ax2.XLim=[0,2.5];
ax2.YLim(1)=0;
ax2.YTick=0:2;
ax2.XTick=0:4;
lg2=legend(ax2,'Location','North',ltx{:});
text(ax2,1.5,0.4,'p$(1234)_{1/4}$ stable',ltx{:})
text(ax2,1.4,0.8,'p$(12)_0(34)_0$ stable',ltx{:})
text(ax2,0.15,1.75,'p$(1234)_0$ stable','Rotation',90,ltx{:})
%%
exportgraphics(figure(2),'../pics/coupled_oscbif2.pdf');
