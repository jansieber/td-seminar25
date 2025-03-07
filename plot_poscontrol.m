%% plot for TDS seminar
clear
load('poscontrol_results.mat');
fn={'fontsize',18};
ltx={'interpreter','latex','fontsize',24};
lww={'linewidth',5};
lw={'linewidth',2};
clr=lines();
figure(1);clf;ax=gca;hold(ax,'on');
gettau=@(b)arrayfun(@(p)p.parameter(ip.tau0),b.point);
getk=@(b)arrayfun(@(p)p.parameter(ip.k),b.point);
gettaui=@(b,rg)arrayfun(@(p)p.parameter(ip.tau0),b.point(rg));
getki=@(b,rg)arrayfun(@(p)p.parameter(ip.k),b.point(rg));
ph=plot(ax,gettau(hopf_wbifs),getk(hopf_wbifs),'-',lww{:},'color',clr(1,:),....
    'DisplayName','Hopf');
pfp=plot(ax,gettau(POfold),getk(POfold),'-',lww{:},'color',clr(2,:),....
    'DisplayName','per.orb. fold');
ps1=plot(ax,gettau(grazing),getk(grazing),'-',lw{:},'color',clr(3,:),....
    'DisplayName','sonic barrier');
pgh=plot(ax,gettaui(hopf_wbifs,bif2ind),getki(hopf_wbifs,bif2ind),'ko',lw{:},'MarkerFaceColor',clr(2,:),....
    'MarkerSize',8,'DisplayName','gen. Hopf');
legend(ax,'Location','NorthEast',fn{:});
set(ax,fn{:},lw{:},'box','on',...
    'FontName','Courier','FontWeight','bold');
ax.YTick=0:2:6;
ax.XTick=0:2;
xl=xlabel(ax,'$\tau_0$',ltx{:});
xl.Position(1:2)=[1.5,-0.05];
yl=ylabel(ax,'$k$',ltx{:});
yl.Rotation=0;
yl.Position(1:2)=[-0.08,5];
yl.VerticalAlignment='middle';
%%
exportgraphics(figure(1),'../pics/poscontrolbif.pdf')
%%
figure(2);clf;ax2=gca;hold(ax2,'on');
plot(ax2,gettau(hopf_wbifs),hopftests.genh(1,:),lw{:});
yline(ax2,0,lw{:})
xl=xlabel(ax2,'$\tau_0$',ltx{:});
yl=ylabel(ax2,'$k$',ltx{:});
yl.Rotation=0;
yl.Position(1:2)=[-0.08,0.05];
yl.VerticalAlignment='middle';
ax2.YTick=-0.1:0.1:0.1;
set(ax2,fn{:},lw{:},'box','on',...
    'FontName','Courier','FontWeight','bold');
%%
exportgraphics(figure(2),'../pics/poscontrol_l1.pdf')
