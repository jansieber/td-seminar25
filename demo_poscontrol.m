%% Demos for TDS25 seminar
clear
ddebiftool_path([pwd(),'/ddebiftool_25_03_07']);
%% Define system and 1st derivative
f=@(x,s,xref,tau0,k,c)[...
    k.*xref-k.*s(2,:).*c/2;...      %delays: tau0, s(t)
    c.*s(1,:)-x(3,:)-x(1,:)];
df=@(x,s,xref,tau0,k,c,dx,ds,dxref,dtau0,dk,dc)[...
    dk.*xref+k.*dxref-dk.*s(2,:).*c/2-k.*ds(2,:).*c/2-k.*s(2,:).*dc/2;...
    dc.*s(1,:)+c.*ds(1,:)-dx(3,:)-dx(1,:)];
tau={@(x,s,xref,tau0,k,c)tau0, @(x,s,xref,tau0,k,c)s(1,:)};
dtau={@(x,s,xref,tau0,k,c,dx,ds,dxref,dtau0,dk,dc)dtau0,...
    @(x,s,xref,tau0,k,c,dx,ds,dxref,dtau0,dk,dc)ds(1,:)};
parnames={'xref','tau0','k','c'};
ip=struct_par(parnames{:});
fcn=set_funcs('sys_rhs',{1:2,[ip.xref,ip.tau0,ip.k,ip.c],f},...
    'sys_dirderi',{df},...
    'sys_ntau',@()2,'sys_tau',tau,'sys_dirdtau',{dtau},'lhs_matrix',diag([1,0]),...
    'x_vectorized',true,'p_vectorized',true);
%%  Define initial parameters and states
[parini([ip.tau0,ip.xref,ip.c,ip.k]), xini]=deal(...
    [       0,       1,     2,   6], [1;1]);
%% Continue the equilibrium branch (with online plotting)
eqs=SetupStst(fcn,'x',xini,'parameter',parini,'contpar',ip.tau0,'step',0.01,...
    'max_step',[ip.tau0,0.05;ip.k,0.1],'min_bound',[ip.tau0,0;ip.k,0],'max_bound',[ip.tau0,2;ip.k,6]);
figure(1);clf;ax1=gca;xlabel(ax1,'tau0');ylabel(ax1,'x (eq)');set(ax1,'FontSize',16);
eqs=br_contn(fcn,eqs,20,'ax',ax1);
[eqs_wbifs,~,indhopf,~]=LocateSpecialPoints(fcn,eqs);
%% Continue Hopf bifurcations in (tau0,k) plane
[hopf,suc]=SetupHopf(fcn,eqs_wbifs,indhopf(1),'contpar',[ip.tau0,ip.k],'dir',ip.k,'step',-0.1);
figure(2);clf;ax2=gca;xlabel(ax2,'tau0');ylabel(ax2,'k');set(ax2,'FontSize',16);
hopf=br_contn(fcn,hopf,100,'ax',ax2);
[hopf_wbifs,hopftests,bif2ind]=LocateSpecialPoints(fcn,hopf);
%% Branch off at generalized Hopf bifurcation
C1info=BranchFromCodimension2(fcn,hopf_wbifs,'step',1e-1,'max_step',[ip.tau0,0.05;ip.k,0.1],'min_bound',[ip.tau0,0;ip.k,0],'max_bound',[ip.tau0,2;ip.k,6]);
%%
POfold=br_contn(C1info(2).funcs,C1info(2).branch,50,'ax',ax2);
POfold=br_stabl(C1info(2).funcs,POfold);
%%
figure(3);clf;ax3=gca;hold(ax3,'on');xlabel(ax3,'tau0');ylabel(ax3,'k');set(ax3,'FontSize',16);
Plot2dBranch(hopf_wbifs,'ax',ax3);
Plot2dBranch(POfold,'ax',ax3,'funcs',C1info(2).funcs);
%% Branch off at Hopf bifurcation for periodic orbits
[per,suc]=SetupPsol(fcn,hopf_wbifs,60,'contpar',ip.tau0,'max_step',[0,0.1],...
    'degree',4,'intervals',20,'matrix','sparse','eigmatrix','sparse',...
    'plot_measure',{@(p)p.parameter(ip.tau0),@(p)max(p.profile(1,:))},'collocation_parameters','force_smooth');
figure(1);
per=br_contn(fcn,per,40,'plotaxis',ax1);
%% plot delays
figure(4);clf;ax4=gca;xlabel(ax4,'t (period)');ylabel(ax4,'delay s(t))');set(ax4,'FontSize',16);hold(ax4,'on');
for i=1:length(per.point)
    plot(ax4,per.point(i).mesh*per.point(i).period,per.point(i).profile(2,:),'.-');
end
per_extremum_t=dde_coll_roots(per.point(end),[0,1],'diff',1);
tmin=[0,1]*dde_coll_eva(per.point(end),per_extremum_t(2));
plot(ax4,per_extremum_t(2)*per.point(end).period,tmin,'ko','MarkerFaceColor','k','MarkerSize',8);
yline(ax4,0,'LineWidth',1);
ylim(ax4,[-0.1,1.7]);
axis(ax4,'equal');grid(ax4,'on');
ind_slope1=find(arrayfun(@(p)~isempty(dde_coll_roots(p,[0,1],'diff',1,'res',p.period)),per.point(2:end)),1,'first')+1;
tsl1=dde_coll_roots(per.point(ind_slope1),[0,1],'diff',2);
plot(ax4,tsl1(end)*per.point(ind_slope1).period,[0,1]*dde_coll_eva(per.point(ind_slope1),tsl1(end)),'bo','MarkerFaceColor','b','MarkerSize',8)
%% Continue locus where periodic orbits have point with max delay slope 1
ie=ip;
npar=length(parnames);
[ie.t0,ie.val]=deal(npar+1,npar+2);
graze0=per;
graze0.point=per.point(ind_slope1);
graze0.point.parameter([ie.t0,ie.val])=[tsl1(end),graze0.point(1).period];
min_cond=@(p,pref)dde_extreme_cond(p,[0,1],ie.val,ie.t0,'diff',1,'res',p.period);
[gfuncs,grazing,suc]=ChangeBranchParameters(fcn,graze0,1,...
    'contpar',[ie.tau0,ie.k,ie.t0],...
    'usercond',{min_cond},'outputfuncs',true,...
    'plot_measure',[],'extra_condition',true);
figure(2);
grazing=br_contn(gfuncs,grazing,80,'ax',ax2);
grazing=br_rvers(grazing);
grazing=br_contn(gfuncs,grazing,40,'ax',ax2);
%% plot delays along sonic barrier curve, indicate sonic barrier
lw={'linewidth',2};
figure(5);clf;ax5=gca;xlabel(ax5,'t (period)');ylabel(ax5,'delay s(t))');set(ax5,'FontSize',16);hold(ax5,'on');
for i=1:length(grazing.point)
    plot(ax5,grazing.point(i).mesh*grazing.point(i).period,grazing.point(i).profile(2,:),'-',lw{:});
end
for i=1:length(grazing.point)
    plot(ax5,mod(grazing.point(i).parameter(ie.t0),1)*grazing.point(i).period,...
        [0,1]*dde_coll_eva(grazing.point(i),mod(grazing.point(i).parameter(ie.t0),1)),'ko',lw{:});    
end
axis(ax5,'equal');
axis(ax5,'tight');
ylim(ax5,[0,2])
%%
per_extremum_t=dde_coll_roots(per.point(end),[0,1],'diff',1);
tmin=[0,1]*dde_coll_eva(per.point(end),per_extremum_t(2));
plot(ax4,per_extremum_t(2)*per.point(end).period,tmin,'ko','MarkerFaceColor','k','MarkerSize',8);
yline(ax4,0,'LineWidth',1);
ylim(ax4,[-0.1,1.7]);
axis(ax4,'equal');grid(ax4,'on');
ind_slope1=find(arrayfun(@(p)~isempty(dde_coll_roots(p,[0,1],'diff',1,'res',p.period)),per.point(2:end)),1,'first')+1;
tsl1=dde_coll_roots(per.point(ind_slope1),[0,1],'diff',2);
plot(ax4,tsl1(end)*per.point(ind_slope1).period,[0,1]*dde_coll_eva(per.point(ind_slope1),tsl1(end)),'bo','MarkerFaceColor','b','MarkerSize',8)
%%
figure(3);
Plot2dBranch(grazing,'funcs',gfuncs);
%%
save('poscontrol_results.mat');
%%