%% Four fully coupled delay oscillators - bifurcations of steady states
% We set the coupling parameter $c_\mathrm{ext}<0$ here. 
%
%%
clear
ddebiftool_path([pwd(),'/ddebiftool_25_03_07']);
format compact
format short g
%% Define system of equations of 4 coupled scalar oscillators
% The parameter $\delta$ breaks the symmetry $x\mapsto -x$.
%
% $$ x_i(t)=\delta-a x_i(t-\tau_s)-x_i^3-c_\mathrm{ext}\sum_{j\neq i}x_j(t-\tau_c)-x_i(t-tau_c)$$
%
n_osc=4;
parnames={'a','c_ext','delta','tau_s','tau_c'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
Cmat=ones(n_osc,n_osc)-n_osc*eye(n_osc);
on=ones(n_osc,1);
f=@(x0,x1,x2,a,c,d)d(on,:)-a(on,:).*x1-x0.^3+c(on,:).*(Cmat*x2);
df={@(x0,x1,x2,a,c,d,dx0,dx1,dx2,da,dc,dd)...
    dd(on,:)-da(on,:).*x1-a(on,:).*dx1-3*x0.^2.*dx0+dc(on,:).*(Cmat*x2)+c(on,:).*(Cmat*dx2),...
    @(x0,x1,x2,d,a,c,dx0,dx1,dx2,da,dc,dd)...
    -2*da(on,:).*dx1-6*x0.*dx0.^2+2*dc(on,:).*(Cmat*dx2)};
funcs=set_funcs('sys_rhs',{'delay',[ip.a,ip.c_ext,ip.delta],f},'sys_dirderi',df,...)
    'sys_tau',@()[ip.tau_s,ip.tau_c],'x_vectorized',true,'p_vectorized',true);
%% Initial guess and parameters
% For symmetry breaking of equilibria we look at the case
% $c_\mathrm{ext}<0$, where most emerging equilibria are stable for small delay.
x0=0.05*ones(n_osc,1);
par([ip.a,    ip.c_ext, ip.delta, ip.tau_s,  ip.tau_c])=...
    [n_osc+1,   -1,         0,      0.025,     0.02];
res=funcs.wrap_rhs(x0(:,[1,1,1]),par);
par(ip.delta)=-res(1);
res0=norm(funcs.wrap_rhs(x0(:,[1,1,1]),par))
%% Parameter ranges
parbd={'min_bound',[ip.a,-n_osc;ip.c_ext,-2;ip.tau_s,0; ip.tau_c,0],...
    'max_bound',[ip.a,n_osc+6;ip.c_ext,3;ip.tau_s,3; ip.tau_c,5],...
    'max_step',[ip.c_ext,0.05; ip.tau_s,0.1;ip.tau_c,0.1; 0,0.1]};
%% Test detection and branching off for equilibrium symmetry breaking
% Track equilibria for decreasing a, which should encounter a symmetry
% breaking. We impose full symchronization to avoid singularities along the
% branch. Delays do not play a role here.
ex1234=permfix('stst','eq1234_0',n_osc,'x','perm',[1,2,3,4]);
[fun_1234,eq_1234,suc]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',ip.a,'step',-0.1,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',ex1234,'outputfuncs',true);
figure(1);clf;ax1=gca;set(ax1,'fontsize',16);
xlabel(ax1,'a');ylabel(ax1,'x');
eq_1234=br_contn(fun_1234,eq_1234,100,'ax',ax1);
[eq_1234ref,eq_1234tests,eq_1234bifs,eq_1234bifind]=...
    MonitorChange(fun_1234,eq_1234,'min_iterations',5);
fprintf('1st Symmetry breaking at a=%g\neigenvalues, eigenvectors:',eq_1234bifs(1).parameter(ip.a));
disp(eq_1234bifs(1).stability.l0)
disp(eq_1234bifs(1).stability.v)
%% Find a few symmetries that one can branch off into
% check_sb below should return 1 zero singular value to satisfy algebraic
% branching lemma
check_sb=@(cnd)SetupFold(funcs,eq_1234ref,eq_1234bifind(1),...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
% 1=2=3
ev123=permfix('stst','p123',n_osc,'v','perm',[1,2,3]);
sv_123=check_sb(ev123)'
%1=2,3=4
ev12_34=[permfix('stst','p12',n_osc,'v','perm',[1,2]),...
         permfix('stst','p34',n_osc,'v','perm',[3,4])];
sv_12_34=check_sb(ev12_34)'
%1=2 only is not enough: we have two zero singular values
ev12=permfix('stst','p12',n_osc,'v','perm',[1,2]);
sv_12=check_sb(ev12)'
%1=2=3=4 is too much: we have no zero singular values
ev1234=permfix('stst','p1234',n_osc,'v','perm',[1,2,3,4]);
sv_1234=check_sb(ev1234)'
%% Continue symmetry breaking in 2 parameters
cnd_sb=[permfix('stst','px1234',n_osc,'x','perm',[1,2,3,4]);...
    permfix('stst','pv123',n_osc,'v','perm',[1,2,3])];
[fun_sb,eq_sb,suc]=SetupFold(funcs,eq_1234ref,eq_1234bifind(1),...
    'initcond',cnd_sb,'usercond',cnd_sb,'outputfuncs',true,...
    'contpar',[ip.a,ip.c_ext],'dir',ip.a,'step',0.1);
fprintf('Initialize two-parameter continuation of  symmetry breaking\n');
disp(eq_sb)
figure(2);clf;ax2=gca;set(ax2,'fontsize',16);
xlabel(ax2,'a');ylabel(ax2,'c_{ext}');
eq_sb=br_contn(fun_sb,eq_sb,100,'ax',ax2);
eq_sb=br_rvers(eq_sb);
eq_sb=br_contn(fun_sb,eq_sb,100,'ax',ax2);
[eq_sbref,eq_sbtests,eq_sbifs,eq_sbbifind]=MonitorChange(fun_sb,eq_sb,...
    'printlevel',1,'min_iterations',5,'output','branch',...
    'locate_trivial',@(p)[0,0,0]);
%% Branch off in two different directions: 1st 1=2=3
cnd_sb=[permfix('stst','px1234',n_osc,'x','perm',[1,2,3,4]);...
    permfix('stst','pv123',n_osc,'v','perm',[1,2,3])];
cnd_123=permfix('stst','px123',n_osc,'x','perm',[1,2,3]);
sbargs={'SetupFold.usercond',cnd_sb,'SetupFold.initcond',cnd_sb};
[fun_eq123,eq_123,suc]=SetupStstFrom_fold(funcs,eq_1234ref,eq_1234bifind(1),...
    'usercond',cnd_123,'outputfuncs',true,...
    'step',0.1,sbargs{:});
fprintf('Branch off from symmetry breaking, keeping (123) symmetry\n');
disp(eq_123)
figure(1);
eq_123=br_contn(fun_eq123,eq_123,100,'ax',ax1);
eq_123=br_rvers(eq_123);
eq_123=br_contn(fun_eq123,eq_123,100,'ax',ax1);
%% The symmetry gaining is detected twice below. The other bifurcations are genuine
[eq_123ref,eq_123tests,eq_123bifs,eq_123bifind]=...
    MonitorChange(fun_eq123,eq_123,'min_iterations',5);
%% Branch off in two different directions: 2nd 1=2, 3=4
cnd_sb=[permfix('stst','px1234',n_osc,'x','perm',[1,2,3,4]),...
        permfix('stst','p12',n_osc,'v','perm',[1,2]),...
        permfix('stst','p34',n_osc,'v','perm',[3,4])];
cnd_12_34=[permfix('stst','p12',n_osc,'x','perm',[1,2]),...
           permfix('stst','p34',n_osc,'x','perm',[3,4])];
sbargs={'SetupFold.usercond',cnd_sb,'SetupFold.initcond',cnd_sb};
[fun_eq12_34,eq_12_34,suc]=SetupStstFrom_fold(funcs,eq_1234ref,eq_1234bifind(1),...
    'usercond',cnd_12_34,'outputfuncs',true,...
    'step',0.1,sbargs{:});
fprintf('Branch off from symmetry breaking, keeping (12)(34) symmetry\n');
disp(eq_12_34);
figure(1);
eq_12_34=br_contn(fun_eq12_34,eq_12_34,100,'ax',ax1);
eq_12_34=br_rvers(eq_12_34);
eq_12_34=br_contn(fun_eq12_34,eq_12_34,100,'ax',ax1);
%% The symmetry gaining is detected twice below. The other bifurcations are genuine
[eq_12_34ref,eq_12_34tests,eq_12_34bifs,eq_12_34bifind]=...
    MonitorChange(fun_eq12_34,eq_12_34,'min_iterations',5);
%% Continue secondary symmetry breaking in 2 parameters
cnd_sb2=[permfix('stst','px123',n_osc,'x','perm',[1,2,3]);...
    permfix('stst','pv12',n_osc,'v','perm',[1,2])];
[fun_sb2,eq_sb2,suc]=SetupFold(funcs,eq_123ref,eq_123bifind(1),...
    'initcond',cnd_sb2,'usercond',cnd_sb2,'outputfuncs',true,...
    'contpar',[ip.a,ip.c_ext],'dir',ip.a,'step',0.1);
fprintf('Initialize two-parameter continuation of secondary symmetry breaking\n');
disp(eq_sb2)
figure(2);
eq_sb2=br_contn(fun_sb2,eq_sb2,100,'ax',ax2);
eq_sb2=br_rvers(eq_sb2);
eq_sb2=br_contn(fun_sb2,eq_sb2,100,'ax',ax2);
[eq_sb2ref,eq_sb2tests,eq_sb2bifs,eq_sb2bifind]=MonitorChange(fun_sb2,eq_sb2,...
    'printlevel',1,'min_iterations',5,'output','branch',...
    'locate_trivial',@(p)[0,0]);
%% Continue 2nd secondary symmetry breaking in 2 parameters
[fun_sb3,eq_sb3,suc]=SetupFold(funcs,eq_123ref,eq_123bifind(end),...
    'initcond',cnd_sb2,'usercond',cnd_sb2,'outputfuncs',true,...
    'contpar',[ip.a,ip.c_ext],'dir',ip.a,'step',0.1);
fprintf('Initialize two-parameter continuation of secondary symmetry breaking\n');
disp(eq_sb3)
figure(2);
eq_sb3=br_contn(fun_sb3,eq_sb3,100,'ax',ax2);
eq_sb3=br_rvers(eq_sb3);
eq_sb3=br_contn(fun_sb2,eq_sb3,100,'ax',ax2);
[eq_sb3ref,eq_sb3tests,eq_sb3bifs,eq_sb3bifind]=MonitorChange(fun_sb3,eq_sb3,...
    'printlevel',1,'min_iterations',5,'output','branch',...
    'locate_trivial',@(p)[0,0]);
%% Continue third symmetry breaking in 2 parameters
cnd_sb3=[permfix('stst','px12',n_osc,'x','perm',[1,2]),...
         permfix('stst','px34',n_osc,'x','perm',[3,4])];
[fun_sb4,eq_sb4,suc]=SetupFold(funcs,eq_12_34ref,eq_12_34bifind(1),...
    'initcond',cnd_sb3,'usercond',cnd_sb3,'outputfuncs',true,...
    'contpar',[ip.a,ip.c_ext],'dir',ip.a,'step',0.1);
fprintf('Initialize two-parameter continuation of secondary symmetry breaking\n');
disp(eq_sb4)
figure(2);
eq_sb4=br_contn(fun_sb4,eq_sb4,100,'ax',ax2);
eq_sb4=br_rvers(eq_sb4);
eq_sb4=br_contn(fun_sb4,eq_sb4,100,'ax',ax2);
[eq_sb4ref,eq_sb4tests,eq_sb4bifs,eq_sb4bifind]=MonitorChange(fun_sb4,eq_sb4,...
    'printlevel',1,'min_iterations',5,'output','branch',...
    'locate_trivial',@(p)0);
%% 1d bifurcation diagram
figure(3);clf;ax3=gca;hold(ax3,'on');grid(ax3,'on');set(ax3,'fontsize',16);
args={'x',@(p)p.parameter(ip.a),'y',@(p)p.x(1),'ax',ax3};
Plot2dBranch(eq_1234ref,'DisplayName','(1234)','color',[0,0,0.5],args{:});
Plot2dBranch(eq_123ref,'DisplayName','(123)','color',[0,0.5,0],args{:});
Plot2dBranch(eq_12_34ref,'DisplayName','(12)(34)','color',[0.5,0,0],args{:});
xlabel(ax3,'a');ylabel(ax3,'x');
legend(ax3,'Location','EastOutside');
%% 2d bifurcation diagram
figure(4);clf;ax4=gca;hold(ax4,'on');grid(ax4,'on');set(ax4,'fontsize',16);
args={'x',@(p)p.parameter(ip.a),'y',@(p)p.parameter(ip.c_ext),'ax',ax4,'exclude_trivial',true};
Plot2dBranch(eq_sbref,'DisplayName','SB(1234)->3','color',[0,0,0.5],args{:},'locate_trivial',@(p)[0,0,0]);
Plot2dBranch(eq_sb2ref,'DisplayName','SB(123)->2','color',[0,0.5,0],args{:},'locate_trivial',@(p)[0,0]);
Plot2dBranch(eq_sb3ref,'DisplayName','SB(123)->2','color',[0,0.5,0],args{:},'locate_trivial',@(p)[0,0]);
Plot2dBranch(eq_sb4ref,'DisplayName','SB(12)(34)->1','color',[0.5,0,0],args{:},'locate_trivial',@(p)0);
xlabel(ax4,'a');ylabel(ax4,'c_{ext}')
%%
save('S4_demo_stst_results.mat')