%% Four fully coupled delay oscillators - Hopf and periodic orbit bifurcations
% We set the coupling parameter $c_\mathrm{ext}>0$ here. 
%
%%
clear;
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
% For periodic orbits we look at the case of diffusive coupling with
% $c_\mathrm{ext}>0$, where the near-origin equilibrium is stable for small delays.
x0=0.2*ones(n_osc,1);
par([ip.a, ip.c_ext, ip.delta, ip.tau_s,  ip.tau_c])=...
    [1,        1,         0,      1.4,     0.02];
res=funcs.wrap_rhs(x0(:,[1,1,1]),par);
par(ip.delta)=-res(1);
res0=norm(funcs.wrap_rhs(x0(:,[1,1,1]),par))
%% parameter ranges
parbd={'min_bound',[ip.a,0;ip.c_ext,0;   ip.tau_s,0;  ip.tau_c,0],...
       'max_bound',[ip.a,2;ip.c_ext,3;   ip.tau_s,3;  ip.tau_c,4.5],...
       'max_step', [       ip.c_ext,0.05;ip.tau_s,0.1;ip.tau_c,0.1; 0,0.1]};
%% Test detection and branching off from Hopf bifurcation
% Track equilibria for increasing $\tau_c$, which should encounter an
% equivariant Hopf bifurcation. We impose full symchronization to avoid
% singularities along the branch.
ex1234=permfix('stst','eq1234_0',n_osc,'x','perm',[1,2,3,4]);
[fun_1234,eq_1234,suc]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',ip.tau_c,'step',0.1,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',ex1234,'outputfuncs',true);
figure(1);clf;ax1=gca;hold(ax1,'on');set(ax1,'fontsize',16);
xlabel(ax1,'\tau_c');ylabel(ax1,'max(mean(x1))');
eq_1234=br_contn(fun_1234,eq_1234,10,'ax',ax1);
[eq_1234ref,eq_1234tests,eq_1234bifs,eq_1234bifind]=...
    MonitorChange(fun_1234,eq_1234,'min_iterations',5);
fprintf('1st Hopf bifurcation at tau_c=%g\neigenvalues, eigenvectors:',eq_1234bifs(1).parameter(ip.tau_c));
disp(round(eq_1234bifs(1).stability.l0,5))
Plot2dBranch(eq_1234ref,'funcs',fun_1234,'max_nunst',2,'ax',ax1);
legend(ax1,'Fontsize',16);
ylim(ax1,[0,0.5]);
%% Find a few symmetries that one can branch off into
% The function check_hopf below should return 2 zero singular value to
% satisfy the equivariant Hopf bifurcation theorem
check_hopf=@(cnd)SetupHopf(funcs,eq_1234ref,eq_1234bifind,...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
%% Permutation (1234)_0 (time shift 0)
ev1234_0=permfix('stst','p1234_1_4',n_osc,'v',...
    'perm',[1,2,3,4],'rotation',[0,1]);
sv_1234_0=check_hopf(ev1234_0)'
%% Permutation (12)_0 only is not enough: we have four zero singular values for time shifts
%0 or 1/2. Time shift 1/4 would be ok
ev12=permfix('stst','p12',n_osc,'v',...
    'perm',[1,2],'rotation',[0,1]);
sv_12=check_hopf(ev12)'
%% Permutation (1234)_{1/4} (time shift 1/4)
ev1234_1_4=permfix('stst','p1234_1_4',n_osc,...
    'v','perm',[1,2,3,4],'rotation',[1,4]);
sv_1234_1_4=check_hopf(ev1234_1_4)'
%% Permutation  (123)_0 (time shift 0)
ev123=permfix('stst','pv123',n_osc,'v',...
    'perm',[1,2,3],'rotation',[0,1]);
sv_123=check_hopf(ev123)'
%% Permutation (12)_0(34)_0:  1=2 (time shift 0),3=4 (time shift 0), 1=3 (time shift 1/2)
ev12_34=[permfix('stst','pv12_0_1',n_osc,'v',...
    'perm',[1,2],'rotation',[0,1]),...
         permfix('stst','pv34_0_1',n_osc,'v',...
         'perm',[3,4],'rotation',[0,1])];
sv_12_34=check_hopf(ev12_34)'
%% Permutation (12)_{1/2}(34)_0: 1=2 (time shift 1/2),3=4 (time shift 0)
ev12s_34=[permfix('stst','p12_1_2',n_osc,'v',...
    'perm',[1,2],'rotation',[1,2]),...
          permfix('stst','p34_0_1',n_osc,'v',...
          'perm',[3,4],'rotation',[0,1])];
sv_12s_34=check_hopf(ev12s_34)';
%% Permutation (12)_{1/4}: 1=2 (time shift 1/4) is sufficient
ev12_1_4=permfix('stst','p12',n_osc,'v',...
    'perm',[1,2],'rotation',[1,4]);
sv_12_1_4=check_hopf(ev12_1_4)'
%% Branch off toward different periodic solutions
hopfsyms={ev1234_1_4,ev123,ev12s_34,ev12_34};
psol_symmetries.p1234_1_4      =permfix('psol','p1234_1_4',n_osc,'x','perm',[1,2,3,4],'shift',[1,4]);
%
psol_symmetries.p123_0_1       =permfix('psol','p123_0_1', n_osc,'x','perm',[1,2,3],  'shift',[0,1]);
%
psol_symmetries.p12_0_1_34_1_2=[permfix('psol','p12_1_2',  n_osc,'x','perm',[1,2],    'shift',[1,2]),...
                                permfix('psol','p34_0_1',  n_osc,'x','perm',[3,4],    'shift',[0,1])];
%
psol_symmetries.p12_0_1_34_0_1=[permfix('psol','p12_0_1',  n_osc,'x','perm',[1,2],    'shift',[0,1]),...
                                permfix('psol','p34_0_1',  n_osc,'x','perm',[3,4],    'shift',[0,1]),...
                                permfix('psol','p13_1_2',  n_osc,'x','perm',[1,3],    'shift',[1,2])];
symnames=fieldnames(psol_symmetries);
psymrealname.p1234_1_4='p$(1234)_{1/4}$';
psymrealname.p123_0_1='p$(123)_0$';
psymrealname.p12_0_1_34_1_2='p$(12)_0(34)_{1/2}$';
psymrealname.p12_0_1_34_0_1='p$(12)_0(34)_0$';
clear sucpsol
psargs={'extracolumns','auto','extra_condition',true,... % important for imposing symmetry
    'matrix','sparse','eigmatrix','sparse','degree',4,'intervals',60,... % determines accuracy
    'max_step',[0,0.1]}; % stepsize along branch
figure(1);
for i=1:length(symnames)
    sm=symnames{i};
    fprintf('Symmetry: %s\n',sm);
    hopfargs={'SetupHopf.usercond',hopfsyms{i},'SetupHopf.initcond',hopfsyms{i}};
    [funpsol.(sm),psol.(sm),sucpsol.(sm)]=SetupPsol(funcs,eq_1234ref,eq_1234bifind(1),...
        'outputfuncs',true,hopfargs{:},psargs{:},...
        'usercond',psol_symmetries.(sm),'initcond',psol_symmetries.(sm),...
        'plot_measure',{@(p)p.parameter(ip.tau_c),@(p)max(mean(p.profile,1))});
    disp(sucpsol.(sm))
    psol.(sm)=br_contn(funpsol.(sm),psol.(sm),200,'ax',ax1);
    [psolref.(sm),psol_tests.(sm),psolbifs.(sm),psolbifind.(sm)]=...
        MonitorChange(funpsol.(sm),psol.(sm),'min_iterations',5,'range',2:length(psol.(sm).point));
end
%% Continue equivariant Hopf bifurcation in tau_c and tau_s
[ehopffun,ehopf,suc]=SetupHopf(funcs,eq_1234ref,eq_1234bifind,...
    'contpar',[ip.tau_c,ip.tau_s],'dir',ip.tau_c,'step',1e-2,...
    'initcond',ev1234_1_4,'usercond',ev1234_1_4,'outputfuncs',true)
figure(2);clf;ax2=gca;hold(ax2,'on');set(ax2,'fontsize',16);
xlabel(ax2,'\tau_c');ylabel(ax2,'\tau_s');
ehopf=br_contn(ehopffun,ehopf,100,'ax',ax2);
ehopf=br_rvers(ehopf);
ehopf=br_contn(ehopffun,ehopf,100,'ax',ax2);
ehopfstab={'locate_trivial',@(p)1i*p.omega*[1,1,1,-1,-1,-1]};
[ehopfref,ehopf_tests,ehopfbifs,ehopfbifind]=MonitorChange(ehopffun,ehopf,...
    'min_iterations',5,ehopfstab{:});
%% Continue 1st symmetry breaking of (1234)_{1/4} periodic orbit
% First check if not imposing any conditions on the nullvector gives a 1d
% critical subspace
sm=symnames{1};
poev1args={'usercond',psol_symmetries.(sm),'initcond',psol_symmetries.(sm)};
indbif=psolbifind.(sm);
indbif=indbif(end);
svsb=SetupPOEV1(funcs,psolref.(sm),indbif,poev1args{:},'output','sv','nulldim',3)
%% Yes, so try and initialize
poev0=psolref.(sm);
poev0.point=p_remesh(poev0.point(indbif),4,100);
[pffun,pfbr,suc]=SetupPOEV1(funcs,poev0,1,'contpar',[ip.tau_c,ip.tau_s],...
    'dir',ip.tau_c,'step',1e-2,poev1args{:},'plot_measure',[],...
    'max_step',[0,0.5],'use_tangent',true,'min_angle',0.5,...
    'minimal_accuracy',1e-8,'skip_jacobian',1,'newton_max_iterations',10);
%% Continue symmetry breaking
figure(2);
pfbr=br_contn(pffun,pfbr,100,'ax',ax2);
pfbr=br_rvers(pfbr);
pfbr=br_contn(pffun,pfbr,300,'ax',ax2);
[pfbrref,pfbr_tests,pfbrbifs,pfbrbifind]=MonitorChange(pffun,pfbr,'min_iterations',5);
%% continue equilibrium in tau_s, crossing symmetric Hopf bifurcation
ex1234=permfix('stst','eq1234_0',n_osc,'x','perm',[1,2,3,4]);
[fun_1234s,eq_1234s,suc]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',ip.tau_s,'step',0.1,parbd{:},...
    'extracolumns','auto','extra_condition',true,...
    'usercond',ex1234,'outputfuncs',true);
figure(1);clf;ax1=gca;
xlabel(ax1,'\tau_s');ylabel(ax1,'max(mean(x))');
eq_1234s=br_contn(fun_1234s,eq_1234s,10,'ax',ax1);
[eq_1234sref,eq_1234stests,eq_1234sbifs,eq_1234sbifind]=...
    MonitorChange(fun_1234s,eq_1234s,'min_iterations',5);
fprintf('1st Hopf bifurcation at tau_s=%g\neigenvalues, eigenvectors:',eq_1234sbifs(1).parameter(ip.tau_s));
disp(round(eq_1234sbifs(1).stability.l0,5))
%% Continue symmetric Hopf bifurcation in tau_c and tau_s
[shopffun,shopf,suc]=SetupHopf(funcs,eq_1234sref,eq_1234sbifind,...
    'contpar',[ip.tau_c,ip.tau_s],'dir',ip.tau_c,'step',1e-2,...
    'initcond',ex1234,'usercond',ex1234,'outputfuncs',true)
figure(2);
shopf=br_contn(shopffun,shopf,100,'ax',ax2);
shopf=br_rvers(shopf);
shopf=br_contn(shopffun,shopf,100,'ax',ax2);
shopfstab={'locate_trivial',@(p)1i*p.omega*[1,-1]};
[shopfref,shopf_tests,shopfbifs,shopfbifind]=MonitorChange(shopffun,shopf,...
    'min_iterations',5,shopfstab{:});
%% Branch off fully symmetric periodic orbit
fprintf('Full symmetry: (1234)_0\n');
hopfargs={'SetupHopf.usercond',ex1234,'SetupHopf.initcond',ex1234};
px1234=permfix('psol','pp1234_0_1',n_osc,'x','perm',[1,2,3,4],'rotation',[0,1]);
[funspsol,spsol,sucspsol]=SetupPsol(funcs,eq_1234sref,eq_1234sbifind(1),...
        'outputfuncs',true,hopfargs{:},psargs{:},...
        'usercond',px1234,'initcond',px1234);
disp(sucspsol)
figure(1);
spsol=br_contn(funspsol,spsol,200,'ax',ax1);
[spsolref,spsol_tests,spsolbifs,spsolbifind]=...
    MonitorChange(funspsol,spsol,'min_iterations',5,'range',2:length(spsol.point));
%% step near tau_s=2 and then increase tau_c, expecting an equivariant torus bifurcation
it=find(arrayfun(@(p)p.parameter(ip.tau_s),spsol.point)>2,1,'first');
[spsolc,suc]=ChangeBranchParameters(funspsol,spsol,it,'contpar',ip.tau_c,'step',1e-2);
figure(1);clf;ax1=gca;
xlabel(ax1,'\tau_c');ylabel(ax1,'max(mean(x))');
spsolc=br_contn(funspsol,spsolc,200,'ax',ax1);
[spsolcref,spsolc_tests,spsolcbifs,spsolcbifind]=...
    MonitorChange(funspsol,spsolc,'min_iterations',5);
%% Impose additional symmetries on critical Floquet eigenfunctions
% we try a few ones below, looking for symmetries that leave a 2d nullspace
% not imposing any symmetry leaves a 6d nullspace
indbif=spsolcbifind(1);
torargs0={'usercond',px1234,'initcond',px1234};
svtor0=SetupTorusBifurcation(funcs,spsolcref,indbif,torargs0{:},'output','sv','nulldim',7)
% (1234) with rotation by 1/4 circle
pv1234=permfix('psol','pv1234_1_4',n_osc,'v','perm',[1,2,3,4],'rotation',[1,4]);
torargs1={'usercond',[px1234,pv1234],'initcond',[px1234,pv1234]};
svtor1=SetupTorusBifurcation(funcs,spsolcref,indbif,torargs1{:},'output','sv','nulldim',7)
% (12) with no rotation, (34) with 1/2 rotation 
pv12=permfix('psol','p12_0_1',n_osc,'v','perm',[1,2]);
pv34s=permfix('psol','p12_0_1',n_osc,'v','perm',[3,4],'rotation',[1,2]);
torargs2={'usercond',[px1234,pv1234],'initcond',[pv12,pv34s]};
svtor2=SetupTorusBifurcation(funcs,spsolcref,indbif,torargs2{:},'output','sv','nulldim',7)
%% track equivariant torus bifurcation in 2 parameters
[trfuncs,trbr,suc]=SetupTorusBifurcation(funcs,spsolcref,indbif,torargs1{:},...
    'contpar',[ip.tau_c,ip.tau_s],'step',-1e-2,'dir',ip.tau_c,...
    'TorusInit.initmethod','svd');
figure(2);
trbr=br_contn(trfuncs,trbr,70,'ax',ax2);
[trbr,trnunst]=br_stabl(trfuncs,trbr);
%%
pn=fieldnames(psolref);
for i=1:length(pn)
    psol.(pn{i})=br_remove_extracolumns(psol.(pn{i}));
    psolref.(pn{i})=br_remove_extracolumns(psolref.(pn{i}));
end
ehopf=br_remove_extracolumns(ehopf);
pfbr=br_remove_extracolumns(pfbr);
eq_1234s=br_remove_extracolumns(eq_1234s);
shopf=br_remove_extracolumns(shopf);
spsol=br_remove_extracolumns(spsol);
spsolc=br_remove_extracolumns(spsolc);
trbr=br_remove_extracolumns(trbr);
%%
save('S4_demo_psol_results.mat')