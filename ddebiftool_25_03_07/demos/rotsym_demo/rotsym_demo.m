%% Demo of extension for rotational symmetry using Lang-Kobayashi equations
%
% The equation (modelling a laser subject to delayed coherent optical
% feedback) is given as
% 
% $$E'(t)=[1+i\alpha]n(t)E(t)+\eta\exp(i\phi)E(t-\tau), \qquad
%   n'(t)=\epsilon[p-n(t)-(2n(t)+1)\bar E(t)E(t)]$$
% 
% The main bifurcation parameters will be $\eta$ and $\phi$.
%
% *Warning* The functions for the extended systems determining relative
% equilibria (rotating waves) and relative periodic orbits (modulated
% waves), have only support for directional user-provided system derivatves
% up to order 1, and have no support for state-dependent delays!
%
% <html>
% $Id: rotsym_demo.m 369 2019-08-27 00:07:02Z jansieber $
% </html>
%
%% Load DDE-Biftool and extension into Path
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_extra_rotsym/'],...
    [base,'ddebiftool_utilities']);
%% Problem definition using |set_rotfuncs|
% In addition to the user-defined functions |set_rotfuncs| needs the
% matrix |A| generating the rotation and (optional) the rotation as a
% function $\phi\mapsto \exp(A\phi)$. Then the system is assumed to have
% rotational symmetry $exp(A\phi)$ where |A| is anti-symmetric.
A=[0,-1,0; 1,0,0; 0,0,0];
expA=@(phi)reshape([...
     cos(phi);sin(phi);0*phi;...
    -sin(phi);cos(phi);0*phi;...
    0*phi;0*phi;1+0*phi],3,3,[]);
%% Initial values of parameters and parameter indices
parnames={...
    'pump',...    % injection current
    'eta',...     % feedback strength
    'phi',...     % feedback phase
    'tau',...     %  delay
    'alpha',...   % alpha factor
    'epsilon',... % carrier relaxation time
    'omega'};     % rotation velocity
cs=[parnames;num2cell(1:length(parnames))];
ip=struct(cs{:});
par([ip.alpha,ip.pump,ip.epsilon,ip.eta,ip.phi,ip.tau])=...
    [    4,      0.1,    5e-3,     5e-3,   0,     100];
%% Right-hand side and call to |set_rotfuncs|
f=@(x,p)LangKobayashi(x(1,1,:)+1i*x(2,1,:),x(1,2,:)+1i*x(2,2,:),x(3,1,:),p,ip);
%rfuncs=set_rotfuncs('sys_rhs',f,'rotation',A,'exp_rotation',expA,...
%    'sys_tau',@()ip.tau,'x_vectorized',true);
sfuncs=set_symfuncs(@sym_LangKobayashi,'sys_tau',@()ip.tau);
rfuncs=set_rotfuncs('wrap_rhs',dde_sym2fun(@sym_LangKobayashi,'rhs'),...
    'wrap_dirderi',dde_sym2fun(@sym_LangKobayashi,'rhs_dirderi'),...
    'rotation',A,'exp_rotation',expA,...
    'sys_tau',@()ip.tau,'x_vectorized',true,'expA_vectorized',true);
%% Extended functions in rotating coordinates for rotating waves
% A rotating wave (relative equilibrium, RE) is a solution of the form
%
% $$ x(t)=\exp(A\omega t)x_0. $$
%
% The extended functions in |rfuncs| will always treat the user-provided
% system in rotating coordinates.
%
% Original system (by user):
%
% $$x'(t)=f(x(t),x(t-\tau_1),...,x(t-\tau_n))$$
%
% Rotating coordinates: $x(t)=\exp(A\omega t)y(t)$:
%
% $$y'(t)=-A \omega y(t)+f(y(t),\exp(-A\omega\tau_1)y(t-\tau_1),\ldots,
%    \exp(-A\omega\tau_n)y(t-\tau_n))$$
%
% The rotation speed is chosen such that the rotating wave
% $x(t)=\exp(A\omega t)x_0$ is turned into an equilibrium $y(t)=y_0$. This
% is achieved by solving for equilibria of the $y$ equation and adding a
% |sys_cond| (file |rot_cond.m|) to povide an equation determining
% $\omega$. For rotating waves this is
%
% $$ y_\mathrm{ref}^TAy=0 $$
%
% As DDE-Biftool does not give the user's |sys_cond| access to reference
% points, |rot_cond| returns residual 0 and Jacobian $y_0^TA$.
%% Initial guess
% The extension for rotating and modulated waves is *not* able to cope with
% invariant equilibria. Thus, we generate a non-trivial rotating wave as
% our initial guess. For the laser the rotating waves correspond to
% stationary lasing (on state) $E(t)=E_0\exp(i\omega t)$, $n=n_0$.
[E0,n0,par]=LK_init(par,ip);
bd={'max_bound',[ip.eta,1e-2;ip.phi,8],'min_bound',[ip.eta,0; ip.phi,-4]};
opt_inputs=[bd,{'extra_condition',1,'print_residual_info',1}];
stab_inputs={'exclude_trivial',true,'pointtype_list',@rot_pointtype_list};
%% Relative equilibria varying phase |phi| of the delayed feedback
% The standard convenience function |SetupStst| works, but one must add the
% continuation parameter index |length(parameter)| to the index list of
% continuation parameters.
%
% *Warning* DDE-Biftool assumes that the rotation speed is the last
% parameter in the parameter vector!
%
rw_phas=SetupStst(rfuncs,'contpar',[ip.phi,ip.omega],'corpar',ip.omega,...
    'x',[E0;n0],'parameter',par,opt_inputs{:},...
    'max_step',[ip.phi,0.2]);
figure(2);clf;ax=gca;
rw_phas=br_contn(rfuncs,rw_phas,200,'plotaxis',ax);
%% Linear Stability of relative equilibria
% The standard convenience function |GetStability| works. However, all
% relative equilibria have a trivial eigenvalue 0, corresponding to phase
% shift. If $x_0$ is a relative equilibrium (with rotation speed $\omega$)
% then $\exp(A\rho)x_0$ is also a relative equilibrium with the same
% rotation speed $\omega$. One can adapt |GetStability| by setting the
% optional flag |exclude_trivial| to |true| and providing a function for
% locating the trivial eigenvalue: |'locate_trivial',@(p)0|.
[rw_phas,rw_phas_nunst,dom,defect]=br_stabl(rfuncs,rw_phas,stab_inputs{:});
% plot with stability information
p2=arrayfun(@(p)p.parameter(ip.phi),rw_phas.point);
A2=arrayfun(@(p)norm(p.x(1:2)),rw_phas.point);
n2=arrayfun(@(p)norm(p.x(3)),rw_phas.point);
om2=arrayfun(@(p)p.parameter(ip.omega),rw_phas.point);
tdeco={'fontsize',14,'fontweight','bold'};
ldeco={'interpreter','LaTeX','fontsize',18};
figure(1);clf;ax1=gca;hold(ax1,'on');
[eqlg]=Plot2dBranch(rw_phas,stab_inputs{:},'y',@(p)norm(p.x(1:2)));
% detect Hopf & fold bifurcations
ind_hopf=find(abs(diff(rw_phas_nunst))==2);
ind_fold=find(abs(diff(rw_phas_nunst))==1);
plot(ax1,p2(ind_hopf),A2(ind_hopf),'kx','linewidth',2,'DisplayName','Hopf');
plot(ax1,p2(ind_fold),A2(ind_fold),'k*','linewidth',2,'DisplayName','fold');
hold(ax1,'off');
set(ax1,tdeco{:});
xlabel(ax1,'$\phi$',ldeco{:});
ylabel(ax1,'$|E|$',ldeco{:});
%% Modulated waves (Relative periodic orbits, RPOs)
% Modulated waves are solutions of the form
%
% $$ x(t)=\exp(A\omega t)x_0(t)$$
%
% where $x0(t)=x0(t-T)$ for all $t$ and some period $T$. That is, $x(t)$ is
% quasi-periodic, but can be turned into a periodic solution in rotating
% coordinates. The transformation to rotating coordinates is the same as
% for relative equilibria, but the additional condition (|rot_cond|) is
%
% $$\int_0^1 y_\mathrm{ref}^TAy(t) \mathrm{d} t=0$$
%
% where $y_\mathrm{ref}(t)$ is a reference solution. Since DDE-Biftool does
% not give access to reference solutions in user-defined conditions,
% |rot_cond| returns residual 0 and Jacobian $y(t)^TA$.
%% RPOs branching off at 2nd Hopf of REs
% The initialization works with the standard routine. Again, the rotation
% speed needs to be added to the list of continuation parameters. We use
% the optional sparse matrix computation for (minor) spped-up. Since in
% this example RPOs do not get very spiky it is beneficial to use
% high-order interpolation. for high-order interpolation the polynomial
% should not be stored on evenly spaced mesh points. This is indicated with
% the option |'submesh'| |'cheb'|, using Chebyshev nodes in each
% subinterval. For the collocation points the default is |'gauss'|
% (Legendre points). One may also choose |'cheb'| for
% |'collocation_parameters'|. For this problem one may even choose a single
% collocation interval (|'intervals',1|) and a large degree (|'degree',50|)
% with |'submesh','cheb'|. In that case |'matrix', 'sparse'| does not make
% sense, since the matrix will be full.
[rw_phas_per,suc]=SetupPsol(rfuncs,rw_phas,ind_hopf(2),opt_inputs{:},...
    'max_step',[0,0.5;ip.phi,0.1],'print_residual_info',1,'radius',0.02,...
    'intervals',30,'submesh','cheb','degree',3,'matrix','sparse','excludefreqs',0);
if ~suc
    error('Psol initialization failed');
end
figure(2);clf;ax=gca;
rw_phas_per=br_contn(rfuncs,rw_phas_per,120,'plotaxis',ax);
%% Stability of RPOs
% Similar to REs, the RPOs have an additional trivial Floquet multiplier 1.
% That is, overall, RPOs have always a double Floquet multiplier 1. To
% exclude the two Floquet mulitpliers closest to unity from the stability,
% set the optional flag |'exclude_trivial'| to |true| and provide for the
% optional argument |'locate_trivial'| the function |@(p)[1,1]|.
% MonitorChange refines the branch near stability changes ain a range of
% points.
[rw_phas_per_ref,rw_phas_per_nunst,rw_phas_per_bifs,rw_phas_per_bifind]=MonitorChange(rfuncs,rw_phas_per,...
    'range',2:length(rw_phas_per.point),'printlevel',1,stab_inputs{:},'print_residual_info',0);
[rw_phas_per_nunst,dom_per,defect,rw_phas_per_ref.point]=...
    GetStability(rw_phas_per_ref,stab_inputs{:},'funcs',rfuncs);
% plot with stability info
pp=arrayfun(@(p)p.parameter(ip.phi),rw_phas_per_ref.point);
Epow=@(x)sqrt(sum(x(1:2,:).^2,1));
Amax=arrayfun(@(p)max(Epow(p.profile)),rw_phas_per_ref.point);
Amin=arrayfun(@(p)min(Epow(p.profile)),rw_phas_per_ref.point);
omp=arrayfun(@(p)p.parameter(ip.omega),rw_phas_per_ref.point);
% detect period doublin, torus & fold bifurcations
ind_tr=find(abs(diff(rw_phas_per_nunst))==2);
ind_pf=find(diff(rw_phas_per_nunst==1) & real(dom_per(1:end-1))>0);
ind_pd=find(diff(rw_phas_per_nunst==1) & real(dom_per(1:end-1))<0);
figure(1);hold(ax1,'on')
plot(ax1,pp(ind_tr),[Amax(ind_tr);Amin(ind_tr)]','ks','linewidth',2,'DisplayName','torus');
plot(ax1,pp(ind_pf),[Amax(ind_pf);Amin(ind_pf)]','kd','linewidth',2,'DisplayName','POfold');
plot(ax1,pp(ind_pd),[Amax(ind_pd);Amin(ind_pd)]','kp','linewidth',2,'DisplayName','period doubling');
Plot2dBranch(rw_phas_per_ref,'ax',ax1,stab_inputs{:},'y',@(p)max(Epow(p.profile)));
Plot2dBranch(rw_phas_per_ref,'ax',ax1,stab_inputs{:},'y',@(p)min(Epow(p.profile)));
hold(ax1,'off');
axis(ax1,'tight');
%% Plot of "phase portraits" of relative periodic orbits
figure(3);clf;ax3=gca;hold(ax3,'on');
for i=1:length(rw_phas_per.point)
    plot(ax3,Epow(rw_phas_per.point(i).profile),rw_phas_per.point(i).profile(3,:),'.-');
end
hold(ax3,'off');
grid(ax3,'on');
axis(ax3,'tight');
set(ax3,tdeco{:});
xlabel(ax3,'$n$',ldeco{:});
ylabel(ax3,'$|E|$',ldeco{:});
%% Continuation of fold of relative equilibria
% The functions for the extended system are generated by |SetupRWFold|. The
% standard routine |SetupPOfold| has to be modified because the extended
% condition |rot_cond| has to be applied to the derivative, too. Also the
% rotation speed needs to be added to the list of continuation parameters.
% The extended system for fold continuation of REs has one additional
% artificial continuation parameter.
ind_fold=find(abs(diff(rw_phas_nunst))==1);
[foldfuncs,fold1branch,suc]=SetupRWFold(rfuncs,rw_phas,ind_fold(1),...
    'contpar',[ip.phi,ip.eta,ip.omega],opt_inputs{:},...
    'print_residual_info',1,'dir',ip.eta,...
    'max_step',[ip.phi,0.1; ip.eta,0.01]);
%%
figure(2);clf;ax=gca;
fold1branch=br_contn(foldfuncs,fold1branch,200,'plotaxis',ax);
fold1branch=br_rvers(fold1branch);
fold1branch=br_contn(foldfuncs,fold1branch,200,'plotaxis',ax);
[fold1branch,nunst1_fold]=br_stabl(foldfuncs,fold1branch,stab_inputs{:});
%% Continuation of 2nd fold of relative equilibria
[foldfuncs,fold2branch,suc]=SetupRWFold(rfuncs,rw_phas,ind_fold(2),...
    'contpar',[ip.phi,ip.eta,ip.omega],opt_inputs{:},...
    'print_residual_info',1,'dir',ip.eta,...
    'max_step',[ip.phi,0.1; ip.eta,0.01]);
figure(2);
fold2branch=br_contn(foldfuncs,fold2branch,200,'plotaxis',ax);
fold2branch=br_rvers(fold2branch);
fold2branch=br_contn(foldfuncs,fold2branch,200,'plotaxis',ax);
[fold2branch,nunst2_fold]=br_stabl(foldfuncs,fold2branch,stab_inputs{:});
%% Continuation of Hopf bifurcations of relative equilibria
% For Hopf bifurcation continuation the standard routine |SetupHopf| works
% without modification (|SetupRWHopf| is a simple wrapper).
[h1branch,suc]=SetupRWHopf(rfuncs,rw_phas,ind_hopf(1),...
    'contpar',[ip.phi,ip.eta,ip.omega],opt_inputs{:},...
    'print_residual_info',1,'dir',ip.eta);
figure(2);
h1branch=br_contn(rfuncs,h1branch,100,'plotaxis',ax);
h1branch=br_rvers(h1branch);
h1branch=br_contn(rfuncs,h1branch,100,'plotaxis',ax);
[h1branch,nunst1_hopf]=br_stabl(rfuncs,h1branch,stab_inputs{:});
%% Plot all bifurcations of relative equilibria
figure(4);clf;ax4=gca;hold(ax4,'on');
Plot2dBranch({fold1branch,fold2branch},'funcs',foldfuncs,stab_inputs{:});
Plot2dBranch(h1branch,'funcs',rfuncs,stab_inputs{:});
axis(ax4,[-2,8,0,0.009]);
grid(ax4,'on');
set(get(ax4,'legend'),'Location','Best');
set(ax4,tdeco{:});
xlabel(ax4,'$\phi$',ldeco{:});
ylabel(ax4,'$\eta$',ldeco{:});
%% Period doubling of relative POs
% The standard initialization works (wrapped to give it sensible name). The
% rotation speed needs to be added to the list of continuation parameters.
% We set a sharper Newton tolerance (|'minimal_accuracy'|) than default
% since the default (|1e-6|) produced spurious solutions for the fold of
% RPOs.
ind_pd=find(diff(rw_phas_per_nunst==1)&real(dom_per(1:end-1))<0);
[pdfuncs,pdbr,suc]=SetupMWPeriodDoubling(rfuncs,rw_phas_per_ref,ind_pd,...
    'contpar',[ip.phi,ip.eta,ip.omega],opt_inputs{:},...
    'print_residual_info',1,'dir',ip.eta,'minimal_accuracy',1e-8);
%% Continuation and stability along branch
% Here the stability includes the Floquet multiplier -1. Alternatively
% include -1 into the list in |'locate_trivial'|:
% |'locate_trivial',@(p)[1,1,-1]|.
figure(2);
pdbr=br_contn(pdfuncs,pdbr,100,'plotaxis',ax);
pdbr=br_rvers(pdbr);
pdbr=br_contn(pdfuncs,pdbr,120,'plotaxis',ax);
[pdbr,nunst_pd,dom,triv_defect]=br_stabl(pdfuncs,pdbr,stab_inputs{:});
fprintf('max error of Floquet mult close to -1 or 1: %g\n',max(abs(triv_defect)));
%% Fold of relative POs
% Folds of RPOs require a modified initialization routine since the
% condition |rot_cond| needs to be applied to the derivative, too. This
% requires again the introduction of an artificial parameter (automatically
% appended to parameter vector).  The user has to add rotation speed to the
% list of continuation parameters.
% We set a sharper Newton tolerance (|'minimal_accuracy'|) than default
% since the default (|1e-6|) produced spurious solutions. We also trial a
% tangent-based pseudo-arclength continuation with a large maximal stepsize
% but a small permitted deviation from straight angle ('minimal_angle'
% close to 1). This options is only active when 'use_tangent' is true.
[~,ind_mwf]=max(arrayfun(@(p)p.parameter(ip.phi),rw_phas_per_bifs));
ind_mwf=rw_phas_per_bifind(ind_mwf);
[pfoldfuncs,mwfoldbr]=SetupMWFold(rfuncs,rw_phas_per_ref,ind_mwf,...
    'contpar',[ip.phi,ip.eta,ip.omega],opt_inputs{:},'step',-0.01,...
    'print_residual_info',1,'dir',ip.phi,'use_tangent',true,...
    'minimal_accuracy',1e-8,'minimal_angle',0.5,'max_step',[0,0.2]);
%%
figure(2);
mwfoldbr=br_contn(pfoldfuncs,mwfoldbr,200,'plotaxis',ax);
mwfoldbr=br_rvers(mwfoldbr);
mwfoldbr=br_contn(pfoldfuncs,mwfoldbr,200,'plotaxis',ax);
%% Stability of orbits at fold of RPOs
% Here the stability includes the Floquet multiplier -1. Alternatively
% include  another 1 into the list in |'locate_trivial'|:
% |'locate_trivial',@(p)[1,1,1]|.
[mwfoldbr,nunst_mwf,dom,triv_defect]=br_stabl(pfoldfuncs,mwfoldbr,stab_inputs{:});
fprintf('max error of Floquet mult close to 1: %g\n',max(abs(triv_defect)));
%% Re-plot one-parameter and two-parameter bifurcation diagrams
% using |Plot2dBranch|
figure(1);clf;ax1=gca;hold(ax1,'on');
plot_inp={'funcs',rfuncs,'pointtype_list',@rot_pointtype_list,'ax',ax1};
Plot2dBranch(rw_phas,plot_inp{:});
Plot2dBranch(rw_phas_per_ref,plot_inp{:});
plot(ax1,p2(ind_hopf),om2(ind_hopf),'kx','linewidth',2,'DisplayName','Hopf');
plot(ax1,p2(ind_fold),om2(ind_fold),'k*','linewidth',2,'DisplayName','fold');
plot(ax1,pp(ind_tr),omp(ind_tr),'ks','linewidth',2,'DisplayName','torus');
plot(ax1,pp(ind_pf),omp(ind_pf),'kd','linewidth',2,'DisplayName','POfold');
plot(ax1,pp(ind_pd),omp(ind_pd),'kp','linewidth',2,'DisplayName','period doubling');
xlabel(ax1,'$\phi$',ldeco{:});
ylabel(ax1,'$\omega$',ldeco{:});
set(legend(ax1),'location','SouthEast');
title(ax1,sprintf('One-parameter bifurcation diagram, $\\eta$=%g',par(ip.eta)),ldeco{:});
set(ax1,'fontsize',18);
figure(2);clf;ax2=gca;hold(ax2,'on');
plot_inp={'pointtype_list',@rot_pointtype_list,'ax',ax2};
Plot2dBranch(fold1branch,'funcs',foldfuncs,plot_inp{:});
Plot2dBranch(fold2branch,'funcs',foldfuncs,plot_inp{:});
Plot2dBranch(h1branch,'funcs',rfuncs,plot_inp{:});
Plot2dBranch(pdbr,'funcs',pdfuncs,plot_inp{:});
Plot2dBranch(mwfoldbr,'funcs',pfoldfuncs,plot_inp{:});
xlabel(ax2,'$\phi$',ldeco{:});
ylabel(ax2,'$\eta$',ldeco{:});
set(legend(ax2),'location','SouthEast');
title(ax2,'Two-parameter bifurcation diagram',ldeco{:});
set(ax2,'fontsize',18);
%%
save('LKbifs.mat');
