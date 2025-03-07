%% An extreme example
% Consider the linear DDE
% 
% $$x'(t)=-x(t)+y(t-\tau)\mbox{,\quad} y'(t)=-y(t)\mbox{.}$$
%
% This DDE has only one finite eigenvalue |-1| of double algebraic
% multiplicity. All other eigenvalues have disappeared at minus infinity,
% indpendent of the delay |tau|. For small delays |tau| all methods
% identify the spectrum correctly, but for larger delays spurious
% eigenvalues occupy part of the complex plane |-1 < Re z < 0| and small
% (in modulus) imaginary parts. This problem can not be remedied by
% refining the discretization.
clear;
ddebiftool_path([pwd(),'/ddebiftool_25_03_07']);
format compact
format short g
evbase=-1;
A=cat(3,eye(2)*evbase,diag(1,1));
tau=20;
lw={'linewidth',2};
mth=getfield(df_mthod('stst','cheb'),'stability');
% set method parameters such that matrix is large and
% all eigenvalues of the discretized system are reported
mth.minimal_real_part=-1.2;
mth.max_number_of_eigenvalues=Inf;
mth.remove_unconverged_roots=0;
mth.maxsize=2050;
mth.discard_accuracy_factor=Inf;
stability=dde_stst_eig_cheb(A,tau,'method',mth);
[~,ix]=sort(abs(stability.l0-evbase),'ascend');
figure(1);clf;hold on;
pc2=plot(real(stability.l0(ix(1:2))),imag(stability.l0(ix(1:2))),...
    'bo',lw{:},'DisplayName','2 ev closest to -1');
po=plot(real(stability.l0(ix(3:end))),imag(stability.l0(ix(3:end))),...
    'rx',lw{:},'DisplayName','spurious');
title(sprintf('pseudospectral discretization, tau=%d',tau));
xlabel('Re');
ylabel('Im');
grid on
legend([pc2,po]);
hold off
set(gca,lw{:},'fontsize',18,'fontname','courier','fontweight','bold');