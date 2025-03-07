function y=rhs_coupled_oscillators(deg,ip,n_osc,x,p,dx,dp)
Cmat=ones(n_osc,n_osc)-n_osc*eye(n_osc);
on=ones(n_osc,1);
tau=@(x,i)reshape(x(:,i+1,:),n_osc,[]);
[x0,x1,x2]=deal(tau(x,0),tau(x,1),tau(x,2));
p=reshape(p,size(p,2),size(p,3));
[d,a,c]=deal(p(ip.delta,:),p(ip.a,:),p(ip.c_ext,:));
if deg==0
    y=d(on,:)-a(on,:).*x1-x0.^3+c(on,:).*(Cmat*x2);
    return
end
dp=reshape(dp,size(dp,2),size(dp,3));
[dx0,dx1,dx2]=deal(tau(dx,0),tau(dx,1),tau(dx,2));
[dd,da,dc]=deal(dp(ip.delta,:),dp(ip.a,:),dp(ip.c_ext,:));
if deg==1
    y=dd(on,:)-da(on,:).*x1-a(on,:).*dx1-...
    3*x0.^2.*dx0+dc(on,:).*(Cmat*x2)+c(on,:).*(Cmat*dx2);
    return
end
y=-2*da(on,:).*dx1-6*x0.*dx0.^2+2*dc(on,:).*(Cmat*dx2);
end
