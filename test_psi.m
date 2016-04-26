function dQ0 = test_psi(u, tau, x0, h0,H,A,zak,k0,E0R)

Q=cost(u,tau,x0,h0,H,A,zak,k0,E0R);
ep=1e-6;
for i=1:length(x0);
    x0_=x0;
    x0_(i)=x0(i)+ep;
    Q_=cost(u,tau,x0_,h0,H,A,zak,k0,E0R);
    dQ0(i,1)=(Q_-Q)/ep;
end
end