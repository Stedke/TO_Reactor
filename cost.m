function Q = cost(u, tau, x0, h0,H,A,zak,k0,E0R) %x0 to stan aktualny
x=x0;
dtau=diff(tau);
n= ceil(dtau/h0);
for j=1:length(dtau)
    h=dtau(j)/n(j);
    h2=h/2;h3=h/3;h6=h/6;
    for i=1:n(j)
        k1=rhs(x,u(:,j),H,A,zak,k0,E0R);
        x1=x+h2*transpose(k1);
        k2=rhs(x1,u(:,j),H,A,zak,k0,E0R); 
        x2=x+h2*transpose(k2);
        k3=rhs(x2,u(:,j),H,A,zak,k0,E0R);
        x3=x+h*transpose(k3);
        k4=rhs(x3,u(:,j),H,A,zak,k0,E0R); 
        x=x+h3*(transpose(k2)+transpose(k3))+h6*(transpose(k1)+transpose(k4));
    end
end

%wartosci z pupy:
T=0.4;
%takie same jak w rob:
ka=10^3;

Q=T+0.5*ka*(x(end,1)-2.13959274764266)^2+0.5*ka*(x(end,2)-1.09030127640364)^2+0.5*ka*(x(end,3)-387.35)^2+0.5*ka*(x(end,4)-386.0655084902178)^2;
end
