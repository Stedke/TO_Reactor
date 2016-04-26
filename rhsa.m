function [y] = rhsa(z,u,H,A,zak,k0,E0R)
%z(1) do z(5)     stan
%z(6) do z(10)    psi
%z(11) do z(12)   dQ
%H(1) do H(3)     duze H(policzone)
%A(1) do A(2)     duze A
%zak(1) do zak(2) z
%k0(1) do k0(3)   k10 do k30
%E0R(1) do E0R(3) E0/R do E3/R

    k=exp(log(k0)-E0R*(1/z(3)));
    poch=exp(log(k0)-E0R*(1/z(3)));
    
    y = [u(1)*(zak(1)-z(1))-k(1)*z(1)-k(3)*z(1)^2;
        -u(1)*z(2)+k(1)*z(1)-k(2)*z(2);
        u(1)*(zak(2)-z(3))-(k(1)*z(1)*H(1)+k(2)*z(2)*H(2)+k(3)*z(1)^2*H(3))+A(1)*(z(4)-z(3));
        -u(2)+A(2)*(z(3)-z(4));
        u(2)^2;

        -((-u(1)-k(1)-2*k(3)*z(1))*z(6)+k(1)*z(7)+(-k(1)*H(1)-2*z(1)*k(3)*H(3))*z(8));
        -((-u(1)-k(2))*z(7)-k(2)*H(2)*z(8));
        -((-(poch(1)*(z(1)*E0R(1)*(1/z(3)^2)))-(poch(3)*(z(1)^2*E0R(3)*(1/z(3)^2))))*z(6)+((poch(1)*z(1)*(E0R(1)*(1/z(3)^2)))-(poch(2)*z(2)*(E0R(2)*(1/z(3)^2))))*z(7)+(-u(1)-A(1)-(poch(1)*(z(1)*H(1)*E0R(1)*(1/z(3)^2)))-(poch(2)*(z(2)*H(2)*E0R(2)*(1/z(3)^2)))-(poch(3)*(z(1)^2*H(3)*E0R(3)*(1/z(3)^2))))*z(8)+A(2)*z(9));
        -(A(1)*z(8)-A(2)*z(9));
        0;
        
        (zak(1)-z(1))*z(6)-z(2)*z(7)+(zak(2)-z(3))*z(8);
        -z(9);
    ];
%bledy w rozniczkowaniu (rzad,kolumna) z worda:
%(1,3)
%(2,3)
%(3,3)
%(3,1)
end