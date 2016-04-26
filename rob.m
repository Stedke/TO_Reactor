% Main 
clear all; close all;
clc
format COMPACT;
format LONG e; 

%%parametry
k0=[1.287*10^12 1.287*10^12 9.043*10^6]; %k3 dzielone przez 1000 bo litry
E0R=[9758.3 9758.3 8560];
Hrab=4200;
Hrbc=-11000;
Hrad=-41850;
Ro=0.9342;
Cp=3.01;
kw=4032;
Ar=0.215;
Vr=0.01;
mk=5;
Cpk=2;

zak=[5.1 293.15];
H=[(Hrab/(Ro*Cp))*10^(-3) (Hrbc/(Ro*Cp))*10^(-3) (Hrad/(Ro*Cp))*10^(-3)]; %dzielenie przez 10^3 bo litry
A=[((kw*Ar)/(Ro*Cp*Vr))*10^(-3) (kw*Ar)/(mk*Cpk)];

%% DANE 
x0 = [5.1 0 293.15 293.15 0]; %war poczatkowe - wektor poziomy
T = 0.4;      %czas Symulacji - zmienic tez w cost.m
tau = linspace(0,T,11); % wektor pionowy, zabiera du¿o czasu

param1=3;   %zmiana wartosci pierwszego sterowania (wplywu substratu)
param2=0;  %zmiana wartosci drugiegoo sterowania (odp³ywu ciep³a)

u1=param1*ones(1,10); %wektor kolumnowy wype³niony jedynkami
u2=param2*ones(1,10); %wektor kolumnowy wype³niony jedynkami

u = [u1;u2];  % ci¹g wartoœci sterowania 

%% Dane do algorytmu symulacyjnego
h0=0.0001; % Krok dyskretyzacji RK4 , podstawowy

dtau = diff (tau);
n = ceil (dtau/h0); % ile razy sie zmiesci w przedziale strukturalnym
cn=cumsum( [1 n] ) ;  % a=[1 2 3 4] cumsum(a)->[1 3 6 10]

x = zeros(cn(end), 5);  % 2- tyle ile zmiennych stanu
x(1,1)=x0(1);
x(1,2)=x0(2);
x(1,3)=x0(3);
x(1,4)=x0(4);
t = zeros(cn(end), 1); % czas

%% proces rozwiazywania
% 1 petla - po przedzialach strukturalnych (o stalej wartosci)
for j = 1: length(tau)-1   
    h = dtau(j)/n(j); % h - krok zmienny
    h2 = h/2;
    h3 = h/3;
    h6 = h/6;
    for i = cn(j) : cn(j+1)-1
    % specyficzna czesc rk4
    k1=h*rhs([x(i,1);x(i,2);x(i,3);x(i,4);x(i,5)],u(:,j),H,A,zak,k0,E0R);
    xtemp=[x(i,1);x(i,2);x(i,3);x(i,4);x(i,5)]+k1./2.0;
    ttemp=t(i)+h2;
    k2=h*rhs(xtemp,u(:,j),H,A,zak,k0,E0R);
    xtemp=[x(i,1);x(i,2);x(i,3);x(i,4);x(i,5)]+k2./2.0;
    ttemp=h2;
    k3=h*rhs(xtemp,u(:,j),H,A,zak,k0,E0R);
    xtemp=[x(i,1);x(i,2);x(i,3);x(i,4);x(i,5)]+k3;
    ttemp=t(i)+h;
    k4=h*rhs(xtemp,u(:,j),H,A,zak,k0,E0R);
    
    x(i+1,1)=x(i,1)+(1.0 / 6.0)*(k1(1) + 2.0*k2(1) + 2.0*k3(1) + k4(1));
    x(i+1,2)=x(i,2)+(1.0 / 6.0)*(k1(2) + 2.0*k2(2) + 2.0*k3(2) + k4(2));
    x(i+1,3)=x(i,3)+(1.0 / 6.0)*(k1(3) + 2.0*k2(3) + 2.0*k3(3) + k4(3));
    x(i+1,4)=x(i,4)+(1.0 / 6.0)*(k1(4) + 2.0*k2(4) + 2.0*k3(4) + k4(4));
    x(i+1,5)=x(i,5)+(1.0 / 6.0)*(k1(5) + 2.0*k2(5) + 2.0*k3(5) + k4(5));
    t(i+1)=t(i)+h;
    end
end

%%Petla w tyl
ka=10^3; %wspolczynnik dla kary - zmienic tez w cost.m

psi=zeros(cn(end),5);
psi(end,:)= [-ka*(x(end,1)-2.13959274764266) -ka*(x(end,2)-1.09030127640364) -ka*(x(end,3)-387.35) -ka*(x(end,4)-386.0655084902178) -1];


for j=length(dtau):-1:1
    h=dtau(j)/n(j);
    h2=h/2; h3=h/3; h6=h/6;
    
    for i=cn(j+1):-1:cn(j)+1
        z=[x(i,:) psi(i,:)];
        k1=rhsa(z,u(:,j),H,A,zak,k0,E0R); 
        z1=z-h2*transpose(k1);
        k2=rhsa(z1,u(:,j),H,A,zak,k0,E0R); 
        z2=z-h2*transpose(k2);
        k3=rhsa(z2,u(:,j),H,A,zak,k0,E0R);
        z3=z-h*transpose(k3);
        k4=rhsa(z3,u(:,j),H,A,zak,k0,E0R); 
        z=z-h3*(transpose(k2)+transpose(k3))-h6*(transpose(k1)+transpose(k4));
        psi(i-1,:)=z(6:10);
    end
end

dQ0=test_psi(u, tau, x0, h0,H,A,zak,k0,E0R);
disp([dQ0 -psi(1,:)']);

figure('units','normalized','outerposition',[0 0 1 1])
figure(1)

subplot(2,2,1);
plot(t,x(:,1));
title('x1 Ca');
xlabel('t[h]');
ylabel('Ca[mol]');
xlim([0 T]);
grid on;

subplot(2,2,2);
plot(t,x(:,2));
title('x2 Cb');
xlabel('t[h]');
ylabel('Cb[mol]');
xlim([0 T]);
grid on;

subplot(2,2,3);
plot(t,x(:,3));
title('x3 T');
xlabel('t[h]');
ylabel('T[K]');
xlim([0 T]);
grid on;

subplot(2,2,4);
plot(t,x(:,4));
title('x4 Tk');
xlabel('t[h]');
ylabel('T[K]');
xlim([0 T]);
grid on;

figure('units','normalized','outerposition',[0 0 1 1])
figure(2);

subplot(2,2,1);
plot(t,psi(:,1));
title('psi1');
xlabel('');
ylabel('');
xlim([0 T]);
grid on;

subplot(2,2,2);
plot(t,psi(:,2));
title('psi2');
xlabel('');
ylabel('');
xlim([0 T]);
grid on;

subplot(2,2,3);
plot(t,psi(:,3));
title('psi3');
xlabel('');
ylabel('');
xlim([0 T]);
grid on;

subplot(2,2,4);
plot(t,psi(:,4));
title('psi4');
xlabel('');
ylabel('');
xlim([0 T]);
grid on;

