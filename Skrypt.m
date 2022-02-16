%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Skrypt z realizacją miniprojektu. Numer zadania: 7.2.3(d).%
%Autorzy: Mikołaj Zapotoczny                %
%Figure1: Charakterystyki dynamiczne za pomocą równań stanu%
%i transmitancji                                           %
%Figure2: Charakterystyki statyczne                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%Wartości nominalne%
TzewN=-20;
TkzN=30;
Tw1N=20;
Tw2N=15;
qkN=20000;
cp=1000;
rop=1.2;

%Objętości%
V1=5*5*3;
V2=5*5*3;

%Pojemności%
Cv1=cp*rop*V1;
Cv2=cp*rop*V2;

%Przepływ%
fpN=qkN/(cp*rop*(TkzN-Tw1N));

%Wartości początkowe%
Tzew0=   TzewN     ;
Tkz0 =   TkzN      ;
fp0  =   fpN       ;

%Skoki%
dfp0   =0*fpN;
dTkz0  =0;
dTzew0 =0;

%Zmienne pomocnicze%
xpN=cp*rop*fpN;
xp0=cp*rop*fp0;

%Wzory na temperatury przejściowe%
Ks1=0.6*qkN/(Tw1N-TzewN);
K0 =(qkN - Ks1*(Tw1N-TzewN))/(Tw1N-Tw2N);
Ks2=(K0*(Tw1N-Tw2N)+xpN*Tw1N-xpN*Tw2N)/(Tw2N-TzewN);

%Wzory na wartości początkowe Tw1 i Tw2%
Tw10=(((xp0*Tkz0 + Ks1*Tzew0)*(K0+xp0+Ks2))+ Ks2*Tzew0*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));
Tw20=(Tw10*(K0+xp0) + Ks2*Tzew0)/(K0 + xp0 + Ks2);


%%%%%%%%%%%%%%%%%%%%%%%Równania stanu%%%%%%%%%%%%%%%%%%%%%%%
A=[((-xp0-Ks1-K0)/Cv1) , (K0/Cv1) ; ((K0+xp0)/Cv2) , ((-K0-xp0-Ks2)/Cv2)];
B=[ (xp0/Cv1) ,    (Ks1/Cv1)   ;  0,        (Ks2/Cv2)     ];
C=[1,0;0,1];
D=zeros(2);

ob2=ss(A,B,C,D);
u0=[Tkz0;Tzew0];
du=[dTkz0,dTzew0];
x0=-(A^(-1))*B*u0;
[y,t]=step(ob2);


%%%%%%%%%%%%%%%%%%%%%%%Transmitancja%%%%%%%%%%%%%%%%%%%%%%%%
s=tf('s');
M=(Cv1*Cv2)*s^2+(Cv1*Ks2+Cv1*K0+Cv1*xp0+K0*Cv2+Ks1*Cv2+xp0*Cv2)*s+(xp0^2+xp0*Ks2+xp0*Ks1+K0*Ks1+Ks1*Ks2+K0*xp0+K0*Ks2);

G11=(Cv2*xp0*s + (xp0+K0+Ks2)*xp0)/M;
G12=(Ks1*Cv1*s + Ks1*(xp0+Ks2+K0)+K0*Ks2)/M;
G21=((K0+xp0)*xp0)/M;
G22=(Cv1*Ks2*s + Ks2*(xp0+Ks1+K0)+K0*Ks1+xp0*Ks1)/M;

G=[G11,G12;G21,G22];
[y,t]=step(G);


%%%%%%%%%%%%%%%%%%%%%%%%%Wykresy%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(221), hold on, grid on, plot(t,Tw10+y(:,1,1)*du(1),'g-'), plot(t,x0(1)+y(:,1,1)*du(1),'k--');
title ('Tw10, dTkz0');
xlabel ('Czas [t]');
ylabel ('Temperatura [^oC]');
legend ('Transmitancja','Równania stanu')

subplot(222), hold on, grid on, plot(t,Tw10+y(:,1,2)*du(2),'g-'), plot(t,x0(1)+y(:,1,2)*du(2),'k--');
title ('Tw10, dTzew0');
xlabel ('Czas [t]');
ylabel ('Temperatura [^oC]');
legend ('Transmitancja','Równania stanu')

subplot(223), hold on, grid on, plot(t,Tw20+y(:,2,1)*du(1),'g-'), plot(t,x0(2)+y(:,2,1)*du(1),'k--');
title ('Tw20, dTkz0');
xlabel ('Czas [t]');
ylabel ('Temperatura [^oC]');
legend ('Transmitancja','Równania stanu')

subplot(224), hold on, grid on, plot(t,Tw20+y(:,2,2)*du(2),'g-'), plot(t,x0(2)+y(:,2,2)*du(2),'k--');
title ('Tw20, dTzew0');
xlabel ('Czas [t]');
ylabel ('Temperatura [^oC]');
legend ('Transmitancja','Równania stanu')


%%%%%%%%%%%%%%Charakterystyki statyczne%%%%%%%%%%%%%%%%%%%%
figure(2);
subplot (321);
Tzew1=[-25:0.1:25];
Tw10=(((xp0*Tkz0 + Ks1*Tzew1)*(K0+xp0+Ks2))+ Ks2*Tzew1*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));

plot(Tzew1,Tw10,'k');
grid on;
hold on;
title ('Tzew0 od Tw1');
xlabel ('Tzew0');
ylabel ('Tw1');
plot (TzewN,Tw1N ,'ro');
legend ('Charakterystyka statyczna','Punkt nominalny','Location','southeast');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot (323);
Tkz1=[-15:0.1:35];
Tw10=(((xp0*Tkz1 + Ks1*Tzew0)*(K0+xp0+Ks2))+ Ks2*Tzew0*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));

plot(Tkz1,Tw10,'k');
grid on;
hold on;
title ('Tkz0 od Tw1');
xlabel ('Tkz0');
ylabel ('Tw1');
plot (TkzN,Tw1N ,'ro');
legend ('Charakterystyka statyczna','Punkt nominalny','Location','southeast');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot (325);
fp1=[0:0.1:50];
xp1=cp*rop*fp1;

k=size(fp1,2);
for i=1:k
Tw10(i)=(((xp1(i)*Tkz0 + Ks1*Tzew0).*(K0+xp1(i)+Ks2))+ Ks2*Tzew0*K0)/((xp1(i)+Ks1+K0).*(K0+xp1(i)+Ks2)-(K0^2+xp1(i)*K0));
end

plot(fp1,Tw10,'k-');
grid on;
hold on;
title ('fp0 od Tw1');
xlabel ('fp0');
ylabel ('Tw1');
plot (fpN,Tw1N,'ro');
legend ('Charakterystyka statyczna','Punkt nominalny','Location','southeast');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot (322);
Tzew1=[-25:0.1:25];
Tw10=(((xp0*Tkz0 + Ks1*Tzew1)*(K0+xp0+Ks2))+ Ks2*Tzew1*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));
Tw20=(Tw10*(K0+xp0) + Ks2*Tzew1)/(K0 + xp0 + Ks2);

plot(Tzew1,Tw20,'k');
grid on;
hold on;
title ('Tzew0 od Tw2');
xlabel ('Tzew0');
ylabel ('Tw2');
plot (TzewN,Tw2N,'ro');
legend ('Charakterystyka statyczna','Punkt nominalny','Location','southeast');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot (324);
Tkz1=[-15:0.1:35];
Tw10=(((xp0*Tkz1 + Ks1*Tzew0)*(K0+xp0+Ks2))+ Ks2*Tzew0*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));
Tw20=(Tw10*(K0+xp0) + Ks2*Tzew0)/(K0 + xp0 + Ks2);

plot(Tkz1,Tw20,'k');
grid on;
hold on;
title ('Tkz0 od Tw2');
xlabel ('Tkz0');
ylabel ('Tw2');
plot (TkzN,Tw2N,'ro');
legend ('Charakterystyka statyczna','Punkt nominalny','Location','southeast');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot (326);
fp1=[0:0.1:50];
xp1=cp*rop*fp1;

k=size(fp1,2);
for i=1:k
Tw10(i)=(((xp1(i)*Tkz0 + Ks1*Tzew0).*(K0+xp1(i)+Ks2))+ Ks2*Tzew0*K0)/((xp1(i)+Ks1+K0).*(K0+xp1(i)+Ks2)-(K0^2+xp1(i)*K0));
Tw20(i)=(Tw10(i)*(K0+xp1(i)) + Ks2*Tzew0)/(K0 + xp1(i) + Ks2);
end

plot(fp1,Tw20,'k-');
grid on;
hold on;
title ('fp0 od Tw2');
xlabel ('fp0');
ylabel ('Tw2');
plot (fpN,Tw2N,'ro');
legend ('Charakterystyka statyczna','Punkt nominalny','Location','southeast');





