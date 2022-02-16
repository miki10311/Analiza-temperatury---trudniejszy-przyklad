%7.2.3(d)
clear all;
close all;

steptime=50;

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

xpN=cp*rop*fpN;
xp0=cp*rop*fp0;

%Wzory na Temperatury przejściowe%
Ks1=0.6*qkN/(Tw1N-TzewN);
K0 =(qkN - Ks1*(Tw1N-TzewN))/(Tw1N-Tw2N);
Ks2=(K0*(Tw1N-Tw2N)+xpN*Tw1N-xpN*Tw2N)/(Tw2N-TzewN);

%Wzory na wartości początkowe Tw1 i Tw2%
Tw10=(((xp0*Tkz0 + Ks1*Tzew0)*(K0+xp0+Ks2))+ Ks2*Tzew0*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));
Tw20=(Tw10*(K0+xp0) + Ks2*Tzew0)/(K0 + xp0 + Ks2);

figure(1)
sim("mini_projekt",1000);
plot(ans.tout, ans.Tw1);
grid on;
title('Tw1');
xlabel('czas [s]');
ylabel('temperatura [^oC]');
legend('Tw1');


figure(2)
sim("mini_projekt",1000);
plot(ans.tout, ans.Tw2);
grid on;
title('Tw2');
xlabel('czas [s]');
ylabel('temperatura [^oC]');
legend('Tw2');


%%%%%%%%%%%%%%%%%bloczek state%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[((-xp0-Ks1-K0)/Cv1) , (K0/Cv1) ; ((K0+xp0)/Cv2) , ((-K0-xp0-Ks2)/Cv2)];
B=[ (xp0/Cv1) ,    (Ks1/Cv1)   ;  0,        (Ks2/Cv2)     ];
C=[1,0;0,1];
D=zeros(2);

u=[Tkz0;Tzew0];
Tw=-(A^(-1))*B*u;

figure(3)
sim("state",1000);
plot(ans.tout, ans.Tw1);
grid on;
title('Tw1 - state');
xlabel('czas [s]');
ylabel('temperatura [^oC]');
legend('Tw1');


figure(4)
sim("state",1000);
plot(ans.tout, ans.Tw2);
grid on;
title('Tw2 - state');
xlabel('czas [s]');
ylabel('temperatura [^oC]');
legend('Tw2');




%%%%%%%%%%%%%%%%%transmitancja%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M1=Cv1*Cv2;
M2=(Cv1*Ks2+Cv1*K0+Cv1*xp0+K0*Cv2+Ks1*Cv2+xp0*Cv2);
M3=(xp0^2+xp0*Ks2+xp0*Ks1+K0*Ks1+Ks1*Ks2+K0*xp0+K0*Ks2);

G11a=Cv2*xp0;
G11b=(xp0+K0+Ks2)*xp0;

G12a=Ks1*Cv1;
G12b=Ks1*(xp0+Ks2+K0)+K0*Ks2;

G21a=(K0+xp0)*xp0;

G22a=Cv1*Ks2;
G22b=Ks2*(xp0+Ks1+K0)+K0*Ks1+xp0*Ks1;


figure(5)
sim("transmitancja",1000);
plot(ans.tout, ans.Tw1);
grid on;
title('Tw1 - transmitancja');
xlabel('czas [s]');
ylabel('temperatura [^oC]');
legend('Tw1');


figure(6)
sim("transmitancja",1000);
plot(ans.tout, ans.Tw2);
grid on;
title('Tw2 - transmitancja');
xlabel('czas [s]');
ylabel('temperatura [^oC]');
legend('Tw2');


%%%%%%%%%%%%%%%%razem%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)
hold on;
sim("mini_projekt",1000);
plot(ans.tout, ans.Tw1, 'r-');

sim("state",1000);
plot(ans.tout, ans.Tw1, 'g--');


sim("transmitancja",1000);
plot(ans.tout, ans.Tw1, 'b.');
grid on;
title('Tw1');
xlabel('czas [s]');
ylabel('temperatura [^oC]');
legend('nlin', 'state', 'transmitacja');


figure(8)
hold on;
sim("mini_projekt",1000);
plot(ans.tout, ans.Tw2, 'r-');

sim("state",1000);
plot(ans.tout, ans.Tw2, 'g--');


sim("transmitancja",1000);
plot(ans.tout, ans.Tw2, 'b.');
grid on;
title('Tw2');
xlabel('czas [s]');
ylabel('temperatura [^oC]');
legend('nlin', 'state', 'transmitacja');


%%%%%%%charakterystyki statyczne%%%%%%%%%%%%%%%%%%%%%

subplot (321) ;
Tzew1 =[-25:0.1:25];
Tw10=(((xp0*Tkz0 + Ks1*Tzew1)*(K0+xp0+Ks2))+ Ks2*Tzew1*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));
Tw20=(Tw10*(K0+xp0) + Ks2*Tzew1)/(K0 + xp0 + Ks2);

plot(Tzew1, Tw10,'k');
grid on;
hold on;
title('Tzew0 od Tw1');
xlabel ('Tzew0 ');
ylabel ('Tw1 ');
plot (TzewN , Tw1N , 'ro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot (323) ;
Tkz1 =[-15:0.1:35];
Tw10=(((xp0*Tkz1 + Ks1*Tzew0)*(K0+xp0+Ks2))+ Ks2*Tzew0*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));
Tw20=(Tw10*(K0+xp0) + Ks2*Tzew0)/(K0 + xp0 + Ks2);

plot(Tkz1, Tw10,'k');
grid on;
hold on;
title('Tkz0 od Tw1');
xlabel ('Tkz0 ');
ylabel ('Tw1 ');
plot (TkzN , Tw1N , 'ro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot (325) ;
fp1 =[0:0.1:50];
xp1=cp*rop*fp1;

k=size(fp1,2);

for i=1:k
Tw10(i)=(((xp1(i)*Tkz0 + Ks1*Tzew0).*(K0+xp1(i)+Ks2))+ Ks2*Tzew0*K0)/((xp1(i)+Ks1+K0).*(K0+xp1(i)+Ks2)-(K0^2+xp1(i)*K0));
Tw20(i)=(Tw10(i)*(K0+xp1(i)) + Ks2*Tzew0)/(K0 + xp1(i) + Ks2);
end

plot(fp1, Tw10,'k-');
grid on;
hold on;
title('fp0 od Tw1');
xlabel ('fp0 ');
ylabel ('Tw1 ');
plot (fpN , Tw1N , 'ro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot (322) ;
Tzew1 =[-25:0.1:25];
Tw10=(((xp0*Tkz0 + Ks1*Tzew1)*(K0+xp0+Ks2))+ Ks2*Tzew1*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));
Tw20=(Tw10*(K0+xp0) + Ks2*Tzew1)/(K0 + xp0 + Ks2);

plot(Tzew1, Tw20,'k');
grid on;
hold on;
title('Tzew0 od Tw2');
xlabel ('Tzew0 ');
ylabel ('Tw2 ');
plot (TzewN , Tw2N , 'ro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot (324) ;
Tkz1 =[-15:0.1:35];
Tw10=(((xp0*Tkz1 + Ks1*Tzew0)*(K0+xp0+Ks2))+ Ks2*Tzew0*K0)/((xp0+Ks1+K0)*(K0+xp0+Ks2)-(K0^2+xp0*K0));
Tw20=(Tw10*(K0+xp0) + Ks2*Tzew0)/(K0 + xp0 + Ks2);

plot(Tkz1, Tw20,'k');
grid on;
hold on;
title('Tkz0 od Tw2');
xlabel ('Tkz0 ');
ylabel ('Tw2 ');
plot (TkzN , Tw2N , 'ro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot (326) ;
fp1 =[0:0.1:50];
xp1=cp*rop*fp1;

k=size(fp1,2);

for i=1:k
Tw10(i)=(((xp1(i)*Tkz0 + Ks1*Tzew0).*(K0+xp1(i)+Ks2))+ Ks2*Tzew0*K0)/((xp1(i)+Ks1+K0).*(K0+xp1(i)+Ks2)-(K0^2+xp1(i)*K0));
Tw20(i)=(Tw10(i)*(K0+xp1(i)) + Ks2*Tzew0)/(K0 + xp1(i) + Ks2);
end

plot(fp1, Tw20,'k-');
grid on;
hold on;
title('fp0 od Tw2');
xlabel ('fp0 ');
ylabel ('Tw2 ');
plot (fpN , Tw2N , 'ro');