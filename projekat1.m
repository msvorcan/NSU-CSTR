clear 
close all
clc

%inicijalizacija konstanti sistema

kw = 10^-14;
Ca = 10^-6;
Cb = 10^-6;
V = 30;
Fa = 0.016667;
y_e = 7;

%% Odredjivanje ravnoteznih stanja


syms x1 x2 u
x1dot = Fa/V * (Ca - x1) - 1 / V* x1 * u;
x2dot = -Fa/V * x2 + 1/V * (Cb - x2) * u;
y = -log10(sqrt((x2-x1)^2 / 4 + kw) - (x2 - x1) / 2);

[x1e,x2e,ue] = solve([x1dot==0,x2dot==0,y==y_e],[x1,x2,u]);
x1e = eval(x1e);
x2e = eval(x2e);
uep = eval (ue);

%% Provera stabilnosti rs , i linearizacija

syms x1 x2 u
x1dot = Fa/V * (Ca - x1) - 1 / V* x1 * u;
x2dot = -Fa/V * x2 + 1/V * (Cb - x2) * u;
y = -log10(sqrt((x2-x1)^2 / 4 + kw) - (x2 - x1) / 2);
A = [diff(x1dot,x1), diff(x1dot,x2); diff(x2dot,x1) , diff(x2dot,x2)];
B = [diff(x1dot,u); diff(x2dot,u)];
C = [diff(y,x1) , diff(y,x2)];
e = eig(A);
u = uep;
x1 = x1e;
x2 = x2e;
epom = eval(e);
A1 = eval(A);
B1 = eval(B);
C1 = eval(C);
[b,a] = ss2tf(A1,B1,C1,0);
G = tf(b,a);                          
figure
margin(G);
figure
bode(G);

%% Testiranje sistema u otvorenoj sprezi


Ca = 10^-6;
for pom = 1:2 
    out = sim('projekat1_otvorena_sprega.slx');
    figure
    subplot(2,1,1)
    plot(out.t,out.u)
    xlabel('vreme[s]')
    ylabel('protok baze[1/s]')
    title('upravljanje')
    subplot(2,1,2)
    plot(out.t,out.y)
    ylim([0 14])
    xlabel('vreme[s]')
    ylabel('PH vrednost')
    title('izlaz sistema')
    figure
    subplot(3,1,1)
    plot(out.t,out.x1)
    xlabel('vreme[s]')
    ylabel('koncentracija kiseline[mol/l]')
    title('promenljiva stanja x1')
    subplot(3,1,2)
    plot(out.t,out.x2)
    xlabel('vreme[s]')
    ylabel('koncentracija baze[mol/l]')
    title('promenljiva stanja x2')
    subplot(3,1,3)
    plot(out.x1,out.x2)
    xlabel('koncentracija kiseline[mol/l]')
    ylabel('koncentracija baze[mol/l]')
    title('fazorska ravan')
end
%%
pom = 1;
figure(1)
hold all
for pom1 = 0.8:0.1:1.2 
    Ca = pom1*10^-6;
    out = sim('projekat1_otvorena_sprega.slx');
    subplot(2,1,1)
    hold all
    plot(out.t,out.u)
    xlabel('vreme[s]')
    ylabel('protok baze[1/s]')
    title('upravljanje')
    subplot(2,1,2)
    hold all
    plot(out.t,out.y)
    ylim([0 14])
    xlabel('vreme[s]')
    ylabel('PH vrednost')
    title('izlaz sistema')
    disp(mean(out.y))
end
 legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')
 figure(2)
 for pom1 = 0.8:0.1:1.2 
    Ca = pom1*10^-6;
    out = sim('projekat1_otvorena_sprega.slx');
    subplot(3,1,1)
    hold all
    plot(out.t,out.x1)
    xlabel('vreme[s]')
    ylabel('koncentracija kiseline[mol/l]')
    title('promenljiva stanja x1')
    subplot(3,1,2)
    hold all
    plot(out.t,out.x2)
    xlabel('vreme[s]')
    ylabel('koncentracija baze[mol/l]')
    title('promenljiva stanja x2')
    subplot(3,1,3)
    hold all
    plot(out.x1,out.x2)
    xlabel('koncentracija kiseline[mol/l]')
    ylabel('koncentracija baze[mol/l]')
    title('fazorska ravan')
 end
legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')

%% Projektovanje kontrolera na bazi inverzije dinamike
Ca = 10^-6;

K = 0.01;
Ti =  900.09;
Tt =  0.01*Ti;

sim('projekat1zatvorenasprega.slx');
figure
plot(tout, y)
ylim([6 8])
figure
plot(tout, u)

%% Testiranje robusnosti PI kontrolera na bazi inverzije dinamike
K = 0.01;
Ti =  900.09;
Tt =  0.01*Ti;

figure
hold on
koncentracije = [0.8*Ca 0.9*Ca Ca 1.1*Ca 1.2*Ca];

for i = 1 : 5
    Ca = koncentracije(i);
    sim('projekat1zatvorenasprega.slx');
    plot(tout, y)
end
legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')
Ca = 10^-6;
figure
hold on
koncentracije = [0.8*Ca 0.9*Ca Ca 1.1*Ca 1.2*Ca];

for i = 1 : 5
    Ca = koncentracije(i);
    sim('projekat1zatvorenasprega.slx');
    plot(tout, u)
end
legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')
Ca = 10^-6;
%% Testiranje otpornosti na sum

K = 0.01;
Tt =  0.01*Ti;
figure
sim('projekat1zatvorenasprega.slx');
plot(tout, y)
ylim([6 8])
figure

plot(tout, u)

%% PID regulacija

out = sim('projekat1_PID_init.slx');
% figure
% plot(out.t,out.yd_ZN)
% figure
% plot(out.tout, out.y_ZN)
% hold all
du = uep;
dy = out.y_ZN(end) -7;
pomtau = out.y_ZN > (0.1*dy + 7);
pomT = out.y_ZN > (0.63*dy + 7);
k = find(pomtau,1,'first');
p = find(pomT,1,'first');
k = k-1000000;
p = p-1000000;
tau = k / 1000;
T = (p-k)/1000;
% [N, Npos] = max(out.yd_ZN(1000001 : end));
% xa = 0 : 0.001 : 5000;
% ya = out.y_ZN(1000000 + Npos) + N * (xa - (1000000 + Npos) * 0.001);
% plot(xa, ya)
% 
% pom = ones(1, length(xa));
% plot(xa, 7*pom)
% plot(xa, out.y_ZN(end)*pom)
% % xlim([49500 51000])
% % ylim([6.9 7.15])
% Tp = 1841.3865 - 1000;
% tau = 0.001;
% 
%xlim([9000 15000])
mi=0.4;


%% Testiranje PI-a
Kp = mi*T*du/dy/tau
Ti = T
Td = T/8
Tf = Td/20
Tt = 0.01*Ti

out = sim('PID_nova_metoda.slx');
figure
plot(tout, y)
ylim([6 8])
figure
plot(tout, u)
