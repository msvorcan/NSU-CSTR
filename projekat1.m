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

%% Provera stabilnosti R.S. i linearizacija

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
    sim('projekat1zatvorenasprega.slx');
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
    
%% Provera robusnosti sistema za parametar Ca u otvorenoj sprezi
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
    ylim([6.5 7.5])
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

%% Grafici
    %out = sim('projekat1zatvorenasprega.slx');
    figure
    subplot(2,1,1)
    plot(t, u)
    xlabel('vreme[s]')
    ylabel('protok baze[1/s]')
    title('upravljanje')
    %xlim([0 1000])
    subplot(2,1,2)
    plot(t, y)
    ylim([6 8])
    %xlim([0 1000])
    xlabel('vreme[s]')
    ylabel('PH vrednost')
    title('izlaz sistema')
    figure
    subplot(3,1,1)
    plot(t,x1)
    xlabel('vreme[s]')
    ylabel('koncentracija kiseline[mol/l]')
    title('promenljiva stanja x1')
    subplot(3,1,2)
    plot(t, x2)
    xlabel('vreme[s]')
    ylabel('koncentracija baze[mol/l]')
    title('promenljiva stanja x2')
    subplot(3,1,3)
    plot(x1, x2)
    xlabel('koncentracija kiseline[mol/l]')
    ylabel('koncentracija baze[mol/l]')
    title('fazorska ravan')

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
xlabel('t[s]')
ylabel('pH vrednost')
title('Testiranje robusnosti izlaza')

Ca = 10^-6;
figure
hold on
koncentracije = [0.8*Ca 0.9*Ca Ca 1.1*Ca 1.2*Ca];

for i = 1 : 5
    Ca = koncentracije(i);
    sim('projekat1zatvorenasprega.slx');
    plot(tout, u)
end

xlabel('t[s]')
ylabel('protok baze[l/s]')
title('Testiranje robusnosti upravljanja')

legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')
Ca = 10^-6;

for pom1 = 0.8:0.1:1.2 
    Ca = pom1*10^-6;
    out1 = sim('projekat1zatvorenasprega.slx');
    subplot(3,1,1)
    hold all
    plot(t, x1)
    xlabel('vreme[s]')
    ylabel('koncentracija kiseline[mol/l]')
    title('promenljiva stanja x1')
    subplot(3,1,2)
    hold all
    plot(t, x2)
    xlabel('vreme[s]')
    ylabel('koncentracija baze[mol/l]')
    title('promenljiva stanja x2')
    subplot(3,1,3)
    hold all
    plot(x1, x2)
    xlabel('koncentracija kiseline[mol/l]')
    ylabel('koncentracija baze[mol/l]')
    title('fazorska ravan')
 end
legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')
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

mi=0.4;


%% Testiranje PI-a
Kp = mi*T*du/dy/tau
Ti = T
Td = T/8
Tf = Td/20
Tt = 0.01*Ti

out = sim('PID_nova_metoda.slx');
figure
plot(t, y)
ylim([6 8])
figure
plot(t, u)

out = sim('PID_nova_metoda.slx');
    figure
    subplot(2,1,1)
    plot(t, u)
    xlabel('vreme[s]')
    ylabel('protok baze[1/s]')
    title('upravljanje')
    %xlim([0 1000])
    subplot(2,1,2)
    plot(t, y)
    ylim([6 8])
    %xlim([0 1000])
    xlabel('vreme[s]')
    ylabel('PH vrednost')
    title('izlaz sistema')
    figure
    subplot(3,1,1)
    plot(t,x1)
    xlabel('vreme[s]')
    ylabel('koncentracija kiseline[mol/l]')
    title('promenljiva stanja x1')
    subplot(3,1,2)
    plot(t, x2)
    xlabel('vreme[s]')
    ylabel('koncentracija baze[mol/l]')
    title('promenljiva stanja x2')
    subplot(3,1,3)
    plot(x1, x2)
    xlabel('koncentracija kiseline[mol/l]')
    ylabel('koncentracija baze[mol/l]')
    title('fazorska ravan')

%% Testiranje robusnosti PI kontrolera na bazi odskocnog odziva

figure
hold on
koncentracije = [0.8*Ca 0.9*Ca Ca 1.1*Ca 1.2*Ca];

for i = 1 : 5
    Ca = koncentracije(i);
    sim('PID_nova_metoda.slx');
    plot(tout, y)
end
legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')
xlabel('t[s]')
ylabel('pH vrednost')
title('Testiranje robusnosti izlaza')

Ca = 10^-6;
figure
hold on
koncentracije = [0.8*Ca 0.9*Ca Ca 1.1*Ca 1.2*Ca];

for i = 1 : 5
    Ca = koncentracije(i);
    sim('PID_nova_metoda.slx');
    plot(tout, u)
end

xlabel('t[s]')
ylabel('protok baze[l/s]')
title('Testiranje robusnosti upravljanja')

legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')
Ca = 10^-6;

for pom1 = 0.8:0.1:1.2 
    Ca = pom1*10^-6;
    out1 = sim('PID_nova_metoda.slx');
    subplot(3,1,1)
    hold all
    plot(t, x1)
    xlabel('vreme[s]')
    ylabel('koncentracija kiseline[mol/l]')
    title('promenljiva stanja x1')
    subplot(3,1,2)
    hold all
    plot(t, x2)
    xlabel('vreme[s]')
    ylabel('koncentracija baze[mol/l]')
    title('promenljiva stanja x2')
    subplot(3,1,3)
    hold all
    plot(x1, x2)
    xlabel('koncentracija kiseline[mol/l]')
    ylabel('koncentracija baze[mol/l]')
    title('fazorska ravan')
 end
legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')
