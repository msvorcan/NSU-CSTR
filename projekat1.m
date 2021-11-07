clear 
close all
clc

kw = 10^-14;
Ca = 10^-6;
Cb = 10^-6;
V = 30;
Fa = 0.016667;
y_e = 7;


syms x1 x2 u
x1dot = Fa/V * (Ca - x1) - 1 / V* x1 * u;
x2dot = -Fa/V * x2 + 1/V * (Cb - x2) * u;
y = -log10(sqrt((x2-x1)^2 / 4 + kw) - (x2 - x1) / 2);
% [x1e,x2e] = solve([x1dot==0,x2dot==0],[x1,x2]);
% x1 = x1e;
% x2 = x2e;
% pom = eval(y);
% ue = solve(pom == y_e,u);
% ue = eval(ue);
% u = ue;
% x1e = eval(x1e);
% x2e = eval(x2e);
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

%% Teswtiranje sistema u otvorenoj sprezi


Ca = 10^-6;
for pom = 1:2 
    out = sim('projekat1_otvorena_sprega.slx');
    figure
    subplot(2,1,1)
    plot(out.t,out.u)
    subplot(2,1,2)
    plot(out.t,out.y)
    ylim([0 14])
    figure
    subplot(3,1,1)
    plot(out.t,out.x1)
    subplot(3,1,2)
    plot(out.t,out.x2)
    subplot(3,1,3)
    plot(out.x1,out.x2)
end

figure
hold all
for pom1 = 0.8:0.1:1.2 
    Ca = pom1*10^-6;
    out = sim('projekat1_otvorena_sprega.slx');
    plot(out.t,out.y)
    ylim([0 14])
end
 legend('80%Ca','90%Ca','Ca','110%Ca','120%Ca')

%% Projektovanje kontrolera na bazi inverzije dinamike
Ca = 10^-6;

sim('projekat1zatvorenasprega.slx');
figure
plot(y)

