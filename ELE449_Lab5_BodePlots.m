%%% Transfer function
close all;
% Make 3 bode plots
%input volt to out
D = 0.5;
Vo = 12;
C = 470e-6;
Rc = 45e-3;
R = 0.8;
L = 23e-6;
Q = R*sqrt(C/L);
Vs = 24;
Rl = 50e-3;

%Gvs(s) = Vo(s)/vs(s);
%Gvs(s) = D*(1+(s/w_esr))/(1+(s/Qwo)+(s^2/wo^2));

%Gvd(s) = Vo(s)
wesr = 1/C*Rc
wo = sqrt(1/L*C)
wz = Rl/L

num1 = [D/wesr D];
den1 = [1/wo^2 1/(Q*wo) 1];
Gvs = tf(num1, den1)
figure;

bodeplot (Gvs)
title(['Gvs'])
grid on
hold

num2 = [Vs/wesr Vs];
den2 = [1/wo^2 1/(Q*wo) 1];
Gvd = tf(num2, den2);
figure;
bodeplot (Gvd);
title(['Gvd'])
grid on
hold

num3 = [Rl/wz 1/wesr Rl];
den3 = [1/wo^2 1/(Q*wo) 1];
Zp = tf(num3, den3);
figure;
bodeplot (Zp)
title(['Zp'])
grid on
hold


