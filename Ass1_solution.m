%% Program to design the compensator for a given system 
% Design a Compensator for a unity feedback system with
% open loop transfer function G(s)= K/(s(s+1)(0.5s +1)) to satisfy the
% following specifications:
% (i) velocity error constant K_v=5, 
% (ii) Phase margin =40 degrees. 
% (iii)Gain Margin = 10dB
 
% Name of the Student:Soham Karak
% Roll number:200102107
 
close all
clear all
clc
% Calculate the value of gain K from staedy state value e_ss and use it. 
num=[5];
den=[0.5 1.5 1 0];
G=tf(num,den)
%% Bode plot of the uncompensated system 
figure(1)
% To check PM and GM of the uncompensated system with gain adjusted           
%          5
% ---------------------
% 0.5 s^3 + 1.5 s^2 + s

title('Bode plot of uncompensated system')
[Gm,Pm,Wcg,Wcp] = margin(G)
bode(G), grid on 
margin(G)
%% Phase margin and System Parameters calculations
% PM desired = 40
% Tolerance error value = 5(would be tuned later depending on scenario)
PM  = 43;  % this refers to the value above -180 degree line
w_m = 0.56;
Mag  = 17.5;
beta = 10^(Mag/20)

%% Compensator design
gamma = 10;
w_cz = w_m/gamma;
T = 1/w_cz;
p = 1/(beta*T);
num = [T 1];
den = [beta*T 1];
C = tf(num,den)
bode(C)
Gt = G*C
bode(G,'--r', Gt,'-'), grid on
% 
legend('Uncompensated system', 'Compensated system')
title('Comparison of compensated and uncompensated bode plots')
[Gm1,Pm1,Wcg1,Wcp1] = margin(Gt)

%% Lead compensator design

alpha  = 1/beta;
phi = asin(1-alpha/1+alpha)



