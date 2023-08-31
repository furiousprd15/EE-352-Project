%% Program to design the compensator for a given system 
% Design a Compensator for a unity feedback system with
% open loop transfer function G(s)= K/(s(s+1)(o.5s +1)) to satisfy the
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

%% Lead compensator Design
Pmd=40;                       % Desired Phase Margin (PM)
Phi_m = Pmd - Pm +10  
if Phi_m>60   % check if angle needed is greater than 60
   phm=Phi_m/2
else
   phm=Phi_m
end
% Maximum phase lead angle (degree) with tolerance of 5 degrees 
Phi_mr=phm*(pi/180)        % Maximum phase lead angle (radian)
alpha=(1-sin(Phi_mr))/(1+sin(Phi_mr))
Mc=20*log10(alpha)          % Magnitude to be compensated in db
%% Input: 2D array of x and y values
[mag,phase,wout] = bode(G);


% % Input: Value of magnitude for which phase needs to be found
% mag_val = Mc;
% 
% % Find the index of the closest value of magnitude
% [~, idx] = min(abs(mag - mag_val));
% 
% % Extract the corresponding value of phase
% phase_val = phase(idx);
% 
% % Print the value of phase
% disp(['The value of phase for magnitude = ', num2str(mag_val), ' is ', num2str(phase_val), ' radians.']);

 %% 
 
% Locate the frequency in Figure(1) for Mc
wm=3.01;
T = 1/(wm*sqrt(alpha))
z=1/T;                         % Zero of lead compensator
p=1/(T*alpha);                  % Pole of lead compensator
gain=alpha;
numc=[T 1];
denc=[T*alpha 1];
Gc=tf(numc,denc)
figure(2)
bode(Gc)
% Total forward transfer function of the compensated system
Gt=Gc*Gc*G

% Comparison of compensated and uncompensated bode plots
figure(3)
bode(G,'--r', Gt,'-'), grid on

legend('Uncompensated system', 'Compensated system')
title('Comparison of compensated and uncompensated bode plots')
[Gm,Pm,Wcg,Wcp] = margin(Gt)
%% Since H(s)=1, the feedback transfer function 
Gc1u=feedback(G/5,1);         % Closed loop TF of uncompensated system
Gclc=feedback(Gt,1);        % Closed loop TF of compensated system
 
% Comparison of compensated and uncompensated step responses
figure(3)
subplot(2,1,1)
step(Gc1u, '--r'); grid on  
title('Uncompensated system step response')

subplot(2,1,2)
step(Gclc, '-'); grid on  
title('Compensated system step response')

%% Comparison of compensated and uncompensated ramp responses
t=0:0.1:10;
figure(4)
a = 1;
ramp = a*t;
[y1,t] = lsim(Gc1u,ramp,t)
hold on
plot(t,y1,'-'), grid on
hold on
plot(t,t)
hold on
[y2,t] = lsim(Gclc,ramp,t)
plot(t,y2,'-'), grid on
legend('Uncompensated system ramp response','Ramp','Compensated system ramp response')
hold off



