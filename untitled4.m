

%%                      Advanced Control Systems
%                          Submitted By:
%                   Soham Karak, 200102107,Grp 16 
% Question 2:Design a PID controller for the syetem: 1/s(s+1)(s+5), obtain unit step response 
% with MATLAB/ Simulink, and compare the responses of systems.
% Clearing the workspace and command window

close all
clc
clear all

% Defining the plant Tf
num = [1];  % numerator coefficients
den = [1 6 5 0];  % denominator coefficients
G = tf(num, den);

% controlSystemDesigner(G)
%% Marginal stability point and Zeigler nichols parameters 

% % Plot the root locus
figure(1)
rlocus(G)
%Obtain the Kcritical value
gm = margin(G)
K  = gm;
num = [1];  % numerator coefficients
den = [1 6 5 30];  % denominator coefficients
G_ = tf(num, den);
step(G_,5)
% Obtain the period of sustained oscillations to get Pcr= 2.81
%% Zeigler Nichols PID tuning

G_Kp = tf(1,[1 6 5 K])
% step(G_Kp,10)
% hold on
Kc = 30;
Pc = 2.81

%Min overshoot PID model

Kp = 0.6*Kc
Ti = (0.5*Pc)
Td = (0.333*Pc)

Ki = Kp/Ti
Kd = Kp*Td
% Transfer Function of PID model
C = tf([Kd Kp Ki],[1 0])

%% PID Transfer Function obtained from root locus pole placement

s = tf('s');
C_updated = (4*(1 + 1.5*s)*(1 + 1.5*s))/s
pid(C_updated)

 %            1          
 % Kp + Ki * --- + Kd * s
 %            s          
 % with Kp = 13.2, Ki = 4, Kd = 10.9
 
%% Comparision of the step responses of both the PID controllers


G_final = feedback(C*G,1)
G_updated = feedback(C_updated*G,1)

% Transient time parametres of the PID models
S = stepinfo(G_final)
S1= stepinfo(G_updated)


% Plot the step responses
figure(2)
hold on
step(G_final) 
grid on
hold on
step(G_updated) 
grid on
hold on
step(feedback(G,1))
legend('Ziegler Nichols method','Root Locus method','Original closed loop response')
hold off

