%%                  Advanced Control System Term Project
%                           Submitted by:
%                         Soham Karak, 200102107
%                          Group 16, Question 1

%% Program to design the compensator for a given system (Cacscaded Lag compensator)
% Design a Compensator for a unity feedback system with
% open loop transfer function G(s)= K/(s(s+1)(0.5s +1)) to satisfy the
% following specifications:
% (i) velocity error constant K_v=5, 
% (ii) Phase margin =40 degrees. 
% (iii)Gain Margin = 10dB


% Clear the workspace and the command window

close all
clear all
clc

% Specifying the following requirement variables
PM_desired = 40;
GM_desired = 10;
Kv=1;


%Designing the uncompensated sytem plant Tf
disp(['Uncompensated system Transfer function']);
% Calculate the value of gain K from steady state value e_ss and use it as the uncompensaated system tf gain. 
num=[5];
den=[0.5 1.5 1 0];
G=tf(num,den)

% Plotting the Bode plot of the gain adjusted uncompensated system 
[gm, pm, wgc, wpc] = margin(G)
figure(1);
bode(G)
title({"Uncompensated system with gain adjusted Kv=5sec^-1",sprintf('Phase Margin = %0.4f deg, Gain Margin = %0.4f dB',pm, gm)})
grid on;

% Finding the desired value of phase by adding a error tolerance value.
% The frequency at which the desired phase value was obtained(wValue), is marked
% and the corresponding value of magnitude(mag_dB) for the frequency is obtained.
% We then design a lag compensator such that the compensator brings down
% the magintude plot by mag_db such that the phase crossover frequency
% becomes wvalue. The compensator poles and zeros are wpc/10 and
% wpc/10*beta

required_PM=-180+PM_desired+10
[mag, phase,w] = bode(G);
phase_dB = squeeze(phase);
wValue = interp1(phase_dB, w, required_PM, 'spline');
disp(['The frequency value corresponding to phase margin = ' num2str(required_PM) ' degrees is ' num2str(wValue) ' rad/s']);
mag_dB = 20*log10(squeeze(mag));
w_at = interp1(w, mag_dB, wValue, 'spline');
disp(['The gain at frequency = ' num2str(wValue) ' degrees is ' num2str(w_at) ' dB']);

% Lag compensator design
% The pole is palced at 0.1 wgc so that the phase shift added by the lag
% compensator would be negligible

beta =10^(w_at/20)
zero=wValue/10;
pole=zero/beta;
T1=10/wValue
T2=beta*T1
disp(['Zero of the lead compensator is at s = ' num2str(zero)]);
disp(['Pole of the lead compensator is at s = ' num2str(pole)]);

%% Design of 1st Lag compensator

%Creating the lag compensator transfer function and the overall closed loop
%tf to see whether we have satisfied the design requirements
disp(['Compensator Transfer function']);
Gc=tf([T1 1],[T2 1])
[gm, pm, wgc, wpc] = margin(Gc);
figure;
bode(Gc);
title({'Compensator',sprintf('Phase Margin = %0.4f deg, Gain Margin = %0.4f dB', pm, gm)});
grid on;
final_TF1 = G*Gc

% Plotting the Bode plot of the Partially-compensated system 

[gm, pm, wgc, wpc] = margin(final_TF1)
figure(2);
bode(final_TF1)
margin(final_TF1)
title({'Compensated System',sprintf('Phase Margin = %0.4f deg, Gain Margin = %0.4f dB', pm, gm)});
grid on;
[mag, phase,w] = bode(final_TF1);

%We observe that although the steady state and phase margin requirement is
%fulfilled the gain margin is still inadequate so we plan to cascade
%another lag compensator to add teh remaining gain shift required
%% Design of the 2nd Lag Compensator
% We have choosen that the lag compensator would shift the mag plot by
% another 5 dbs so that the gain margin requirement wouyld be satisfied .
%For that puprose we have taken the addition gain added by the lag compensator to be 8.5 dB's at wgc/10 or the zero  of teh first compensator as well as any frequency there after  

beta  = 10^(8/20)
wz = wgc/10
t = 1/wz
num = [t 1]
den = [beta*t 1]

% Design of 2nd lag compenstaor tf
Gc1 = tf(num,den)
bode(Gc1)
margin(Gc1)
[gm, pm, wgc, wpc] = margin(Gc1)
%% Final Transfer function
final_TF = G*Gc*Gc1

% Plotting the Bode plot of the compensated system 
[gm, pm, wgc, wpc] = margin(final_TF)

figure(3);
bode(final_TF)
margin(final_TF)
title({'Compensated System',sprintf('Phase Margin = %0.4f deg, Gain Margin = %0.4f dB', pm, gm)});
grid on;
[mag, phase,w] = bode(final_TF);
%% Step Response

cl_G=feedback(G,1);     
% Closed loop TF of compensated system
cl_Final_TF=feedback(final_TF,1);

figure(4);
hold on;
subplot(211);step(cl_G);title({'Step response of uncompensated system(closed loop system)'});
grid on;
subplot(212);step(cl_Final_TF);title({'Step response of compensated system(closed loop system)'});
grid on;
hold off;
%% Ramp response

% the ramp signal is created for a duration of 0-10secs
t=0:0.1:10;
alpha=1;
ramp=alpha*t;      
[y1,t]=lsim(cl_G,ramp,t);
[y2,t]=lsim(cl_Final_TF,ramp,t);

figure(5);

subplot(211);
hold on
plot(t,t)
plot(t,y1);title({'Ramp response of uncompensated system(closed loop system)'});
legend('Ramp signal','Uncompensated system')
grid on;

subplot(212); 
hold on 
plot(t,t)
plot(t,y2);title({'Ramp response of compensated system(closed loop system)'});
legend('Ramp signal','compensated system')
grid on;
hold off;