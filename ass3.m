close all
clear all
clc
% Calculate the value of gain K from staedy state value e_ss and use it. 
num=[5];
den=[0.5 1.5 1 0];
G=tf(num,den)
C = tf([10 1],[195 1]);
T = C*G
f = feedback(T,1);
bode(T)
[A,b,c,d] = margin(T)