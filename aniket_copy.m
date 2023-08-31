clc
%For Kv = 5 we define G_uncompensated as 
%G_uncompensated = G1 
num1 = [5];
den1 = [0.5 1.5 1 0] ; 
G1 = tf(num1 , den1 ) ;
%Extracting Gain Margin and Phase Margin 
[Gm, Pm, w_gm, w_pm] = margin(G1);
G1
%Epsilon value for adjustment if required 
epsilon = 5;

%calculating phie required 
phie_m = 40 - Pm + epsilon ;

%design of the lead compensator 
%D(s) = (s + 1/t )/ (s + 1/alpha * t)

alpha = (1- sin(deg2rad(phie_m)) ) / (1+ sin(deg2rad(phie_m)) );
up = 10 * log(1/alpha) ;
%By this much up value the bode plot of uncompensated will shift up 

%So Locating the freq where the uncomp bode plot hits -up value 
%Location has been done through stimulating bode plot 
% figure ;
% bode(G1) 
w_m = 5.43;
var =   w_m * sqrt(alpha) ;
T = 1/var;

%Compensator
num_C = [1 1/T];
den_C = [1 1/(alpha*T)] ;
G_C = tf(num_C , den_C) ;
figure 
margin(G_C)
G_C
%Finding Geq
G_eq = G1 * G_C;
G_eq
%For verfiying the phase margin and Gain margin , we draw the bode plot of
%G_eq
figure
margin(G_eq)
%alpha 

[Gm1, Pm1, w_gm1, w_pm1] = margin(G_eq);
% figure 
% margin(G_eq) 
 epsilon = 12;
%1.11 ,11.4
phie_req = 40 +epsilon  ;
 w_at_135 =  1.11;
 g_un_w = 11.4;
 beta = 10 ^(11.4/20);
wc_z = w_at_135/10 ;
 T1 = 1/wc_z;
wcp = 1/(beta * T1) ;
% 
% 
numc = [T1 1];
 denc = [beta*T 1] ;
 GC = tf(numc , denc) ;
GC
figure 
margin(GC)
 Geq = G_eq * GC;
 figure
 margin(Geq)