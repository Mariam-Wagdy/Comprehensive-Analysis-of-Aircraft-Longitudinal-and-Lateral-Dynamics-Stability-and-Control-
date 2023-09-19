clc;
clear;
%% data normal configuration

Ixx=3582;
Iyy=55802;
Izz=56669;
Ixz=2658;
c=9.55; 

gravity=32.174;
alpha=2.30;
gamma=0;
%% table1
M=0.257;
W=14126;
%% table2
Xu=-0.737;
Zu=-0.204;
Mu=0.000294;
Xw=0.0631;
Zw=-0.570;
Mw=0.00732;
Zwd=0;Zq=0;
Mwd=-0.000304;
Mq=-0.317;
Xde=1.19;%Xds
Zde=-29.7;
Mde=-4.79;
Xdth=0.00228;
Zdth=0.994e-4;%
Mdth=0;

%% table
Yv=-0.178;
Yb=-51.1; %assume Yp=Yb
Yp=-51.1;
Lb=-20.9;
Nb=2.68;
Lp=-1.38;
Np=0.0993;
Lr=1.16;
Nr=-0.157;
Yda=0;
Lda=4.76;
Nda=0.266;
Ydr=0.0317;
Ldr=5.35;
Ndr=-0.923;
%eq
V_tot=287;
% Lv=Lp/V_tot;
% Nv=Np/V_tot;
%Assume
alpha0=2.3;
w0=V_tot*sin(alpha0);
u0=V_tot*cos(alpha0);
theta0=0;
Yr=0;
%%  longitudinal
A_long_full=[Xu                     Xw                    -w0                        -gravity*cos(theta0)
             Zu/(1-Zwd)             Zw/(1-Zwd)            (Zq+u0)/(1-Zwd)            (-gravity*sin(theta0))/(1-Zwd)
             Mu+(Mwd*Zu)/(1-Zwd)    Mw+(Mwd*Zw)/(1-Zwd)   Mq+(Mwd*(Zu+u0))/(1-Zwd)   (-Mwd*gravity*sin(theta0))/(1-Zwd)
             0                      0                     1                           0];
B_long_full=[Xde                         Xdth
             Zde/(1-Zwd)                 Zdth/(1-Zwd)
             Mde+(Mwd*Zde)/(1-Zwd)       Mdth+(Mwd*Zdth)/(1-Zwd)
             0                           0];
C_long_full=eye(4);
D_long_full=zeros(4,2);
long_full_state_space= ss(A_long_full,B_long_full,C_long_full,D_long_full);
[V_long_full,D_long_full] = eig(A_long_full);
[d_long_full,ind_long_full]=sort(diag(D_long_full));

disp('Full System step Info');
damp(A_long_full);

figure;
step(long_full_state_space(:,1))
title('Full System Response to a Step Input');

[yfss,t]=step(long_full_state_space(:,1));

stepinfo(yfss,t);
nondimen(yfss,t,c,u0);

figure;
impulse(long_full_state_space(:,1))
title('Full System Response to an Impulse Input');

[yfsi,t]=impulse(long_full_state_space(:,1));
nondimen(yfsi,t,c,u0);
stepinfo(yfsi,t);

x0=[u0 ;0; 0 ;0];  %%%%%%%%%%

figure;
initial(long_full_state_space,x0)
title('Full System Response to Initial Inputs');
%%%%%%%%%%% formalized

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% long period
A_long_period=[Xu             -gravity
               -Zu/(Zq+u0)     0];
B_long_period=[Xde             Xdth
              -Zde/(Zq+u0)     -Zdth/(Zq+u0)];
C_long_period=eye(2);
D_long_period=zeros(2,2);

long_period_state_space=ss(A_long_period,B_long_period,C_long_period,D_long_period);
[V_long_period,D_long_period] = eig(A_long_period);
[d_long_period,ind_long_period]=sort(diag(D_long_period));

disp('Long Period Appromiation System Info');
damp(A_long_period);
stepinfo(long_period_state_space(:,1))

figure;
step(long_period_state_space(:,1));
title('Long Period Response to Step Input');


figure;
impulse(long_period_state_space(:,1));
title('Long Period Response to Impulse Input');


x0=[0;0]; %%%%%%%%%%%%%
figure;
initial(long_period_state_space,x0)
title('Long Period Response to Initial Conditions');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% short period
A_short_period=[Zw/(1-Zwd)           (Zq+u0)/(1-Zwd)
                Mw+(Mwd*Zw)/(1-Zwd)   Mq+(Mwd*(Zu+u0))/(1-Zwd)];
B_short_period=[ Zde/(1-Zwd)          Zdth/(1-Zwd)
                 Mde+(Mwd*Zde)/(1-Zwd) Mdth+(Mwd*Zdth)/(1-Zwd)];
C_short_period=eye(2);
D_short_period=zeros(2,2);
short_period_state_space=ss(A_short_period,B_short_period,C_short_period,D_short_period);
[V_short_period,D_short_period] = eig(A_short_period);
[d_short_period,ind_short_period]=sort(diag(D_short_period));

disp('Short Period Appromiation System Info');
damp(A_short_period);
S_short_period=stepinfo(short_period_state_space)

figure;
step(short_period_state_space(:,1));
title('Short Period Response to Step Input');

figure;
impulse(short_period_state_space(:,1));
title('Short Period Response to Impulse Input');

figure;
x0=[u0;theta0]; %%%%%%%%%%%%%
initial(short_period_state_space,x0)
title('Short Period Response to Initial Conditions');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% transfer function
%long full
 long_full_tf=tf(long_full_state_space);
 long_full_u_de=      long_full_tf(1,1);
 long_full_w_de=      long_full_tf(2,1);
 long_full_alpha_de=  long_full_w_de/u0;
 long_full_q_de=      long_full_tf(3,1);
 long_full_theta_de=  long_full_tf(4,1);
 long_full_u_dT=      long_full_tf(1,2);
 long_full_w_dT=      long_full_tf(2,2);
 long_full_alpha_dT=  long_full_w_dT/u0;
 long_full_q_dT=      long_full_tf(3,2);
 long_full_theta_dT=  long_full_tf(4,2);

 %long period
 long_period_tf=tf(long_period_state_space);
 long_period_u_de=    long_period_tf(1,1);
 long_period_theta_de=long_period_tf(2,1);
 long_period_u_dT=    long_period_tf(1,2);
 long_period_theta_dT=long_period_tf(2,2);

 %short period  
 short_period_tf=tf(short_period_state_space);
 short_period_u_de=    short_period_tf(1,1);
 short_period_theta_de=short_period_tf(2,1);
 short_period_u_dT=    short_period_tf(1,2);
 short_period_theta_dT=short_period_tf(2,2);
 
