clc;
clear;
%% data normal configuration
mass=14126; %=W
Ixx=3582;
Iyy=55802;
Izz=56669;
Ixz=2658;
b=21.95;

gravity=32.174;
alpha=2.30;
gamma=0;
%% table1
M=0.257;
W=14126;

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
%Assume
alpha0=2.3;
w0=V_tot*sin(alpha0);
u0=V_tot*cos(alpha0);
v0=0;
theta0=0;
Yr=0;
%%
% lateral state space
A_latr_full=[Yv     Yp      Yr-u0         -gravity*cos(theta0)      0
             Lb     Lp      Lr            0                         0
             Nb     Np      Nr            0                         0
             0      1       tan(theta0)   0                         0
             0      0       1/cos(theta0) 0                         0];
B_latr_full=[Yda     Ydr
             Lda     Ldr
             Nda     Ndr
             0       0
             0       0];
C_latr_full=eye(5);
D_latr_full=zeros(5,2);
latr_full_state_space=ss(A_latr_full,B_latr_full,C_latr_full,D_latr_full);
[V,D] = eig(A_latr_full);
[d,ind]=sort(diag(D));

disp('Full System Info');
damp(A_latr_full);
stepinfo(latr_full_state_space)

figure;
step(latr_full_state_space(:,1))
title('Full System Reponse to Step Input');

[yfss,t]=step(latr_full_state_space(:,1));
nondimenlat(yfss,t,b,u0);

figure;
impulse(latr_full_state_space(:,1))
title('Full System Reponse to Impulse Input');

[yfsi,t]=impulse(latr_full_state_space(:,1));
nondimenlat(yfsi,t,b,u0);

figure;
x0=[0 0 0 0 0]; %%%%%%%%%%%%%
initial(latr_full_state_space,x0)
title('Full System Response to Initial Conditions');


%%1 dof spiral
A_spiral=[Lr*Nb/Lb-Nr];
B_spiral=[0];
C_spiral=[1];
D_spiral=[0];

Spiral_state_space=ss(A_spiral,B_spiral,C_spiral,D_spiral);
[V_spiral,D_spiral] = eig(A_spiral);

disp('Spiral  System Info');
damp(A_spiral);
stepinfo(Spiral_state_space)


figure;
step(Spiral_state_space);
title('Spiral  Response to Step Input');



figure;
impulse(Spiral_state_space);
title('Spiral  Response to Impulse Input');

figure;
x0=[0]; %%%%%%%%%%%%%
initial(Spiral_state_space,x0)
title('Spiral Response to Initial Conditions');







% 2D dutch Roll
A_dutch=[Yv    Yr/u0-1;
         Nb      Nr];
B_dutch=[Yda Ydr/u0;
         Nda Ndr]; 
     C_dutch=eye(2);
D_dutch=zeros(2,2);
Dutch_roll_state_space=ss(A_dutch,B_dutch,C_dutch,D_dutch);
[V_dutch,D_dutch] = eig(A_dutch);
[d_dutch,ind_dutch]=sort(diag(D_dutch));

disp('Dutch Roll System Info');
damp(A_dutch);
stepinfo(Dutch_roll_state_space)

figure;
step(Dutch_roll_state_space(:,1));
title('Dutch Roll Response to Step Input');

figure;
impulse(Dutch_roll_state_space(:,1));
title('Dutch Roll Response to Impulse Input');

figure;
x0=[0;0]; %%%%%%%%%%%%%
initial(Dutch_roll_state_space,x0)
title('Dutch Response to Initial Conditions');






%%% %%%%%
%%%%%Rolling%%%%%%%%%
A_roll=[Lp];
B_roll=[Lda];
C_roll=[1];
D_roll=[0];

roll_state_space=ss(A_roll,B_roll,C_roll,D_roll);
[V_roll,D_roll] = eig(A_roll);

disp('Roll System Info');
damp(A_roll);
stepinfo(roll_state_space)

figure;
step(roll_state_space);
title('roll  Response to Step Input');

figure;
impulse(roll_state_space);
title('roll  Response to Impulse Input');

figure
x0=[0]; %%%%%%%%%%%%%
initial(roll_state_space,x0);
title('roll  Response to Initial Conditions');
     
     
     
     
     