%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This is a program to simulate Slind-Mode Control (SMC)
%      and Integral Sliding-Mode Control (ISMC) to control
%      a single-mass manipulator (or inverted pendulum).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main
clear;  % clear all variables
%clf;        % clear figure windows
clc;        % clear console screen

global m1 m2 g l1 l2 

RAD2DEG=180/pi;
g=9.81;

% system parameters
m1= 2; % kg
m2 = 1; % kg
l1= 1; % m
l2= 1; % m

% SMC parameters
Cbar = [100 0;0 100];
epsilon = 0.05;
Mhat = [(m1+m2)*l1^2 m2*l1*l2; m2*l1*l2 m2*l2^2];


% trajectory parameters
A= pi/3; % amplitude of sinusoidal input
w= 2*pi*1; % frequency of sinusoidal inupt (w=2*pi*f)


% Inverse Kinematics
syms L_1 L_2 theta_1 theta_2 XE YE
L1 = 1;
L2 = 1;
XE_RHS = L_1*sin(theta_1) + L_2*sin(theta_2);
YE_RHS = -1*L_1*cos(theta_1)-1*L_2*cos(theta_2);
XE_MLF = matlabFunction(XE_RHS,'Vars',[L_1 L_2 theta_1 theta_2]);
YE_MLF = matlabFunction(YE_RHS,'Vars',[L_1 L_2 theta_1 theta_2]);

XE_EQ = XE == XE_RHS;
YE_EQ = YE == YE_RHS;

S = solve([XE_EQ YE_EQ], [theta_1 theta_2]);
simplify(S.theta_1);
simplify(S.theta_2);

TH1_MLF{1} = matlabFunction(S.theta_1(1),'Vars',[L_1 L_2 XE YE]);
TH1_MLF{2} = matlabFunction(S.theta_1(2),'Vars',[L_1 L_2 XE YE]);
TH2_MLF{1} = matlabFunction(S.theta_2(1),'Vars',[L_1 L_2 XE YE]);
TH2_MLF{2} = matlabFunction(S.theta_2(2),'Vars',[L_1 L_2 XE YE]);
% initial/terminal conditions
xmat = 0.5;
ymat = -1;
y0_1= TH1_MLF{2}(L1,L2,xmat,ymat); % rad
y0_2= TH2_MLF{2}(L1,L2,xmat,ymat); % rad
ylog=[];
y=[y0_1 y0_2 0 0 0 0]';  % [theta1 theta2 dot_theta1 dot_theta2 input_torque1 2]
dt=0.001;
t_fin=1.08;


% The simulation starts right here.
for t=0:dt:t_fin
    % Inverse Kinematics
    xmat = 0.5;
    ymat = -1*cos(pi*t/2);
    
    % desired values 
    ypos1= TH1_MLF{2}(L1,L2,xmat,ymat);  % theta1           
    ypos2= TH2_MLF{2}(L1,L2,xmat,ymat);  % theta2
    if ypos1 <= 0
        ypos1 = 2*pi+ypos1;
    end
    if t == 0
        ypos1_before = y0_1;
        ypos2_before = y0_2;
        yvel1_before = 0;
        yvel2_before = 0;
    end
    yvel1= (ypos1-ypos1_before)/dt;                        % theta_dot_1
    yvel2= (ypos2-ypos2_before)/dt;                        % theta_dot_2
    
    ypos1_before = ypos1;
    ypos2_before = ypos2;
    
    yacc1= (yvel2-yvel1_before)/dt;              % theta_dot_dot_1
    yacc2= (yvel2-yvel2_before)/dt;              % theta_dot_dot_2
    
    yvel1_before = yvel1;
    yvel2_before = yvel2;
    
    err1= y(1)-ypos1;
    err2= y(2)-ypos2;
    err_vel1 = y(3)- yvel1;
    err_vel2 = y(4)- yvel2;
    
    % Control Algorithm here!!
    s= [Cbar(1,1)*(err1)+err_vel1;Cbar(2,2)*(err2)+err_vel2]; % Slidng surface defined
    Chat = [0 m2*l1*l2*abs(y(4)); m2*l1*l2*abs(y(3)) 0];
    K1 = Mhat(1,1)*abs(yacc1-Cbar(1,1)*err_vel1)+Chat(1,1)*abs(y(3)-s(1))+...
        Mhat(1,2)*abs(yacc2-Cbar(2,2)*err_vel2)+Chat(1,2)*abs(y(4)-s(2))+epsilon;
    K2 = Mhat(2,1)*abs(yacc1-Cbar(1,1)*err_vel1)+Chat(2,1)*abs(y(3)-s(1))+...
        Mhat(2,2)*abs(yacc2-Cbar(2,2)*err_vel2)+Chat(2,2)*abs(y(4)-s(2))+epsilon;
    G_theta = [(m1+m2)*l1*sin(y(1))*g; m2*l2*g*sin(y(2))];
    sgns = [sign(s(1)) sign(s(2))];
    u= [G_theta(1)-K1*sgns(1);G_theta(2)-K2*sgns(2)];

    % when the motor power is off, u=0
    % u=0;
    % Here, y(3) is not a state, but is used as the input torque.
    y(5)= u(1);
    y(6)= u(2);
    % The control torque is applied at the system, simulated with
    % the 4/5th-order RK-method ('ode45') for the duration of dt.
    [t_dummy,ydata]=ode45(@manipulator, [0 dt], y);
    % Since ode45 takes adaptive stepping, each call of this produces so many data.
    % To store only the last data, which is associated with time 'dt' 
    % from the time of each call, the last data is extracted.
    y=ydata(end,:);
   
    %  storing data for later processing such as plotting
    ylog=[ylog; t+dt y ypos1 ypos2 yvel1 yvel2 y(1)-ypos1 y(2)-ypos2 y(3)-yvel1 y(4)-yvel2];
    % y=[
    %1   time
    %234567 y = [position1 position2 vel1 vel2 input-torque1 input_torque2]^T 
    %8   desired position1
    %9   desired position1
    %10   desired velocity
    %11   desired velocity 
    %12,13   position error1,2 
    %14,15   velocity error1,2
end


figure(1)
subplot(2,4,1);
plot(ylog(:,1),ylog(:,2), ylog(:,1), ylog(:,8));
xlabel('Time [s]', 'FontSize', 14);
ylabel('\theta(t) [rad]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('SMC: Trajectory \\lambda= %2.1f; K= %4.2f', K1);
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,2);
plot(ylog(:,1),ylog(:,3), ylog(:,1), ylog(:,9));
xlabel('Time [s]', 'FontSize', 14);
ylabel('\theta(t) [rad]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('SMC: Trajectory \\lambda= %2.1f; K= %4.2f', K2);
title(title_str, 'FontSize', 16);
grid;


subplot(2,4,3);
plot(ylog(:,8),ylog(:,10));
xlabel('e (=\theta-\theta_d) [rad]', 'FontSize', 14);
ylabel('de/dt [rad/s]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Phase Plot');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,4);
plot(ylog(:,9),ylog(:,11));
xlabel('e (=\theta-\theta_d) [rad]', 'FontSize', 14);
ylabel('de/dt [rad/s]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Phase Plot');
title(title_str, 'FontSize', 16);
grid;


subplot(2,4,5);
plot(ylog(:,1),ylog(:,12));
xlabel('Time [s]', 'FontSize', 14);
ylabel('e (=\theta-\theta_d) [rad]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Error');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,6);
plot(ylog(:,1),ylog(:,13));
xlabel('Time [s]', 'FontSize', 14);
ylabel('e (=\theta-\theta_d) [rad]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Error');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,7);
plot(ylog(:,1),ylog(:,6));
xlabel('Time [s]', 'FontSize', 14);
ylabel('\tau (t) [N-m]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Torque');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,8);
plot(ylog(:,1),ylog(:,7));
xlabel('Time [s]', 'FontSize', 14);
ylabel('\tau (t) [N-m]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Torque');
title(title_str, 'FontSize', 16);
grid;

beep;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part describes the dynamics of the system to simulate, 
% which is expressed by  
%
%                     dx/dt = f(t, x).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot= manipulator(t, x)
%
% Dynamics is described here, which is described by
% dx/dt = f(x, u, t), where x denotes the state.
%
% In this code, the input torque is treated as a state!
% THus, x= (x(t), v(t), tau(t))^T
%
% THe dynamics to simulate is:
%      ml^2\ddot{z}=-mgl*sin(theta) + tau 
%
% NOTE that input torque is x(3).
% 
global m1 m2 g l1 l2 
xd(1)=x(3);
xd(2)=x(4);
M = [(m1+m2)*l1^2 m2*l1*l2*cos(x(1)-x(2));m2*l1*l2*cos(x(1)-x(2)) m2*l2^2];
C = [0 m2*l1*l2*sin(x(1)-x(2))*x(4)^2;0 -1*m2*l1*l2*sin(x(1)-x(2))*x(3)^2];
G = [(m1+m2)*l1*sin(x(1))*g; m2*l2*g*sin(x(2))];
E = [1 -1;0 1];
xd(5)=0;
xd(6)=0;% this is input u(t); therefore, its derivative is meaningless.
A = (M)\(-C-G+E*[x(5);x(6)]);
xd(3) = A(1);
xd(4) = A(2);
xdot=xd';
end
