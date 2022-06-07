%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      This is a program to simulate PD Control to control
%      a pendulum with double masses manipulator.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main
clear;      % clear all variables
%clf;       % clear figure windows
clc;        % clear console screen

global m1 m2 g l1 l2 Kp1 Kd1 Kp2 Kd2 M E U
global ypos1 ypos2 yvel1 yvel2

% system parameters
m1= 2; % kg
m2 = 1; % kg
l1= 1; % m
l2= 1; % m
g=9.81;
E = [1 -1;0 1];

% PD controller parameters
Kp1 = 200;
Kd1 = 200;
Kp2 = 200;
Kd2 = 200;

% Inverse Kinematics for a desired trajectory
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


TH1_MLF{2} = matlabFunction(S.theta_1(2),'Vars',[L_1 L_2 XE YE]);
TH2_MLF{2} = matlabFunction(S.theta_2(2),'Vars',[L_1 L_2 XE YE]);

% initial/terminal conditions
xmat = 0.5;
ymat = -1;
y0_1= TH1_MLF{2}(L1,L2,xmat,ymat); % rad
y0_2= TH2_MLF{2}(L1,L2,xmat,ymat); % rad
ylog=[];
y=[y0_1 y0_2 0 0 0 0]';  % [theta1 theta2 dot_theta1 dot_theta2 input_torque1 input_torque2]
dt=0.001;
t_fin=2;


% The simulation starts right here.
for t=0:dt:t_fin
    % Trajectory
    xmat = 0.5;
    ymat = -1*cos(pi*t/2);
    
    % desired values 
    ypos1= TH1_MLF{2}(L1,L2,xmat,ymat);  % theta1           
    ypos2= TH2_MLF{2}(L1,L2,xmat,ymat);  % theta2
    
    % Inverse Kinematics 함수 에서 각도가 180도가 넘어가면 음의 각으로 변하는 것을 방지
    if ypos1 <= 0
        ypos1 = 2*pi+ypos1;
    end
    
    
    if t == 0
        ypos1_before = y0_1;
        ypos2_before = y0_2;
    end
    yvel1= (ypos1-ypos1_before)/dt;      % theta_dot_1
    yvel2= (ypos2-ypos2_before)/dt;      % theta_dot_2
    
    ypos1_before = ypos1;
    ypos2_before = ypos2;
   
    % Control Algorithm here!!

    % Here, y(5), y(6) are not a state, but is used as the input torque.
    % The control torques are applied at the system, simulated with
    % the 4/5th-order RK-method ('ode45') for the duration of dt.
    [t_dummy,ydata]=ode45(@manipulator, [0 dt], y);
    % Since ode45 takes adaptive stepping, each call of this produces so many data.
    % To store only the last data, which is associated with time 'dt' 
    % from the time of each call, the last data is extracted.
    y=ydata(end,:);
    % For the real torque values
    Tor = (E)\M*U;
    y(5) = Tor(1);
    y(6) = Tor(2);
    
%     
%     if y(5) > 10000
%         y(5) = y(5) * 10^(-8);
%     end
%     if y(6) > 10000
%         y(6) = y(6) * 10^(-4);
%     end

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
title_str=sprintf('Trajectory1');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,2);
plot(ylog(:,1),ylog(:,3), ylog(:,1), ylog(:,9));
xlabel('Time [s]', 'FontSize', 14);
ylabel('\theta(t) [rad]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Trajectory2');
title(title_str, 'FontSize', 16);
grid;


subplot(2,4,3);
plot(ylog(:,8),ylog(:,10));
xlabel('e (=\theta-\theta_d) [rad]', 'FontSize', 14);
ylabel('de/dt [rad/s]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Phase Plot1');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,4);
plot(ylog(:,9),ylog(:,11));
xlabel('e (=\theta-\theta_d) [rad]', 'FontSize', 14);
ylabel('de/dt [rad/s]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Phase Plot2');
title(title_str, 'FontSize', 16);
grid;


subplot(2,4,5);
plot(ylog(:,1),ylog(:,12));
xlabel('Time [s]', 'FontSize', 14);
ylabel('e (=\theta-\theta_d) [rad]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Error1');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,6);
plot(ylog(:,1),ylog(:,13));
xlabel('Time [s]', 'FontSize', 14);
ylabel('e (=\theta-\theta_d) [rad]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Error2');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,7);
plot(ylog(:,1),ylog(:,6));
xlabel('Time [s]', 'FontSize', 14);
ylabel('\tau (t) [N-m]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Torque1');
title(title_str, 'FontSize', 16);
grid;

subplot(2,4,8);
plot(ylog(:,1),ylog(:,7));
xlabel('Time [s]', 'FontSize', 14);
ylabel('\tau (t) [N-m]', 'FontSize', 14);
%title_str=sprintf('Nonlinear: torque = %2.1f \\times \\tau_0', (per+100)/100)
title_str=sprintf('Torque2');
title(title_str, 'FontSize', 16);
grid;

beep;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part describes the dynamics of the system to simulate.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot= manipulator(t, x)
global m1 m2 g l1 l2 Kp1 Kd1 Kp2 Kd2 M U
global ypos1 ypos2 yvel1 yvel2
xd(1)=x(3);
xd(2)=x(4);
M = [(m1+m2)*l1^2 m2*l1*l2*cos(x(1)-x(2));m2*l1*l2*cos(x(1)-x(2)) m2*l2^2];
C = [0 m2*l1*l2*sin(x(1)-x(2))*x(4)^2;0 -1*m2*l1*l2*sin(x(1)-x(2))*x(3)^2];
G = [(m1+m2)*l1*sin(x(1))*g; m2*l2*g*sin(x(2))];
U = [Kp1*(ypos1-x(1))+Kd1*(yvel1-x(3));Kp2*(ypos2-x(2))+Kd2*(yvel2-x(4))];
xd(5)=0;            % this is input u(t); therefore, its derivative is meaningless.
xd(6)=0;            % this is input u(t); therefore, its derivative is meaningless.
A = (M)\(-C-G)+U;
xd(3) = A(1);
xd(4) = A(2);
xdot=xd';
end
