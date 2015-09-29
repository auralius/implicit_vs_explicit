% Auralius Manurung, 2015
% A system: mass-spring-damper, with initial condition
% We will compare the explicit vs implicit method
% while using ode45 executed at 1kHz as the ground truth.

function test1()
clear all;
close all;
clc;

% Simulation time
t_start = 0;
t_end   = 10;  %final time in seconds.
Ts = 0.1;
time_span =t_start:Ts:t_end;

% This are the paramters for the mass-spring-damper system
k = 100;       % spring stiffness. N/m
m = 5;         % mass, kg
%b = 2*sqrt(k*m);  %critical damping
b = 0;         % When no damping, implit Euler still damps the oscillation!
 
% Inital conditoin
initial_position = 0.1;
initial_speed    = 0; 
x0 = [initial_position  initial_speed];

figure
hold on

% Use ode45 as ground truth (1kHz sampling rate):
[t,x]=ode45(@msd, t_start:0.001:t_end, x0, [], m, b, k);
plot(t,x(:,1), 'r');

[x,t] = implicit_euler(time_span, x0, m, b, k);
plot(t,x(:,1), '--b');

[x,t] = explicit_euler(time_span, x0, m, b, k);
plot(t,x(:,1), '--m');

legend('ode45', 'implicit euler', 'explicit euler')
end
 
%**************************************
% solves m x''+ c x' + k x = 0
% use ode45
%**************************************
function xdot=msd(t,x, m, b, k)    
    xdot_1 = x(2);
    xdot_2 = -(b/m)*x(2) - (k/m)*x(1);

    xdot = [xdot_1 ; xdot_2 ];
end

%*************************************
% Implicit Euler
% By Auralius Manurung
% Reference: D. Baraff and A. Witkin, Large steps in cloth simulation
%            SIGGRAPH ’98
%*************************************
function [x,t] = implicit_euler(time_span, x0, m, b, k)
Ts = time_span(2)-time_span(1);

x(1, :) = x0;
    
dfdx = -k;
dfdv = -b;

x_ = x0(1);
v_ = x0(2);
f_ = -k*x0(1) + -b*x0(2);
for i = 1:length(time_span)-1
    deltav = (Ts/m*(f_+Ts*dfdx*v_)) / (1-Ts/m*dfdv-Ts^2/m*dfdx);
    deltax = Ts*(v_ + deltav);
    x_ = deltax + x_;
    v_ = deltav + v_;
    f_ = -k*x_ + -b*v_;
    
    x(i+1, :) = [x_ v_];
end    

t = time_span;
end

%*************************************
% Explicit Euler
% By Auralius Manurung
%*************************************
function [x,t] = explicit_euler(time_span, x0, m, b, k) 
Ts = time_span(2)-time_span(1);
 
x(1, :) = x0;
    
x_ = x0(1);
v_ = x0(2);
f_ = -k*x0(1) + -b*x0(2);
for i = 1:length(time_span)-1
    deltav = Ts*(f_/m);
    deltax = Ts*(v_ + deltav);
    x_ = deltax + x_;
    v_ = deltav + v_;
    f_ = -k*x_ + -b*v_;
    
    x(i+1, :) = [x_ v_];
end    

t = time_span;
end