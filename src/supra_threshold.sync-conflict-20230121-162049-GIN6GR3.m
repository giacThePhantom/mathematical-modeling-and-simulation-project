% Utils
clear all;
close all;

% Useless?
dt = 0.0002;
sim_time = 60; %s
init_cond = [-65; 0];
global s_array;
global t_array;

[t, v] = ode15s(@neuron, [0 100], init_cond);


t_array = 0:100;

for i=1:length(t_array),
    s_array(end + 1) = compute_s(i);
end
    
subplot(2, 1, 1)
plot(t, v)
subplot(2, 1, 2)
plot(t_array, s_array)

% Functions
function dv = neuron(t, v)
    
    % Parameters
    E_L = -65; %mV
    g_L = 0.02; %uS
    C = 0.2; %nF
    V_th = -55; %mV
    alpha = 10000;
    gamma = 1;
    s = compute_s(t);
    %global s_array
    %global t_array
    disp([v(2) gamma])
    beta_test = exp(-(v(2)^2)/(2*gamma^2));
    disp(beta_test)

    %s_array(length(s_array) + 1) = s;
    %t_array(length(t_array) + 1) = t;

    dv = [
       1/C*(g_L*(E_L-v(1))+s) + alpha * (E_L-v(1)) * beta_test;
       1 - alpha * v(2) * H(v(1), V_th);
    ]; 
end

function s = compute_s(t)
    if t > 20 && t < 40
        s = 100;
    elseif t > 45 && t < 80
        s = 300; 
    else, s = 0;
    end
end

function heaviside = H(V, V_t)
    if V >= V_t
        heaviside = 1;
    else
        heaviside = 0;
    end
end

%figure;
%plot(T, Y(:, 1), 'LineWidth', 3.0);
%xlabel('Time [ms]');
%ylabel('Potential [mV]');
%grid on;


%for t = 1:length(T) - 1,
    %if t > 20 & t <=21,
        %s(t) = 10;
    %elseif t > 15000 & t< 200000 & mod(t, 20000) == 0,
        %s(t) = 5;
    %elseif t > 200000 & t< 250000 & mod(t, 20000) == 0,
        %s(t) = 3;
    %end
%end
%
%
%
%
%for t=1:length(T)-1,
    %T(t+1) = T(t) + dt * ( 1 - (alpha * T(t) * H(Vm(t), V_th)));
    %beta = exp(-(T(t).^2)/(2*gamma.^2));
    %Vm(t+1) = Vm(t) + dt * ((g_L*(E_L-Vm(t)) + s(t)*1000)/C + alpha*(E_L - Vm(t))*beta);
%end;
%
%T = T*1000;
%Vm
%
%subplot(2,1,1);
%plot(T,Vm,'LineWidth', 3.0);
%
%xlabel('Time [ms]');
%
%ylabel('Voltage [mV]');
%
%subplot(2,1,2)
%plot(T, s, 'r', 'LineWidth', 3.0);
%
%xlabel('Time [ms]');
%
%ylabel('Current [nA]');
%
%grid on;
%
%
%% Parameters from https://www.cns.nyu.edu/~david/handouts/integrate-and-fire.pdf
%function res = H(V, V_th)
    %if V >= V_th,
        %res = 1;
    %else
        %res = 0;
    %end
%end
