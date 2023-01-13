clear all;
close all;
sim_time = 60; %s
dt = 0.0002;
T = 0:dt:sim_time; % 1 second simulation
Vm = zeros(sim_time/dt + 1, 1);
E_L = -65; %mV
Vm(1) = -65; %mV
g_L = 0.02; %uS
C = 0.2; %nF
V_th = -55; %mV
alpha = 10000;
gamma = 1;


s = zeros(sim_time/dt + 1, 1);


opts = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
y_0 = [-0.1, -0.4, 15];
[T, Y] = ode15s(@(t, y) neuron(t, y), [0 100], y_0, opts);

figure;
plot(T, Y(:, 1), 'LineWidth', 3.0);
xlabel('Time [ms]');
ylabel('Potential [mV]');
grid on;



function dy = neuron(t, y)

	y = reshape(y, [], 3);
	xi = y(:, 1);
	xi_rest = y(:, 2);
	tau = y(:, 3);

	dy = zeros(size(y));
	dy(:,1) = (xi_rest - xi)/tau;
	dy = dy(:);

end



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
