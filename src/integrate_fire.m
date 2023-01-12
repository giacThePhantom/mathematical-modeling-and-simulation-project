clear all;
close all;
dt = 0.0002;
T = 0:dt:60; % 1 second simulation
Vm = zeros(60/dt + 1, 1);
E_L = -65; %mV
Vm(1) = -65; %mV
g_L = 0.02; %uS
C = 0.2; %nF
V_th = -55; %mV

s = zeros(60/dt + 1, 1);


for t=1:length(T)-1,

		if t>0 & t <= 1,
		    s(t) = 10;
    end

    if Vm(t) >= V_th,

        Vm(t+1) = E_L;

    else,

      Vm(t+1) = Vm(t) + dt * (g_L*(E_L-Vm(t)) + s(t)*1000)/C;
    end;

end;

T = T*1000;

subplot(2,1,1);
plot(T,Vm,'LineWidth', 3.0);

xlabel('Time [ms]');

ylabel('Voltage [mV]');

subplot(2,1,2)
plot(T, s, 'LineWidth', 3.0);

xlabel('Time [ms]');

ylabel('Current [nA]');

grid on;


% Parameters from https://www.cns.nyu.edu/~david/handouts/integrate-and-fire.pdf
