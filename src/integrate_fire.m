clear all;
close all;
sim_time = 60; %s
dt = 0.0002;
T = 0:dt:sim_time; % 1 second simulation
Vm = zeros(sim_time/dt + 1, 1);
E_L = -73; %mV
Vm(1) = E_L; %mV
g_L = 0.025; %uS
C = 0.375; %nF
V_th = -53; %mV
V_R = -90;

s = zeros(sim_time/dt + 1, 1);


for t = 1:length(T) - 1,
    if t > 20 & t <=21,
        s(t) = 37.5;
    elseif t > 15000 & t< 200000 & mod(t, 20000) == 0,
        s(t) = 20;
    elseif t > 200000 & mod(t, 20000) == 0,
        s(t) = 5;
    end
end




for t=1:length(T)-1,
    if Vm(t) >= V_th,

        Vm(t+1) = V_R;

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
plot(T, s, 'r', 'LineWidth', 3.0);

xlabel('Time [ms]');

ylabel('Current [nA]');

grid on;


% Parameters from https://www.cns.nyu.edu/~david/handouts/integrate-and-fire.pdf
