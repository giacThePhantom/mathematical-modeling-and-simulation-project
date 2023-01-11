
dt = 0.0002;
T = 0:dt:1; % 1 second simulation
Vm = zeros(1/dt, 1);
E_L = -75;
Vm(1) = E_L;
g_L = 10;
C = 10;
V_th = -55;

s = zeros(1/dt, 1);


for t=1:length(T)-1,

		if t>20 & t < 200,
			s(t) = 10000;
		end

    if Vm(t) > V_th,

        Vm(t+1) = E_L;

    else,

      Vm(t+1) = Vm(t) + dt * (g_L*(E_L-Vm(t)) + s(t))/C;
    end;

end;

T = T*1000;
figure;
plot(T,Vm,'LineWidth', 3.0);

xlabel('Time [ms]');

ylabel('Voltage [mV]');

grid on;
