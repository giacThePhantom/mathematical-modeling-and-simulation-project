% Utils
clear all;
close all;

dt = 0.0002;
sim_time = 60; %s
p_ampa = 0.875;
p_gaba = 0.0625;
p_nmda = 0.0625;
opts = odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'InitialStep', 1e-8, 'MaxStep', 1e-2);

y_0 = [-73; 0; 0; 0; 0; 0];
global s_array;
global t_array;

[T, Y] = ode15s(@(t, y) neuron(t, y), [0 100], y_0, opts);

t_array = 0:100;

for i=1:length(t_array),
    s_array(end + 1) = s(i);
end

subplot(7, 1, 1)
plot(T, Y(:, 1), 'LineWidth', 3.0)
legend("Voltage")


subplot(7, 1, 2)
plot(T, Y(:, 2), 'g', 'LineWidth', 3.0)
legend("InterSpikeTime")

subplot(7, 1, 3)
plot(T, Y(:, 3), 'm', 'LineWidth', 3.0)
legend("SK")

subplot(7, 1, 4)
plot(T, Y(:, 4), 'y', 'LineWidth', 3.0)
legend("AMPA")

subplot(7, 1, 5)
plot(T, Y(:, 5), 'c', 'LineWidth', 3.0)
legend("GABA")

subplot(7, 1, 6)
plot(T, Y(:, 6), 'k', 'LineWidth', 3.0)
legend("NMDA")

subplot(7, 1, 7)
plot(t_array, s_array, 'r', 'LineWidth', 3.0)
legend("Current")

% Functions
function dy = neuron(t, y)

    % Parameters
    E_L = -73; %mV
    V_R = -90;
    g_L = 0.025; %nS
    C = 0.375; %nF
    V_th = -53; %mV
    alpha_m = 10000;
    gamma_m = 1;
    beta_m = @(T) exp(-(T^2)/(2*gamma_m^2));

    g_sk = 0.128; %nS

    g_ampa = 0.024; %nS
    g_gaba = 0.064; %nS
    g_nmda = 0.08; %nS

    E_sk = -90; %mV
    E_ampa = 0; %mV
    E_gaba = -70; %mV
    E_nmda = 0; %mV

    p_ampa = 0.875;
    p_gaba = 0.0625;
    p_nmda = 0.0625;

    tau_ampa = 2.4; %ms
    tau_gaba = 7; %ms
    tau_nmda = 100; %ms
    tau_sk = 80; %ms

    a = -53;
    b = 100;

    y = reshape(y, [], 6);
    V = y(:, 1);
    T = y(:, 2);
    X_sk = y(:, 3);
    X_ampa = y(:, 4);
    X_gaba = y(:, 5);
    X_nmda = y(:, 6);

   dy = zeros(size(y));
   dy(:, 1) = 1/C*(g_L*(E_L - V) + g_sk * X_sk *(E_sk - V) + g_ampa * X_ampa * (E_ampa - V) + g_gaba * X_gaba * (E_gaba - V) + (g_nmda* X_nmda * (E_nmda - V))/(1+exp(-((V - a)/(b))))) + alpha_m * (V_R-V) * beta_m(T);
   dy(:, 2) = 1 - alpha_m * T * H(V, V_th);
   dy(:, 3) = ((1 - X_sk) * 4 * beta_m(T) - X_sk)/(tau_sk);
   dy(:, 4) = ((1 - X_ampa) * (p_ampa + s(t)) - X_ampa)/(tau_ampa);
   dy(:, 5) = ((1 - X_gaba) * p_gaba - X_gaba)/(tau_gaba);
   dy(:, 6) = ((1 - X_nmda) * p_nmda - X_nmda)/(tau_nmda);

   dy = dy(:);

end

function compute_s = s(t)
    if t > 10 && t < 10.1
        compute_s = 63770;
    elseif t > 15 && t < 25
        compute_s = 3;
    elseif t > 30 && t < 32
        compute_s = 3;
    elseif t > 36 && t < 38
        compute_s = 3;
    elseif t > 42 && t < 45
        compute_s = 3;
    elseif t > 50 && t < 52
        compute_s = 1;
    elseif t > 54 && t < 56
        compute_s = 1;
    elseif t > 58 && t < 60
        compute_s = 1;
    else
        compute_s = 0;
    end
end

function heaviside = H(V, V_t)
    if V >= V_t
        heaviside = 1;
    else
        heaviside = 0;
    end
end
