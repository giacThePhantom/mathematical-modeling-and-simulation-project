% Utils
clear all;
close all;

% Useless?
dt = 0.0002;
sim_time = 60; %s
p_ampa = 0.875; 
p_gaba = 0.0625;
p_nmda = 0.0625;

init_cond = [-65; 0; 0; p_ampa/(1+p_ampa); p_gaba/(1+p_gaba); p_nmda/(1+p_nmda)];
global s_array;
global t_array;

[t, v] = ode15s(@neuron, [0 100], init_cond);

t_array = 0:100;

for i=1:length(t_array),
    s_array(end + 1) = compute_s(i);
end
    
subplot(2, 1, 1)
plot(t, v)
legend("Voltage", "InterSpikeTime", "SK", "AMPA", "GABA", "NMDA")
subplot(2, 1, 2)
plot(t_array, s_array)

% Functions
function dv = neuron(t, v)
    
    % Parameters
    E_L = -73; %mV
    V_R = -90;
    g_L = 25; %nS
    C = 0.375; %nF
    V_th = -53; %mV
    alpha = 0;
    gamma = 1;
    s = compute_s(t);
    beta_test = exp(-(v(2)^2)/(2*gamma^2));

    g_sk = 128; %nS

    g_ampa = 24; %nS
    g_gaba = 64; %nS
    g_nmda = 8; %nS

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

    dv = [
       1/C*(g_L*(E_L - v(1)) + g_sk * v(3) *(E_sk - v(1)) + g_ampa * v(4) * (E_ampa - v(1)) + g_gaba * v(5) * (E_gaba - v(1)) + (g_nmda* v(6) * (E_nmda - v(1)))/(1+exp(-((v(1) - a)/(b))))) + alpha * (V_R-v(1)) * beta_test;
       -alpha * v(2) * H(v(1), V_th);
       ((1 - v(3)) * 4 * beta_test - v(3))/(tau_sk);
       ((1 - v(4)) * (p_ampa + s) - v(4))/(tau_ampa);
       ((1 - v(5)) * p_gaba - v(5))/(tau_gaba);
       ((1 - v(6)) * p_nmda - v(6))/(tau_nmda);      
    ]; 
end

function s = compute_s(t)
    if t > 20 && t < 40
        s = 0;
    elseif t > 45 && t < 80
        s = 0; 
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
