clear all;
close all;

opts = odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'InitialStep', 1e-8, 'MaxStep', 1e-2);
y_0 = [-60, 0, 0];
[T, Y] = ode15s(@(t, y) neuron(t, y), [0 1000], y_0, opts);

subplot(3,1,1)
plot(T, Y(:, 1), 'LineWidth', 3.0);
xlabel('Time [ms]', 'FontSize', 16);
ylabel('Potential [mV]', 'FontSize', 16);

subplot(3,1,2)
plot(T, Y(:, 2), 'g', 'LineWidth', 3.0);
xlabel('Time [ms]', 'FontSize', 16);
ylabel('N', 'FontSize', 16 );

grid on;

subplot(3,1,3)
plot(T, Y(:, 3), 'r', 'LineWidth', 3.0);
xlabel('Time [ms]', 'FontSize', 16);
ylabel('M', 'FontSize', 16 );

grid on;

function dy = neuron(t, y)

    v_k = -84;
    v_l = -60;
    v_ca = 120;
    i_app = 100;
    g_k = 8;
    g_l = 2;
    C = 20;
    v_a = -1.2;
    v_b = 18;
    v_c = 2;
    v_d = 30;
    phi_n = 0.04;
    phi_m = 0.4;
    g_ca = 4.4;

    y = reshape(y, [], 3);
    V = y(:, 1);
    N = y(:, 2);
    M = y(:, 3);

    % beta_m = @(T) exp(-(T^2)/(2*gamma_m^2));
    eps_n = (V-v_c)/v_d;
    eps_m = (V-v_a)/v_b;
    %alpha = (phi*cosh(eps))/(1+exp(2*eps));
    %beta = (phi*cosh(eps))/(1+exp(-2*eps));
    n_inf = (1 + tanh(eps_n))/2;
    tau_n = 1/(phi_n*cosh(eps_n/2));
    m_inf = (1+tanh(eps_m))/2;
    tau_m = 1/(phi_m*cosh(eps_m/2));

    dy = zeros(size(y));

    dy(:, 1) = (i_app-g_ca*M*(V-v_ca) - g_l*(V-v_l) - g_k*N*(V-v_k))/C;
    dy(:, 2) = (n_inf - N)/tau_n;
    dy(:, 3) = (m_inf - M)/tau_m;
    dy = dy(:);
end
