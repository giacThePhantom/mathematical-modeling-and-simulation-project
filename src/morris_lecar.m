clear all;
close all;

opts = odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'InitialStep', 1e-8, 'MaxStep', 1e-2);
y_0 = [-84, 0];
[T, Y] = ode15s(@(t, y) neuron(t, y), [0 100], y_0, opts);

subplot(2,1,1)
plot(T, Y(:, 1), 'LineWidth', 3.0);
xlabel('Time [ms]');
ylabel('Potential [mV]');

subplot(2,1,2)
plot(T, Y(:, 2), 'LineWidth', 3.0);
xlabel('Time [ms]');
ylabel('N');

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
    phi = 0.04;
    g_ca = 4.4;

    y = reshape(y, [], 2);
    V = y(:, 1);
    N = y(:, 2);

    % beta_m = @(T) exp(-(T^2)/(2*gamma_m^2));
    eps = (V-v_c)/v_d;
    %alpha = (phi*cosh(eps))/(1+exp(2*eps));
    %beta = (phi*cosh(eps))/(1+exp(-2*eps));
    n_inf = (1 + tanh(eps))/2;
    tau = 1/(phi*cosh(eps/2));
    m_inf = (1+tanh((V-v_a)/v_b))/2;

    dy = zeros(size(y));

    dy(:, 1) = (i_app-g_ca*m_inf*(V-v_ca) - g_l*(V-v_l) - g_k*N*(V-v_k))/C;
    dy(:, 2) = (n_inf - N)/tau;
    dy = dy(:);
end

function compute_s = s(t)
    if t > 10 && t < 10.1
        compute_s = 6;
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
