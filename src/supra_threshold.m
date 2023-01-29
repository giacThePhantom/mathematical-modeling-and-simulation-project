clear all;
close all;

opts = odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'InitialStep', 1e-8, 'MaxStep', 1e-2);
y_0 = [-73, 0];
[T, Y] = ode15s(@(t, y) neuron(t, y), [0 100], y_0, opts);

t_array = 0:100;
global s_array;

for i=1:length(t_array),
    s_array(end + 1) = s(i);
end
for t = 1:length(T),
    disp(t)
    disp(T(t))
    current(t) = s(T(t));
end

subplot(3,1,1)
plot(T, Y(:, 1), 'LineWidth', 3.0);
xlabel('Time [s]');
ylabel('Potential [mV]');
grid on;
subplot(3,1,2)
plot(T, Y(:, 2), 'g', 'LineWidth', 3.0);
xlabel('Time [s]');
ylabel('T');
grid on;
subplot(3,1,3)
plot(t_array, s_array, 'r', 'LineWidth', 3.0);
xlabel('Time [s]');
ylabel('Current');
grid on;

function dy = neuron(t, y)

    E_L = -73;
    g_L = 0.025;
    C = 0.375;
    V_th = -53;
    alpha_m = 10000;
    gamma_m = 1;
    V_r = -90;


    y = reshape(y, [], 2);
    V = y(:, 1);
    T = y(:, 2);

    beta_m = @(T) exp(-(T^2)/(2*gamma_m^2));

    dy = zeros(size(y));

    dy(:, 1) = (g_L*(E_L-V)+s(t))/C + alpha_m*(V_r-V)*beta_m(T);
    dy(:, 2) = 1 - alpha_m*T*H(V, V_th);
    dy = dy(:);
end

function heaviside = H(V, V_t)
    if V >= V_t
        heaviside = 1;
    else
        heaviside = 0;
    end
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
