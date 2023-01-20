clear all;
close all;
opts = odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'InitialStep', 1e-8, 'MaxStep', 1e-2);
y_0 = [-73, 10];
[T, Y] = ode15s(@(t, y) supra_threshold(t, y), [0 100], y_0, opts);

subplot(2,1,1)
plot(T, Y(:, 1), 'LineWidth', 3.0);
xlabel('Time [ms]');
ylabel('Potential [mV]');
subplot(2,1,2)
plot(T, Y(:, 2), 'LineWidth', 3.0);
xlabel('Time [ms]');
ylabel('T');
grid on;




function dy = supra_threshold(t, y)

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

    H = 1/(1 + exp(-200*(V-V_th)));

    s = 1 * (sign(t-10) - sign(t-11));

    dy = zeros(size(y));

    dy(:, 2) = 1- alpha_m*T*H;
    dy(:, 1) = (g_L*(E_L-V)+s)/C + alpha_m*(V_r-V)*beta_m(T);
    dy = dy(:);
end
