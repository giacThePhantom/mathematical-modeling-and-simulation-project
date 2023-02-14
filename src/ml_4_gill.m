function [V,M,N,t,tau_arr, eps_arr]=ml_4_gill(csv_name, tau_arr_in, eps_arr_in, show_figures)

    params = readtable(csv_name);

    global Mcalcium_tot Mcalcium
    global Npotassium_tot Npotassium
    global vK vL vCa
    global gK gL gCa
    global C
    global tau
    global I

    % Parameters
    tmax = params.tmax;

    Npotassium_tot = params.nktot;
    Mcalcium_tot = params.mcatot;

    phi_m = params.phi_m;
    phi_n = params.phi_n;
    va = params.va;
    vb = params.vb;
    vc = params.vc;
    vd = params.vd;

    vK = params.vk;
    vL = params.vl;
    vCa = params.vca;
    gK = params.gk;
    gL = params.gl;
    gCa = params.gca;
    C = params.c;

    V0 = params.v0; % start at an arbitrary middle voltage
    M0 = params.m0; % start with calcium channels all closed
    N0 = params.n0; % start with half potassium channels open

    I = params.i;

    t0=0;

    % Initial conditions
    t=t0; % global "external" time
    V=V0; % initial voltage
    M=M0; % use M to record the time course
    N=N0; % use N to record the time course

    Mcalcium=M0; % initial state of calcium channel
    Npotassium=N0; % initial state of potassium channel

    % Morris Lecar
    global Iapp; Iapp=@(t)I; % applied current
    global xi_m; xi_m=@(v)(v-va)/vb; % scaled argument for m-gate input voltage
    global minf; minf=@(v)0.5*(1+tanh(xi_m(v))); % m-gate activation function
    global tau_m; tau_m=@(v)1./(phi_m*cosh(xi_m(v)/2)); % m-gate time constant
    global alpha_m; alpha_m=@(v)(minf(v)./tau_m(v));
    global beta_m; beta_m=@(v)((1-minf(v))./tau_m(v));
    global xi_n; xi_n=@(v)(v-vc)/vd; % scaled argument for n-gate input
    global ninf; ninf=@(v)0.5*(1+tanh(xi_n(v))); % n-gate activation function
    global tau_n; tau_n=@(v)1./(phi_n*cosh(xi_n(v)/2));  % n-gate time constant
    global alpha_n; alpha_n=@(v)(ninf(v)./tau_n(v));
    global beta_n; beta_n=@(v)((1-ninf(v))./tau_n(v));

    options=odeset('Events', @nextevent);

    index = 0;
    tau_arr = tau_arr_in;
    eps_arr = eps_arr_in;

    while t(end)<tmax
        index = index + 1;
        if length(tau_arr) < index
            tau_arr = [tau_arr, -log(rand)];
        end

        if length(eps_arr) < index
            eps_arr = [eps_arr, rand];
        end

        tau = tau_arr(index);
        epsilon = eps_arr(index);

        %disp(epsilon)
        U0=[V0;0;0;0;0;M0;N0];
        tspan=[t(end),tmax];
        [tout,Uout,~,~,event_idx]=ode23(@dudtfunc,tspan,U0,options);
        Vout=Uout(:,1); % voltage at time of next event
        Mout=Uout(:,6); % number of calcium channels open at end of next event
        Nout=Uout(:,7); % number of potassium channels open at end of next event
        t=[t,tout'];
        V=[V,Vout'];
        M=[M,Mout'];
        N=[N,Nout'];


        lambda_0 = (alpha_n(V(end)) * (Npotassium_tot-N0)) + (beta_n(V(end)) * N0) + (alpha_m(V(end)) * (Mcalcium_tot-M0)) + (beta_m(V(end)) * M0);
        lambda_1 = alpha_n(V(end)) * (Npotassium_tot-N0);
        lambda_2 = lambda_1 + (beta_n(V(end)) * (N0));
        lambda_3 = lambda_2 + (alpha_m(V(end)) * (Mcalcium_tot-M0));

        q_1 = lambda_1/lambda_0;
        q_2 = lambda_2/lambda_0;
        q_3 = lambda_3/lambda_0;

        if (epsilon < q_1)
            N0=N0+1;
        elseif (epsilon >= q_1 && epsilon < q_2)
            N0=N0-1;   % increment channel state
        elseif (epsilon >= q_2 && epsilon < q_3)
            M0 = M0 + 1;
        else
            M0 = M0 - 1;
        end

        Mcalcium=M0;
        Npotassium=N0;
        if M0>Mcalcium_tot, error('M>Mcalcium_tot'), end
        if M0<0, error('M<0'), end
        if N0>Npotassium_tot, error('N>Npotassium_tot'), end
        if N0<0, error('N<0'), end

        V0=V(end);

    end % while t(end)<tmax

    if show_figures == 1
        % Plot output
        figure
            subplot(5,1,1),
            plot(t,V, 'LineWidth', 3),
            xlabel('Time', 'FontSize', 16),
            ylabel('V', 'FontSize', 16)
            subplot(5,1,2),
            plot(t,M, 'r', 'LineWidth', 3),
            ylabel('M', 'FontSize', 16),
            xlabel('Time', 'FontSize', 16)
            subplot(5,1,3),
            plot(t,N, 'g', 'LineWidth', 3),
            ylabel('N', 'FontSize', 16),
            xlabel('Time', 'FontSize', 16)
            subplot(5,1,4),
            plot(V,M,'.-', 'color', 'r', 'LineWidth', 2),
            xlabel('V', 'FontSize', 16), ylabel('M', 'FontSize', 16)
            subplot(5,1,5),
            plot(V,N,'.-', 'color', 'g', 'LineWidth', 2),
            xlabel('V', 'FontSize', 16),
            ylabel('N', 'FontSize', 16)
            grid on
    end

end

function dudt=dudtfunc(t,u)

    global Iapp % applied current
    global Mcalcium Mcalcium_tot
    global Npotassium Npotassium_tot
    global alpha_m beta_m
    global alpha_n beta_n
    global vK vL vCa
    global gK gL gCa
    global C

    v=u(1); % extract the voltage from the input vect
dudt=[(Iapp(t)-gCa*(Mcalcium/Mcalcium_tot)*(v-vCa)-gL*(v-vL)-gK*(Npotassium/Npotassium_tot)*(v-vK))/C;
       alpha_n(v)*(Npotassium_tot-Npotassium);     % Potassium chan. opening, internal time elapsed
       beta_n(v)*Npotassium; % Potassium chan. closing, internal time elapsed
       alpha_m(v)*(Mcalcium_tot-Mcalcium); % Calcium chan. opening, internal time elapsed
       beta_m(v)*Mcalcium; % Calcium chan. closing, internal time elapsed
        0; % M is constant between events
        0]; % N is constant between events
end

% Define behavior at threshold crossing
function [value,isterminal,direction] = nextevent(~,u)

    global tau

    value=[u(2)+u(3)+u(4)+u(5)-tau];
    isterminal=[1]; % stop and restart integration at crossing
    direction=[1]; % increasing value of the quantity at the trigger
end
