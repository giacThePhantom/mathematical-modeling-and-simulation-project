function [V,M,N,t,tau_mat]=ml_4_rtc(csv_name, tau_mat_in, show_figures),

    params = readtable(csv_name);

    global Mcalcium_tot Mcalcium
    global Npotassium_tot Npotassium
    global vK vL vCa
    global gK gL gCa
    global C
    global tau1 tau2 tau3 tau4
    global T1 T2 T3 T4
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

    T1 = 0;
    T2 = 0;
    T3 = 0;
    T4 = 0;

    Mcalcium=M0; % initial state of calcium channel
    Npotassium=N0; % initial state of potassium channel

    % Functions for Morris-Lecar
    global Iapp; Iapp=@(t)100; % applied current
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

    % ODE options including reset
    options=odeset('Events',@nextevent);

    index = 0;
    tau_mat = tau_mat_in;

    %% Loop over events
    while t(end)<tmax

        index = index + 1;
        if size(tau_mat, 1) < index
            tau_arr = [-log(rand), -log(rand), -log(rand), -log(rand)];
            tau_mat = [tau_mat; tau_arr];
        end

        if index == 1
            tau1 = tau_mat(index, 1);
            tau2 = tau_mat(index, 2);
            tau3 = tau_mat(index, 3);
            tau4 = tau_mat(index, 4);
            index = index + 1;
            if size(tau_mat, 1) < index
                tau_arr = [-log(rand), -log(rand), -log(rand), -log(rand)];
                tau_mat = [tau_mat; tau_arr];
            end
        end

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

        %% Identify which reaction occurred, adjust state, and continue;
        % and set trigger for next event.
        mu=event_idx; % next reaction index
        %fprintf("Fired: %d\n", mu)
        if mu==1 % next reaction is a calcium channel opening
            M0=M0+1; % increment calcium channel state
            tau1=tau1+tau_mat(index, 1); % increment in tau1 is exponentially distributed with mean 1
        elseif mu==2 % next reaction is a calcium channel closing
            M0=M0-1; % decrement calcium channel state
            tau2=tau2+tau_mat(index, 2); % increment in tau2 is exponentially distributed with mean 1
        elseif mu==3 % next reaction is a potassium channel opening
            N0=N0+1;    % increment potassium channel state
            tau3=tau3+tau_mat(index, 3); % increment in tau3 is exponentially distributed with mean 1
        elseif mu==4 % next reaction is a potassium channel closing
            N0=N0-1;    % decrement potassium channel state
            tau4=tau4+tau_mat(index, 4); % increment in tau4 likewise
        end

        Mcalcium=M0;
        Npotassium=N0;
        if M0>Mcalcium_tot, error('M>Mtot'), end
        if M0<0, error('M<0'), end
        if N0>Npotassium_tot, error('N>Ntot'), end
        if N0<0, error('N<0'), end
        T1=T1+Uout(end,2);
        T2=T2+Uout(end,3);
        T3=T3+Uout(end,4);
        T4=T4+Uout(end,5);

        %fprintf("%f %f %f %f \n", T1, T2, T3, T4)

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
end % End of function mlexactboth

    %% Define the RHS for Morris-Lecar

    function dudt=dudtfunc(t,u)

        global Iapp % applied current
        global Mcalcium Mcalcium_tot
        global Npotassium Npotassium_tot
        global alpha_m beta_m
        global alpha_n beta_n
        global vK vL vCa
        global gK gL gCa
        global C


        v=u(1); % extract the voltage from the input vector

        dudt=[
            (Iapp(t)-gCa*(Mcalcium/Mcalcium_tot)*(v-vCa)-gL*(v-vL)-gK*(Npotassium/Npotassium_tot)*(v-vK))/C;
            alpha_m(v)*(Mcalcium_tot-Mcalcium); % Calcium chan. opening, internal time elapsed
            beta_m(v)*Mcalcium; % Calcium chan. closing, internal time elapsed
            alpha_n(v)*(Npotassium_tot-Npotassium); % Potassium chan. opening, internal time elapsed
            beta_n(v)*Npotassium; % Potassium chan. closing, internal time elapsed
            0; % M is constant between events
            0]; % N is constant between events
    end

function [value,isterminal,direction] = nextevent(~,u)

        global tau1 T1 % timing trigger for reaction 1 (Calcium opening)
        global tau2 T2 % timing trigger for reaction 2 (Calcium closing)
        global tau3 T3 % timing trigger for reaction 3 (Potassium opening)
        global tau4 T4 % timing trigger for reaction 4 (Potassium closing)

        %fprintf("%f %f %f %f\n", u(2)-(tau1-T1), u(3)-(tau2-T2), u(4)-(tau3-T3), u(5)-(tau4-T4))

        value=[u(2)-(tau1-T1);u(3)-(tau2-T2);u(4)-(tau3-T3);u(5)-(tau4-T4)];
        isterminal=[1;1;1;1]; % stop and restart integration at crossing
        direction=[1;1;1;1]; % increasing value of the quantity at the trigger
end
