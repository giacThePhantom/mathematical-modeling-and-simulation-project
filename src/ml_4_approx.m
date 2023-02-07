
function [V,M,N,t,Mtot,Ntot]=ml_rtc()

    % Calcium
    global Mcalcium_tot % total number of Calcium channels
    global Mcalcium % number of Calcium channels in conducting state
    global tau1  T1 % dummy variables for time to next Ca-opening event
    global tau2  T2 % dummy variables for time to next Ca-closing event

    Mcalcium_tot = 40;
    Mtot = Mcalcium;

    % Potassium
    global Npotassium_tot % total number of Potassium channels
    global Npotassium % number of Potassium channels in conducting state
    global tau3 T3 % dummy variables for time to next K-opening event
    global tau4 T4 % dummy variables for time to next K-closing event

    Npotassium_tot = 40;
    Ntot = Npotassium_tot;

    tmax=4e3;

    % Parameters
    phi_m = 0.4;
    va = -1.2;
    vb=18;
    vc = 2;
    vd = 30;
    phi_n = 0.04;

    % Initial conditions
    t0=0;
    t=t0; % global "external" time
    tau1=-log(rand); % time of next event on reaction stream 1 ("internal time")
    tau2=-log(rand); % time of next event on reaction stream 2 ("internal time")
    tau3=-log(rand); % time of next event on reaction stream 3 ("internal time")
    tau4=-log(rand); % time of next event on reaction stream 4 ("internal time")
    T1=0; % integrated intensity function for reaction 1
    T2=0; % integrated intensity function for reaction 2
    T3=0; % integrated intensity function for reaction 3
    T4=0; % integrated intensity function for reaction 4
    V0=-50; % start at an arbitrary middle voltage
    V=V0; % initial voltage
    M0=0; % start with calcium channels all closed
    Mcalcium=M0; % initial state of calcium channel
    M=M0; % use M to record the time course
    N0=ceil(Ntot/2); % start with half potassium channels open
    Npotassium=N0; % initial state of potassium channel
    N=N0; % use N to record the time course

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

    %% Loop over events
    while t(end)<tmax
        global t_start
        t_start = t(end);
        U0=[V0;M0;N0];
        tspan=[t(end),tmax];
        [tout,Uout,~,~,event_idx]=ode23(@dudtfunc,tspan,U0,options);
        Vout=Uout(:,1); % voltage at time of next event
        Mout=Uout(:,2); % number of calcium channels open at end of next event
        Nout=Uout(:,3); % number of potassium channels open at end of next event
        t=[t,tout'];
        V=[V,Vout'];
        M=[M,Mout'];
        N=[N,Nout'];

        %% Identify which reaction occurred, adjust state, and continue;
        % and set trigger for next event.
        mu=event_idx; % next reaction index
        if mu==1 % next reaction is a calcium channel opening
            M0=M0+1; % increment calcium channel state
            tau1=tau1-log(rand); % increment in tau1 is exponentially distributed with mean 1
        elseif mu==2 % next reaction is a calcium channel closing
            M0=M0-1; % decrement calcium channel state
            tau2=tau2-log(rand); % increment in tau2 is exponentially distributed with mean 1
        elseif mu==3 % next reaction is a potassium channel opening
            N0=N0+1;    % increment potassium channel state
            tau3=tau3-log(rand); % increment in tau3 is exponentially distributed with mean 1
        elseif mu==4 % next reaction is a potassium channel closing
            N0=N0-1;    % decrement potassium channel state
            tau4=tau4-log(rand); % increment in tau4 likewise
        end
        Mcalcium=M0;
        Npotassium=N0;
        if M0>Mcalcium_tot, error('M>Mtot'), end
        if M0<0, error('M<0'), end
        if N0>Npotassium_tot, error('N>Ntot'), end
        if N0<0, error('N<0'), end

        global alpha_m_t alpha_n_t
        global beta_m_t beta_n_t

        T1=T1+alpha_m_t;
        T2=T2+beta_m_t;
        T3=T3+alpha_n_t;
        T4=T4+beta_n_t;
        V0=V(end);
    end % while t(end)<tmax

    %% Plot output
    figure
        subplot(5,1,1),plot(t,M),ylabel('M'), xlabel('Time')
        subplot(5,1,2),plot(t,N),ylabel('N'), xlabel('Time')
        subplot(5,1,3),plot(t,V),xlabel('Time'), ylabel('V')
        subplot(5,1,4),plot(V,M,'.-'),xlabel('V'), ylabel('M')
        subplot(5,1,5),plot(V,N,'.-'),xlabel('V'), ylabel('N')
        grid on
    end % End of function mlexactboth

    %% Define the RHS for Morris-Lecar

    function dudt=dudtfunc(t,u)
        global Iapp % applied current
        global Mcalcium Mcalcium_tot
        global Npotassium Npotassium_tot
        global alpha_m beta_m
        global alpha_n beta_n
        global alpha_m_t; alpha_m_t = 0;
        global beta_m_t; beta_m_t = 0;
        global alpha_n_t; alpha_n_t = 0;
        global beta_n_t; beta_n_t = 0;
        global t_start
        delta_t = t - t_start;
        v=u(1); % extract the voltage from the input vector
        alpha_m_t = delta_t * alpha_m(v)*(Mcalcium_tot-Mcalcium);
        beta_m_t = delta_t * beta_m(v)*Mcalcium;
            % Potassium chan. opening, internal time elapsed
        alpha_n_t = delta_t * alpha_n(v)*(Npotassium_tot-Npotassium);
            % Potassium chan. closing, internal time elapsed
        beta_n_t = delta_t * beta_n(v)*Npotassium;


        vK = -84;
        vL = -60;
        vCa = 120;
        gK = 8;
        gL = 2;
        C=20;
        gCa = 4.4;

        dudt=[
            (Iapp(t)-gCa*(Mcalcium/Mcalcium_tot)*(v-vCa)-gL*(v-vL)...
            -gK*(Npotassium/Npotassium_tot)*(v-vK))/C;
            % Calcium chan. opening, internal time elapsed
            0; % M is constant between events
            0]; % N is constant between events
    end

function [value,isterminal,direction] = nextevent(~,u)
        global tau1 T1 % timing trigger for reaction 1 (Calcium opening)
        global tau2 T2 % timing trigger for reaction 2 (Calcium closing)
        global tau3 T3 % timing trigger for reaction 3 (Potassium opening)
        global tau4 T4 % timing trigger for reaction 4 (Potassium closing)
        global alpha_m_t beta_m_t
        global alpha_n_t beta_n_t
        alpha_m_t
        value=[alpha_m_t-(tau1-T1);beta_m_t-(tau2-T2);alpha_n_t-(tau3-T3);beta_n_t-(tau4-T4)];
        isterminal=[1;1;1;1]; % stop and restart integration at crossing
        direction=[1;1;1;1]; % increasing value of the quantity at the trigger
end
