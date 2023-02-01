function [V,N,t,Ntot]=ml_gill()
    close all;
    clear all;

    global Mcalcium_tot % total number of Calcium channels
    global Mcalcium % number of Calcium channels in conducting state
    global Npotassium_tot
    global Npotassium
    global tau1  T1
    global tau2  T2

    Npotassium_tot = 40;
    Ntot = Npotassium_tot;
    Mtot = 40;
    Mcalcium_tot = Mtot;
    tmax=4e3;

    % Parameters
    phi_m=0.4;
    va = -1.2; vb=18;
    vc = 2; vd = 30;
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

    % Morris Lecar
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

    options=odeset('Events',@nextevent);

    while t(end)<tmax

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

        epsilon = rand;

        lambda_0 = (alpha_n(V(end)) * (Ntot-N0)) + (beta_n(V(end)) * N0) + (alpha_m(V(end)) * (Mtot-M0)) + (beta_m(V(end)) * M0);
        lambda_1 = alpha_n(V(end)) * (Ntot-N0);
        lambda_2 = lambda_1 + (beta_n(V(end)) * (N0));
        lambda_3 = lambda_2 + (alpha_m(V(end)) * (Mtot-M0));

        q_1 = lambda_1/lambda_0;
        q_2 = lambda_2/lambda_0;
        q_3 = lambda_3/lambda_0;

        if (epsilon < q_1)
            N0=N0+1;
        elseif (epsilon >= q_1 & epsilon < q_2)
            N0=N0-1;   % increment channel state
        elseif (epsilon >= q_2 & epsilon < q_3)
            M0 = M0 + 1;
        else
            M0 = M0 - 1;
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
        V0=V(end);

    end % while t(end)<tmax

    %% Plot output
    figure
    subplot(6,1,1),plot(t,M),ylabel('M','FontSize',20),set(gca,'FontSize',20)
    subplot(6,1,2),plot(t,N),ylabel('N','FontSize',20),set(gca,'FontSize',20)
    subplot(6,1,6),plot(t,V),xlabel('Time','FontSize',20)
        ylabel('V','FontSize',20),set(gca,'FontSize',20)
    subplot(6,1,4:6),plot3(V,M,N,'.-'),xlabel('V','FontSize',20)
        ylabel('M','FontSize',20),zlabel('N','FontSize',20),set(gca,'FontSize',20)
    grid on, rotate3d, shg

end

function dudt=dudtfunc(t,u)

    global Iapp % applied current
    global Mcalcium Mcalcium_tot
    global Npotassium Npotassium_tot
    global alpha_m beta_m
    global alpha_n beta_n

    %% Parameters
    vK = -84; vL = -60; vCa = 120;
    gK =8; gL =2; C=20; gCa = 4.4;

    v=u(1); % extract the voltage from the input vect
dudt=[...  % voltage
    (Iapp(t)-gCa*(Mcalcium/Mcalcium_tot)*(v-vCa)-gL*(v-vL)...
    -gK*(Npotassium/Npotassium_tot)*(v-vK))/C;
    % Calcium chan. opening, internal time elapsed
    alpha_m(v)*(Mcalcium_tot-Mcalcium);
    % Calcium chan. closing, internal time elapsed
    beta_m(v)*Mcalcium;
    % Potassium chan. opening, internal time elapsed
    alpha_n(v)*(Npotassium_tot-Npotassium);
    % Potassium chan. closing, internal time elapsed
    beta_n(v)*Npotassium;
    0; % M is constant between events
    0]; % N is constant between events
end

% Define behavior at threshold crossing
function [value,isterminal,direction] = nextevent(~,u)
    global tau1 T1 % timing trigger for reaction 1 (opening)
    global tau2 T2 % timing trigger for reaction 2 (closing)
    value=[u(2)-tau1];
    isterminal=[1]; % stop and restart integration at crossing
    direction=[1]; % increasing value of the quantity at the trigger
end
