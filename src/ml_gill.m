function [V,N,t,Ntot]=ml_gill()
    
    global Npotassium_tot
    global Npotassium
    global tau1 T1
    global tau2 T2

    Npotassium_tot = 40;
    Ntot = Npotassium_tot;
    tmax=4e3;
    
    % Parameters
    va = -1.2;
    vb = 18;
    vc = 2;
    vd = 30;
    phi = 0.04;

    % Initial conditions
    t0=0;
    t=t0; % global "external" time
    tau1=-log(rand); % time of next event on reaction stream 1 ("internal time")
    T1=0; % integrated intensity function for reaction 1
    T2=0; % integrated intensity function for reaction 2
    V0=-50; % start at an arbitrary middle voltage
    V=V0;
    N0=ceil(Ntot/2); % start with half of channels open
    Npotassium=N0;
    N=N0; % use N to record the time course

    % Morris Lecar
    global Iapp; Iapp=@(t)100; % applied current
    global minf; minf=@(v)0.5*(1+tanh((v-va)/vb)); % m-gate activation
    global xi; xi=@(v)(v-vc)/vd;  % scaled argument for n-gate input
    global ninf; ninf=@(v)0.5*(1+tanh(xi(v)));  % n-gate activation function
    global tau_n; tau_n=@(v)1./(phi*cosh(xi(v)/2)); % n-gate activation t-const
    global alpha; alpha=@(v)(ninf(v)./tau_n(v)); % per capita opening rate
    global beta; beta=@(v)((1-ninf(v))./tau_n(v)); % per capita closing rate

    options=odeset('Events',@nextevent);

    while t(end)<tmax
        
        U0=[V0;0;0;N0];
        tspan=[t(end),tmax];
        [tout,Uout,~,~,event_idx]=ode23(@dudtfunc,tspan,U0,options);
        Vout=Uout(:,1); % voltage at time of next event
        Nout=Uout(:,4); % number of channels open at end of next event
        t=[t,tout'];
        V=[V,Vout'];
        N=[N,Nout'];

        espilon = rand;

        q_k = (alpha(V(end)) * (Ntot-N0))/((alpha(V(end)) * (Ntot-N0)) + (beta(V(end)) * N0));

        if (espilon < q_k)
            N0=N0+1; 
        else
            N0=N0-1;   % increment channel state
        end
    
        Npotassium=N0;
        if N0>Npotassium_tot, error('N>Ntot'), end
        if N0<0, error('N<0'), end
        T1=T1+Uout(end,2);
        T2=T2+Uout(end,3);
        V0=V(end);
    end % while t(end)<tmax
    
    %% Plot output
    figure
    subplot(3,1,1),plot(t,V),xlabel('time'),ylabel('V')
    subplot(3,1,2),plot(t,N),xlabel('time'),ylabel('N')
    subplot(3,1,3),plot(V, N, '-.'),xlabel('V'),ylabel('N')
    %subplot(3,1,3),plot(t,Iapp),xlabel('time'),ylabel('I')
    shg
end

function dudt=dudtfunc(t,u)
       
    global Iapp  % Applied Current
    global minf  % asymptotic target for (deterministic) Calcium channel
    global Npotassium Npotassium_tot  % num. open, total num. of channels
    global alpha beta  % per capita transition rates
    
    vK = -84;
    vL = -60;
    vCa = 120;
    gK = 8;
    gL = 2;
    C = 20;
    gCa = 4.4;

    v=u(1); % extract the voltage from the input vect
    dudt=[(Iapp(t)-gCa*minf(v)*(v-vCa)-gL*(v-vL)-gK*(Npotassium/Npotassium_tot)*(v-vK))/C; % voltage
    alpha(v)*(Npotassium_tot-Npotassium) + beta(v)*Npotassium; % channel opening internal time
    beta(v)*Npotassium;% channel closing internal time
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