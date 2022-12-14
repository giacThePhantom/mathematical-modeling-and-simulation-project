function [V,N,t,Ntot]=ml_gill_exact_konly(tmax,Ntot)

%function [V,N,t,Ntot]=ml_rtc_exact_konly(tmax,Ntot);
%
% Exact solution of Morris Lecar with stochastic potassium channel. 
% Using the random time change algorithm.  We track two reactions: 
%
% Rxn 1: closed -> open (per capita rate alpha)
% Rxn 2: open -> closed (per capita rate beta)
%
% Default Ntot=40, tmax=4000.  
%
% Author: PJT July 2013, Case Western Reserve University.

% Applied current "Iapp" set on line 47 below.

%% Use global variables to represent the channel state
%    and random trigger for Poisson process
global Npotassium_tot
global Npotassium
global tau1  T1
global tau2  T2
% total number of Potassium channels
% number of open Potassium channels
% time to next opening event (internal to reaction 1)
% time to next closing event (internal to reaction 2)

%% Set defaults for input arguments
if nargin < 2, Ntot=40; end
Npotassium_tot=Ntot;
if nargin < 1, tmax=4e3; end

%% Parameters
% Standard Morris-Lecar parameters giving a globally attracting limit cycle
% (if applied current Iapp=100) or a stable fixed point (if Iapp=75).
va = -1.2;
vb=18;
vc = 2;
vd = 30;
phi = 0.04;

%% Functions for Morris-Lecar
global Iapp; Iapp=@(t)100; % applied current.
global minf; minf=@(v)0.5*(1+tanh((v-va)/vb)); % m-gate activation
global xi; xi=@(v)(v-vc)/vd;  % scaled argument for n-gate input
global ninf; ninf=@(v)0.5*(1+tanh(xi(v)));  % n-gate activation function
global tau_n; tau_n=@(v)1./(phi*cosh(xi(v)/2)); % n-gate activation t-const
global alpha; alpha=@(v)(ninf(v)./tau_n(v)); % per capita opening rate
global beta; beta=@(v)((1-ninf(v))./tau_n(v)); % per capita closing rate

%% ODE options including reset
options=odeset('Events',@nextevent);

%% Initialize
t0=0;
t=t0; % global "external" time
tau1=-log(rand); % time of next event on reaction stream 1 ("internal time")
tau2=-log(rand); % time of next event on reaction stream 2 ("internal time")
T1=0; % integrated intensity function for reaction 1
T2=0; % integrated intensity function for reaction 2
V0=-50; % start at an arbitrary middle voltage
V=V0;
N0=ceil(Ntot/2); % start with half of channels open
Npotassium=N0;
N=N0; % use N to record the time course

%% Loop over events
while t(end)<tmax

%% Integrate ODE for voltage, until the next event is triggered
% State vector U for integration contains the following components
%[voltage;
% integral of (# closed)*alpha(v(t));
% integral of (# open)*beta(v(t));
% # open].
U0=[V0;0;0;N0];
tspan=[t(end),tmax];
[tout,Uout,~,~,event_idx]=ode23(@dudtfunc,tspan,U0,options);
Vout=Uout(:,1); % voltage at time of next event
Nout=Uout(:,4); % number of channels open at end of next event
t=[t,tout'];
V=[V,Vout'];
N=[N,Nout'];

%% Identify which reaction occurred, adjust state, and
%  continue, and set trigger for next event.
mu=event_idx; % next reaction index
if mu==1 % next reaction is a channel opening
    N0=N0+1;    % increment channel state
    tau1=tau1-log(rand); % increment in tau1 is exp'l with mean 1
elseif mu==2 % next reaction is another channel closing
    N0=N0-1;    % decrement channel state
    tau2=tau2-log(rand); % increment in tau2 likewise
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
subplot(4,1,1),plot(t,N),xlabel('time'),ylabel('N')
subplot(4,1,2:3),plot(V,N,'.-'),xlabel('V'),ylabel('N')
subplot(4,1,4),plot(t,V),xlabel('time'),ylabel('V')
shg
end

%% Define the RHS for Morris-Lecar
% state vector for integration is as follows:
% u(1) = voltage
% u(2) = integral of activation hazard function
% u(3) = integral of inactivation hazard function
% u(4) = N (number of open potassium channels)
function dudt=dudtfunc(t,u)
global Iapp  % Applied Current
global minf  % asymptotic target for (deterministic) Calcium channel
global Npotassium Npotassium_tot  % num. open, total num. of channels
global alpha beta  % per capita transition rates

%% Parameters
vK = -84;
vL = -60;
vCa = 120;
gK =8;
gL =2;
C=20;
gCa = 4.4;

%% calculate the RHS;
v=u(1); % extract the voltage from the input vector
dudt=[(Iapp(t)-gCa*minf(v)*(v-vCa)-gL*(v-vL)-...
        gK*(Npotassium/Npotassium_tot)*(v-vK))/C; % voltage
    alpha(v)*(Npotassium_tot-Npotassium); % channel opening internal time
    beta(v)*Npotassium;% channel closing internal time
0]; % N is constant between events
end

%% Define behavior at threshold crossing
function [value,isterminal,direction] = nextevent(~,u)
    global tau1 T1 % timing trigger for reaction 1 (opening)
    global tau2 T2 % timing trigger for reaction 2 (closing)
    value=[u(2)-(tau1-T1);u(3)-(tau2-T2)];
    isterminal=[1;1]; % stop and restart integration at crossing
    direction=[1;1]; % increasing value of the quantity at the trigger
end