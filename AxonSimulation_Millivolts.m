%%% Project Part 1

close all
clear all

% Constants
gK = (36);
gNa = (120);
gL = (.3);
eK = -12;
eNa = 115;
eL = 10.6;
Cm = 1;
Vrest = 0; % This variable is the difference between membrane voltage and resting voltage

% Timing parameters (in milliseconds)
t_0 = 0;
t_f = 100;
t_num = 200000;
delta_t = (t_f-t_0)./t_num; % Change for more/fewer time steps

% Create time vector
t = linspace(t_0, t_f, t_num);

% Create Vm vector to correspond to time vector
Vm = zeros(1,length(t));

% Initialize vectors for n,m,h to be same size as time vector
n = zeros(1,length(t));
m = zeros(1,length(t));
h = zeros(1,length(t));

% Define excitation current for Problem 1
% I_in = zeros(1,length(t));

% Define excitation current for Problem 2
% I_in = zeros(1,length(t));
% I_in(1:round(.5/100.*length(t))) = (5);

% Define excitation current for Problem 3
I_in = zeros(1,length(t));
I_in(1:end) = (5);

% Initialize current vectors
I_Na = zeros(1,length(t));
I_K = zeros(1,length(t));
I_L = zeros(1,length(t));
I_ion = I_in-I_K-I_Na-I_L;

% Initialize alphas and betas to be same size as time vector
alpha_m = zeros(1,length(t));
beta_m = zeros(1,length(t));
alpha_n = zeros(1,length(t));
beta_n = zeros(1,length(t));
alpha_h = zeros(1,length(t));
beta_h = zeros(1,length(t));


% Set intial conditions for n,m,h
alpha_m0 = .1.*((25-Vrest)./(exp((25-Vrest)./10)-1));
beta_m0 = 4.*exp(-Vrest./18);
alpha_n0 = .01.*(10-Vrest)./(exp((10-Vrest)./10)-1);
beta_n0 = .125.*exp(-Vrest./80);
alpha_h0 = .07.*exp(-Vrest./20);
beta_h0 = 1./(exp((30-Vrest)./10)+1);

alpha_m(1) = alpha_m0;
beta_m(1) = beta_m0;
alpha_n(1) = alpha_n0;
beta_n(1) = beta_n0;
alpha_h(1) = alpha_h0;
beta_h(1) = beta_h0;

n(1) = alpha_n0/(alpha_n0 + beta_n0);
m(1) = alpha_m0/(alpha_m0 + beta_m0);
h(1) = alpha_h0/(alpha_h0 + beta_h0);

% Set initial conditions for Vm
Vm(1) = Vrest;

% Set intial conditons for current
I_ion(1) = 0;

% Parse through every element of the time/membrane voltage vectors,
% calculating m,n,h
for ind = 1:length(t)-1

    Vm(ind+1) = Vm(ind) + delta_t.*(I_ion(ind)./Cm);
    
    alpha_m(ind+1) = .1.*((25-Vm(ind+1))./(exp((25-Vm(ind+1))./10)-1));
    beta_m(ind+1) = 4.*exp(-Vm(ind+1)./18);
    alpha_n(ind+1) = .01.*((10-Vm(ind+1))./(exp((10-Vm(ind+1))./10)-1));
    beta_n(ind+1) = .125.*exp(-Vm(ind+1)./80);
    alpha_h(ind+1) = .07.*exp(-Vm(ind+1)./20);
    beta_h(ind+1) = 1./(exp((30-Vm(ind+1))./10)+1);
    
    % Rearranged first order differential equation defining m,n,h; solving via
    % Euler's Method
    m(ind+1) = m(ind) + delta_t.*(alpha_m(ind).*(1 - m(ind))-m(ind).*beta_m(ind));
    n(ind+1) = n(ind) + delta_t.*(alpha_n(ind).*(1 - n(ind))-n(ind).*beta_n(ind));
    h(ind+1) = h(ind)   + delta_t.*(alpha_h(ind).*(1 - h(ind))-h(ind).*beta_h(ind));
        
    I_Na(ind+1) = (m(ind+1)).^3.*gNa.*h(ind+1).*(Vm(ind+1)-eNa);
    I_K(ind+1) = (n(ind+1)).^4.*gK.*(Vm(ind+1)-eK);
    I_L(ind+1) = gL.*(Vm(ind+1)-eL);
    I_ion(ind+1) = (I_in(ind+1)-I_K(ind+1)-I_Na(ind+1)-I_L(ind+1));
    
        
end
hold on

%%% Plot Membrane Potential

% plot(t,Vm-70, '-r');
% title('Membrane Potential');
% xlabel('Time (milliSeconds)');
% ylabel('Voltage (milliVolts)');

%%% Plot Channel Conductance

plot(t,(m).^3.*gNa.*h, '-r');
plot(t,(n).^4.*gK, '-g');
title('Conductances of Sodium and Potassium Channels');
xlabel('Time (milliSeconds)');
ylabel('Conductance (milliSiemens/cm^2)');
legend('Na+', 'K+');

