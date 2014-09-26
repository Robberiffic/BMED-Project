%%% Project Part 1

close all
clear all

% Constants
gK = (36e-3).*100.^2;
gNa = (120e-3).*100.^2;
gL = (.3e-3).*100.^2;
eK = -12e-3;
eNa = 115e-3;
eL = 10.6e-3;
Cm = 1e-6.*(100.^2);
Vrest = -70e-3;

% Timing parameters (in milliseconds)
t_0 = 0;
t_f = 100e-3;
t_num = 10000;
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
I_in = zeros(1,length(t));
% Define excitation current for Problem 2
% I_in = zeros(1,length(t));
% I_in(1:round(.5e-3./100e-3).*length(t)) = (5e-6).*(100.^2)
% Define excitation current for Problem 3
% I_in = zeros(1,length(t));

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
alpha_m0 = .1.*((25-1000.*Vrest)./(exp((25-1000.*Vrest)./10)-1));
beta_m0 = 4.*exp(-1000.*Vrest./18);
alpha_n0 = .01.*((10-1000.*Vrest)./(exp((10-1000.*Vrest)./10)-1));
beta_n0 = .125.*exp(-1000.*Vrest./80);
alpha_h0 = .07.*exp(-1000.*Vrest./20);
beta_h0 = 1./(exp((30-1000.*Vrest)./10)+1);

alpha_m(1) = .1.*((25-1000.*Vrest)./(exp((25-1000.*Vrest)./10)-1));
beta_m(1) = 4.*exp(-1000.*Vrest./18);
alpha_n(1) = .01.*((10-1000.*Vrest)./(exp((10-1000.*Vrest)./10)-1));
beta_n(1) = .125.*exp(-1000.*Vrest./80);
alpha_h(1) = .07.*exp(-1000.*Vrest./20);
beta_h(1) = 1./(exp((30-1000.*Vrest)./10)+1);

n(1) = alpha_n0./(alpha_n0+beta_n0);
m(1) = alpha_m0./(alpha_m0+beta_m0);
h(1) = alpha_h0./(alpha_h0+beta_h0);

% Set initial conditions for Vm
Vm(1) = Vrest;

% Set intial conditons for current
% I_Na(1) = m(1).^3.*gNa.*h(1).*(Vm(1)-eNa);
% I_K(1) = n(1).^4.*gK.*h(1).*(Vm(1)-eK);
% I_L(1) = gL.*(Vm(1)-eL);
% I_ion(1) = I_in(1)-I_K(1)-I_Na(1)-I_L(1);

% Parse through every element of the time/membrane voltage vectors,
% calculating m,n,h
for ind = 1:length(t)-1

      
    % Rearranged first order differential equation defining m,n,h; solving via
    % Euler's Method
    m(ind+1) = m(ind);% + delta_t.*(alpha_m(ind).*(1 - m(ind))-m(ind).*beta_m(ind));
    n(ind+1) = n(ind);% + delta_t.*(alpha_n(ind).*(1 - n(ind))-n(ind).*beta_n(ind));
    h(ind+1) = h(ind);% + delta_t.*(alpha_h(ind).*(1 - h(ind))-h(ind).*beta_h(ind));
    
    % Currents
    I_Na(ind) = (m(ind)).^3.*gNa.*h(ind).*(Vm(ind)-eNa);
    I_K(ind) = (n(ind)).^4.*gK.*(Vm(ind)-eK);
    I_L(ind) = gL.*(Vm(ind)-eL);
    I_ion(ind) = (I_in(ind)-I_K(ind)-I_Na(ind)-I_L(ind));
    
    Vm(ind+1) = Vm(ind) + delta_t.*(I_ion(ind)./Cm);
    
    alpha_m(ind+1) = .1.*((25-1000.*Vm(ind+1))./(exp((25-1000.*Vm(ind+1))./10)-1));
    beta_m(ind+1) = 4.*exp(-1000.*Vm(ind+1)./18);
    alpha_n(ind+1) = .01.*((10-1000.*Vm(ind+1))./(exp((10-1000.*Vm(ind+1))./10)-1));
    beta_n(ind+1) = .125.*exp(-1000.*Vm(ind+1)./80);
    alpha_h(ind+1) = .07.*exp(-1000.*Vm(ind+1)./20);
    beta_h(ind+1) = 1./(exp((30-1000.*Vm(ind+1))./10)+1);
    
end

plot(t,Vm);
