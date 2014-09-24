%%% Project Part 1

% Constants
gK = 36e-3;
gNa = 120e-3;
gL = .3e-3;
eK = -12e-3;
eNa = 115e-3;
eL = 10.6e-3;
Vrest = -70e-3;

% Derivatives
dVm = I_ion./Cm;
dm = alpha_m.*(1-m)-beta_m.*m;
dn = alpha_n.*(1-n)-beta_n.*n;
dh = alpha_h.*(1-h)-beta_h.*h;
Euler's Method
y_n1 = y_n + h.*ff;

% Timing parameters
t_0 = 0;
%t_f = 
delta_t = (t_f-t_0)./8000;

% Create time vector
t = linspace(t_0, delta_t, t_f);

% Create Vm vector to correspond to time vector
%

% Initialize vectors for n,m,h to be same size as time vector
n = t;
m = t;
h = t;

% Initialize alphas and betas to be same size as time vector
alpha_m = t;
beta_m = t;
alpha_n = t;
beta_n = t;
alpha_h = t;
beta_h = t;

% Set intial conditions for n,m,h
% n(1) = 
% m(1) =
% h(1) =

% Parse through every element of the time/membrane voltage vectors,
% calculating m,n,h

for ind = 1:length(t)-1

    alpha_m(ind) = .1.*((25-Vm(ind))./(exp((25-Vm(ind))/10)-1));
    beta_m(ind) = 4.*exp(-Vm(ind)./18);
    alpha_n(ind) = .1.*((10-Vm(ind))./(exp((10-Vm(ind))/10)-1));
    beta_n(ind) = 125.*exp(-Vm(ind)./80);
    alpha_h(ind) = .07.*exp(-Vm(ind)./20);
    beta_h(ind) = 1./(exp((30-Vm(ind))./10)+1);
    
    % Rearranged first order differential equation defining m,n,h; solving via
    % Euler's Method
    m(ind+1) = m(ind) + delta_t.*(alpha_m(ind) - m(ind).*(alpha_m(ind)-beta_m(ind)));
    n(ind+1) = n(ind) + delta_t.*(alpha_n(ind) - n(ind).*(alpha_n(ind)-beta_n(ind)));
    h(ind+1) = h(ind) + delta_t.*(alpha_h(ind) - h(ind).*(alpha_h(ind)-beta_h(ind)));
    
end

% Currents
I_Na = m.^3.*gNa.*h.*(Vm-eNa);
I_K = n.^4.*gK.*h.*(Vm-eK);
I_L = gL.*(Vm-eL);
I_ion = I-I_K-I_Na-I_L;

