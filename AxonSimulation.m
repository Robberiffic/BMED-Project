%%% Project Part 1
%%% Create and enter values for variables

% Constants
gK = 36e-3;
gNa = 120e-3;
gL = .3e-3;
eK = -12e-3;
eNa = 115e-3;
eL = 10.6e-3;
Vrest = -70e-3;

% Currents
I_Na = m.^3.*gNa.*h.*(Vm-eNa);
I_K = n.^4.*gK.*h.*(Vm-eK);
I_L = gL.*(Vm-eL);
I_ion = I-I_K-I_Na-I_L;


% Derivatives
dVm = I_ion./Cm;
dm = alpha_m.*(1-m)-beta_m.*m;
dn = alpha_n.*(1-n)-beta_n.*n;
dh = alpha_h.*(1-h)-beta_h.*h;
Euler's Method
y_n1 = y_n + h.*ff;

% Create time vector
t = linspace(t_0, delta_t, t_f);
% Create Vm vector to correspond to time vector
%

% Calculate Gating Variables from Vm Vector
alpha_m = .1.*((25-Vm)./(exp((25-Vm)/10)-1));
beta_m = 4.*exp(-Vm./18);
alpha_n = .1.*((10-Vm)./(exp((10-Vm)/10)-1));
beta_n = 125.*exp(-Vm./80);
alpha_h = .07.*exp(-Vm./20);
beta_h = 1./(exp((30-Vm)./10)+1);

% Set intial conditions for gating variables
% n0 =
% m0 =
% h0 =

% Parse through every element of the time/membrane voltage vectors,
% calculating m,n,h, conductances, and currents

for ind = 0:length(t)

    
    
end



