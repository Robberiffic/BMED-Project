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

% Gating Variables
alpha_m = .1.*((25-Vm)./(exp((25-Vm)/10)-1));
beta_m = 4.*exp(-Vm./18);
alpha_n = .1.*((10-Vm)./(exp((10-Vm)/10)-1));
beta_n = 125.*exp(-Vm./80);
alpha_h = .07.*exp(-Vm./20);
beta_h = 1./(exp((30-Vm)./10)+1);

% Derivatives
dVm = I_ion./Cm;
dm = alpha_m.*(1-m)-beta_m.*m;
dn = alpha_n.*(1-n)-beta_n.*n;
dh = alpha_h.*(1-h)-beta_h.*h;
% Euler's Method
y_n1 = y_n + h.*ff;



