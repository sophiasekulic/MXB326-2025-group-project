%% Semi analytical solution

function S = semi_analytical_solution(Swr,Sor,beta,F,L,N,M,T)
% Inputs
%   Swr: min. saturation of water (constant)
%   Sor: min. saturation of oil (constant)
%   beta: rel. magnitude of viscous to capillary terms
%   F: water-to-oil viscosity ratio
% Outputs
%   S: solution matrix (Nx(M+1))

% initialise - explain chosen domains in report
x = linspace(0,L,N)'; % x: node locations (vector) with x(1) = 0 and x(end) = L 
dt = T/M; % time step duration
t = 0:dt:T; % time vector 
S = ones(N,M) * (1-Swr); % initial condition
gamma = beta*(1-Swr-F*Sor)/(F-1); % rearrange F equation 
v = gamma/beta; % viscosity ratio
alpha = beta^2 * (-((Sor+v)*(1-Swr+v))/(1-Swr-Sor)); % formula given
omega = beta - (alpha/(beta-beta*Swr+gamma)); % formula given

% create symbolic function
syms xt t_sym real 
phiExp = (1/2)*erfc(xt/(2*sqrt(t_sym)))...
     + (1/2)*exp(-omega*xt+omega^2*t_sym)*erfc(xt/(2*sqrt(t_sym))-omega*sqrt(t_sym))...
     - (1/2)*exp(-beta*xt+beta^2*t_sym)*erfc(xt/(2*sqrt(t_sym))-beta*sqrt(t_sym))...
     - (1/2)*exp(-(omega-beta)*xt+(omega-beta)^2*t_sym)*erfc(xt/(2*sqrt(t_sym))-(omega-beta)*sqrt(t_sym))...
     + exp((beta-omega)*xt+(beta-omega)^2*t_sym);
% differentiate wrt xt
dphiExp = diff(phiExp,xt);

% convert sym to string
dphiExp_str = char(dphiExp);
dphiExp_str = strrep(dphiExp_str,'xt','x');
dphiExp_str = strrep(dphiExp_str,'t_sym','t');

% create anonymous function
phiFunc = @(xt, t) ...
    (1/2)*erfc(xt/(2*sqrt(t))) ...
    + (1/2)*exp(-omega*xt + omega^2*t)*erfc(xt/(2*sqrt(t)) - omega*sqrt(t)) ...
    - (1/2)*exp(-beta*xt + beta^2*t)*erfc(xt/(2*sqrt(t)) - beta*sqrt(t)) ...
    - (1/2)*exp(-(omega-beta)*xt + (omega-beta)^2*t)*erfc(xt/(2*sqrt(t)) - (omega-beta)*sqrt(t)) ...
    + exp((beta-omega)*xt + (beta-omega)^2*t);

dphiFunc = str2func(['@(x,t)', dphiExp_str]);

for n = 2:M+1
    for i = 1:N
        % define rhs of nonlinear eqn (15)
        rhs = exp(alpha*x(i));
        t_n = t(n);
        % initial guess for xt using trapezoidal
        if i == 1 
            xt0 = (beta*S(1,n-1)+gamma)*x(1)+omega*t(n-1); % no integral since x(1) = 0
        else 
            xt0 = trapz(x(1:i), beta*S(1:i,n-1)+gamma)+omega*t(n-1);
            xt0 = max(min(xt0, 100), -100);
        end
        % solve phi(xt,t) = exp(alpha*x)
        xt = fzero(@(xt) phiFunc(xt,t_n) - rhs, xt0);
        % evaluate derivative
        dphi_n = dphiFunc(xt,t_n);
        % compute S using semi analytical formula
        S(i,n) = (1/beta)*((alpha*exp(alpha*x(i)))/dphi_n-gamma);
    end
end 
