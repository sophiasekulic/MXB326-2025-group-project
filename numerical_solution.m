%% Numerical Solution

function Sn = numerical_solution(Swr,Sor,F_visc,beta,theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)
% Inputs
%   Swr: minimum saturation of water (constant)
%   Sor: minimum saturation of oil (constant)
%   F_visc: water-to-oil viscosity ratio (constant)
%   beta: relative magnitude of viscous to capillary terms (constant)
%   theta: theta method weighting parameter (theta = 0, 1/2 or 1)
%   sigma: weighted average approximation parameter (sigma = 0, 1/2 or 1)
%   N: number of nodes (positive integer)
%   M: number of time steps (positive integer)
%   T: end time (positive number)
%   x: node locations (vector)
%   BC: boundary condition at x = L, 'zeroflux' or 'zerogradient' (string)
%   jacobian: Jacobian generation options 'analytical','finitecolumn'
%       or 'finitetridiag' (string)
%   ls: line searching options 'simple', 'two-point' or 'three-point' (string)
%   linear_method: methods for solving the linear system at each Newton
%       iteration, 'full', 'sparse' or 'sparsediag'
%   atol: absolute error tolerence for Newton's method (positive scalar)
%   rtol: relative error tolerence for Newton's method (positive scalar)
%   maxiters: maximum number of Newton iterations (positive integer)
% Outputs
%   Sn: solution matrix of size N x M, where Sn(i,j) is the solution
%       S(x,t) evaluated at x = x(i) and t = t(j)

h = diff(x); % node spacings
dt = T/M; % size of time step
V = [h(1)/2; (h(1:N-2)+h(2:N-1))/2; h(N-1)/2]; % control volume lengths

% Physical characteristic parameters
gamma = beta*(1-Swr-F_visc*Sor)/(F_visc-1); % rearranged equation F 
alpha = -(beta^2)*((Sor+gamma/beta)*(1-Swr+gamma/beta))/(1-Swr-Sor); % formula provided

% Capillary-hydraulic properties of the fluid-porous system
f = @(S) alpha/(beta^2)*(1/(1-Swr+gamma/beta) - 1/(S+gamma/beta)); % f(S) provided
g = @(S) 1/((beta*S+gamma)^2); % g(S) provided
fprime = @(S) (alpha/beta^2)*(1/(S+gamma/beta)^2); % f'(S): derivative of f(S)
gprime = @(S) -2*beta/(beta*S+gamma)^3; % g'(S): derivative of g(S)

Sn = zeros(N,M+1); % initialising solution matrix
S0 = ones(N,1)*(1 - Swr); % initial condition
Sn(:,1) = S0; % setting the initial condition in the solution matrix

J = @(S) Jfunc(S,N,h,V,dt,theta,sigma,g,gprime,fprime,BC); %

% Time stepping
for n = 1:M
    F_newton = @(S) Ffunc(S,Sn(:,n),h,V,N,g,f,sigma,BC,theta,dt); % compute F_newton using solution at previous time step
    Sn(:,n+1) = newton_solver(F_newton,Sn(:,n),N,atol,rtol,maxiters,jacobian,linear_method,ls,J); % solution at current time step
end

% Subfunctions

% Gfunc: dS/dt = G(S)
function G = Gfunc(S,h,V,N,g,f,sigma,BC)
% Inputs
%   S: (S1,...,SN) (vector)
%   h: node spacings (vector)
%   V: control volume lengths (vector)
%   N: number of nodes (positive integer)
%   g: capillary-hydraulic properties g(S) (function handle)
%   f: capillary-hydraulic properties f(S) (function handle)
%   sigma: weighted average approximation parameter (sigma = 0, 1/2 or 1)
%   BC: boundary condition at x = L, 'zeroflux' or 'zerogradient' (string)
% Outputs
%   G: vector-valued function G evaluated at S (vector)

G = zeros(N,1);
G(1) = 1/V(1)* (-1+((1-sigma)*g(S(1))+sigma*g(S(2)))*(S(2)-S(1))/h(1)+(1-sigma)*f(S(1))+sigma*f(S(2)));
  for m = 2:N-1
    G(m) = 1/V(m) * (((sigma-1)*g(S(m-1))-sigma*g(S(m)))*(S(m)-S(m-1))/h(m-1)+(sigma-1)*f(S(m-1))-...
        sigma*f(S(m))+((1-sigma)*g(S(m))+sigma*g(S(m+1)))*(S(m+1)-S(m))/h(m)+(1-sigma)*f(S(m))+...
        sigma*f(S(m+1)));
  end
if isequal(BC, 'zeroflux')
G(N) = 1/V(N) * (((sigma-1)*g(S(N-1))-sigma*g(S(N)))*(S(N)-S(N-1))/h(N-1)+(sigma-1)*f(S(N-1))-...
        sigma*f(S(N)));
elseif isequal(BC, 'zerogradient')
G(N) = 1/V(N) * (((sigma-1)*g(S(N-1))-sigma*g(S(N)))*(S(N)-S(N-1))/h(N-1)+(sigma-1)*f(S(N-1))-...
        sigma*f(S(N))+f(S(N)));
end
end

%Ffunc: 0 = F(S)
function F = Ffunc(S,S_old,h,V,N,g,f,sigma,BC,theta,dt)
% Inputs
%   S: (S1,...,SN) (vector)
%   S_old: (S1,...,SN) at t = tn (vector)
%   h: node spacings (vector)
%   V: control volume lengths (vector)
%   N: number of nodes (positive integer)
%   g: capillary-hydraulic properties g(S) (function handle)
%   f: capillary-hydraulic properties f(S) (function handle)
%   sigma: weighted average approximation parameter (sigma = 0, 1/2 or 1)
%   BC: boundary condition at x = L, 'zeroflux' or 'zerogradient' (string)
%   theta: theta method weighting parameter (theta = 0, 1/2 or 1)
%   dt: time stepsize (positive scalar)
% Outputs
%   F: vector-valued function F evaluated at S (vector)

F = S - S_old - dt*theta*Gfunc(S,h,V,N,g,f,sigma,BC) -...
    dt*(1-theta)*Gfunc(S_old,h,V,N,g,f,sigma,BC);
end

%Jfunc: Jacobian of F(S)
function J = Jfunc(S,N,h,V,dt,theta,sigma,g,gprime,fprime,BC)
% Inputs
%   S: (S1,...,SN) (vector)
%   N: number of nodes (positive integer)
%   h: node spacings (vector)
%   V: control volume lengths (vector)
%   dt: time stepsize (positive scalar)
%   theta: theta method weighting parameter (theta = 0, 1/2 or 1)
%   sigma: weighted average approximation parameter (sigma = 0, 1/2 or 1)
%   g: capillary-hydraulic properties g(S) (function handle)
%   gprime: g'(S), derivative of g(S) (function handle)
%   fprime: f'(S), derivative of f(S) (function handle)
%   BC: boundary condition at x = L, 'zeroflux' or 'zerogradient' (string)
% Outputs
%   J: Jacobian matrix of F evaluated at S (matrix of size N by N)

J = zeros(N,N);
J(1,1) = 1-dt*theta/V(1)*((1-sigma)*gprime(S(1))*(S(2)-S(1))/h(1)-((1-sigma)*g(S(1))+sigma*g(S(2)))/h(1)+...
    (1-sigma)*fprime(S(1)));
J(1,2) = -dt*theta/V(1)*(sigma*gprime(S(2))*(S(2)-S(1))/h(1)+((1-sigma)*g(S(1))+sigma*g(S(2)))/h(1)+...
    sigma*fprime(S(2)));
   for i = 2:N-1
    J(i,i-1) = dt*theta/V(i)*((1-sigma)*gprime(S(i-1))*(S(i)-S(i-1))/h(i-1)-...
        ((1-sigma)*g(S(i-1))+sigma*g(S(i)))/h(i-1)+(1-sigma)*fprime(S(i-1)));
    J(i,i) = 1-dt*theta/V(i)*(-sigma*gprime(S(i))*(S(i)-S(i-1))/h(i-1)-...
        ((1-sigma)*g(S(i-1))+sigma*g(S(i)))/h(i-1)-sigma*fprime(S(i))+...
        (1-sigma)*gprime(S(i))*(S(i+1)-S(i))/h(i)-...
        ((1-sigma)*g(S(i))+sigma*g(S(i+1)))/h(i)+(1-sigma)*fprime(S(i)));
    J(i,i+1) = -dt*theta/V(i)*(sigma*gprime(S(i+1))*(S(i+1)-S(i))/h(i)+...
        ((1-sigma)*g(S(i))+sigma*g(S(i+1)))/h(i)+sigma*fprime(S(i+1)));
   end
J(N,N-1) = dt*theta/V(N)*((1-sigma)*gprime(S(N-1))*(S(N)-S(N-1))/h(N-1)-...
        ((1-sigma)*g(S(N-1))+sigma*g(S(N)))/h(N-1)+(1-sigma)*fprime(S(N-1)));
if isequal(BC, 'zeroflux')
J(N,N) = 1 + dt*theta/V(N)*(sigma*gprime(S(N))*(S(N)-S(N-1))/h(N-1)+...
        ((1-sigma)*g(S(N-1))+sigma*g(S(N)))/h(N-1)+sigma*fprime(S(N)));
elseif isequal(BC, 'zerogradient')
J(N,N) = 1 + dt*theta/V(N)*(sigma*gprime(S(N))*(S(N)-S(N-1))/h(N-1)+...
        ((1-sigma)*g(S(N-1))+sigma*g(S(N)))/h(N-1)+sigma*fprime(S(N))+fprime(S(N)));
end
end

% Newton's Method: Solving F(S) = 0
function S = newton_solver(F_newton,S0,N,atol,rtol,maxiters,jacobian,linear_method,ls,J)
% Inputs
%   F_newton: vector-valued function, F_newton(S) (function handle)
%   S0: initial guess (S1,...,SN) (vector)
%   N: number of nodes (positive integer)
%   atol: absolute error tolerence for Newton iteration (positive scalar)
%   rtol: relative error tolerence for Newton iteration (positive scalar)
%   maxiters: maximum number of Newton iterations (positive integer)
%   jacobian: Jacobian generation options 'analytical','finitecolumn'
%       or 'finitetridiag' (string)
%   linear_method: methods for solving the linear system at each Newton
%       iteration, 'sparse' or 'sparsediag'
%   ls: line searching options 'simple', 'two-point' or 'three-point' (string)
%   J: matrix-valued function of analytical Jacobian, J(S) (function handle)
% Outputs
%   S: approximate solution of F(S) = 0 (vector)

I = eye(N); sqrt_eps = sqrt(eps); Jk = zeros(N);
Sk = S0; Fk = F_newton(S0); normF0 = norm(Fk,2); % Initial function evaluation and norm
normFk = normF0; alpha_2p = 1e-4; 

for k = 1:maxiters
    % Update Jacobian
    if isequal(jacobian,'analytical') % Analytical Jacobian using function handle J(S)
       Jk = J(Sk);
    elseif isequal(jacobian,'finitecolumn') % Jacobian using column-wise finite difference
        normS = norm(Sk,2);
           if normS == 0
              epsilon = sqrt_eps;
           else
              epsilon = sqrt_eps*normS;
           end
           for j = 1:N
               Jk(:,j) = (F_newton(Sk + epsilon*I(:,j)) - Fk)/epsilon;
           end
    elseif isequal(jacobian,'finitetridiag') % Jacobian exploiting tridiagonal structure
        for i = 1:3
            shift_vector = zeros(N,1);
            shift_vector((i:3:N)') = 1;
            normS = norm(Sk,2);
            norm_shift_vector = norm(shift_vector,2);
            if normS == 0
                epsilon = sqrt_eps / norm_shift_vector;
            else
                epsilon = sqrt_eps*normS / norm_shift_vector;
            end
            J_shift = (F_newton(Sk + epsilon*shift_vector) - Fk)/epsilon;

            for j = i:3:N
                if j > 1 && j < N
                    Jk(j-1:j+1, j) = J_shift(j-1:j+1);
                elseif j == 1
                    Jk(1:2, 1) = J_shift(1:2);
                elseif j == N
                    Jk(N-1:N, N) = J_shift(N-1:N);
                end
            end
        end
    end

    % Method for storing Jacobian
    % Note: if linear_method = 'full' this if statement is skipped and the
    % full Jacobian matrix will be stored
    if isequal(linear_method,'sparse') % storing the Jacobian in sparse format
        Jk = sparse(Jk);
    elseif isequal(linear_method,'sparsediag') % using tridiagonal matrix algorithm
        ind = [-1 0 1]; % position of diagonals relative to main diagonal
        diags = spdiags(Jk,ind);
        sub = diags(:,1); % sub diagonal of Jacobian
        main = diags(:,2); % main diagonal of Jacobian
        super = diags(:,3); % super diagonal of Jacobian
        Jk = spdiags([sub,main,super],[-1 0 1],N,N);
    end

    % Newton iteration
    deltaSk = -Jk\Fk; % compute Newton step
    lambda = 1; % initialise lambda
    Skt = Sk + lambda*deltaSk; % compute Newton iteration for next iterate
    Fkt = F_newton(Skt); % compute solution of new iterate
    normFkt = norm(Fkt,2); % norm of new solution

    % Line searching
    % Note: if ls = 'none' this if statement is skipped and the code
    % progresses directly to line () and accepts the Newton iteration.

    if isequal(ls,'simple') % simple line searching
        while normFkt >= normFk
            lambda = lambda/2;
            Skt = Sk + lambda*deltaSk;
            Fkt = F_newton(Skt); 
            normFkt = norm(Fkt,2);
            if lambda < 1e-4 % Stop an infinite while loop
                error('Line searching aborted as lambda is less than 1e-4.\n')  
            end
        end
    elseif isequal(ls,'two-point') % two-point line searching
        while normFkt >= (1-alpha_2p*lambda)*normFk
            g0 = normFk^2; 
            glambda = normFkt^2;
            lambda = g0*lambda^2/(glambda + (2*lambda-1)*g0);
            Skt = Sk + lambda*deltaSk;
            Fkt = F_newton(Skt); 
            normFkt = norm(Fkt,2);
            if lambda < 1e-4 % Stop an infinite while loop
                error('Line searching aborted as lambda is less than 1e-4.\n')  
            end
        end
    elseif isequal(ls,'three-point') % three-point line searching
        while normFkt >= (1-alpha_2p*lambda)*normFk
            if lambda == 1
                lambda = lambda/2;
                Skt = Sk + lambda*deltaSk;
                Fkt = F_newton(Skt); 
                normFkt = norm(Fkt,2);
                while normFkt >= (1-alpha_2p*lambda)*normFk
                    lambda = lambda/2;
                    Skt = Sk + lambda*deltaSk;
                    Fkt = F_newton(Skt); 
                    normFkt = norm(Fkt,2);
                end
                lambda1 = lambda;
            end
            g0 = normFk^2;

            Sk_lambda1 = Sk + lambda1*deltaSk;
            F_lambda1 = F_newton(Sk_lambda1);
            normF_lambda1 = norm(F_lambda1,2);
            g_lambda1 = normF_lambda1^2;

            g_lambda2 = normFkt^2;

            C1 = (g_lambda1-g0)/((lambda1-lambda)*lambda1);
            C2 = (g_lambda2-g0)/((lambda1-lambda)*lambda);

            lambda_new = (lambda1*C2-lambda*C1)/(2*C2-2*C1);

            lambda = lambda_new;
            Skt = Sk + lambda*deltaSk;
            Fkt = F_newton(Skt); 
            normFkt = norm(Fkt,2);

            if lambda < 1e-4 % Stop an infinite while loop
                error('Line searching aborted as lambda is less than 1e-4.\n')  
            end
        end
    end

    % Accept Newton iteration
    Sk = Skt;
    Fk = Fkt; 
    normFk = normFkt;
    
    % Stopping criterion
    if normFk < rtol*normF0 + atol
        S = Sk;
        break;
    end

end     

end

end