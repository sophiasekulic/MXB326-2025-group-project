%% MXB326 Group Project: Group 10
% Sheren Zein (n10818120), Sophia Sekulic (n10755861), Zach Eyre (n10818189)
clc, close, clear all

%% ================== 1.0 Numerical Solution ==================
%% 1.1 Proof that boundary conditions (9) and (19) are equivalent

% Let S1 = S(0,t)

% (9): DS/dx(S1) = alpha/beta*(beta*S1+gamma) + omega/beta*(beta*S1+gamma)^2

% (19): g(S1)*DS/dx(S1) + f(S1) = 1
% Rearrange (19): DS/dx(S1) = (1 - f(S1))/g(S1)

% Statements will be equivalent if RHS of (9) = RHS of (19)

syms S1 beta Swr Sor F_visc

gamma = beta*(1-Swr-F_visc*Sor)/(F_visc-1);
alpha = -(beta^2)*((Sor+gamma/beta)*(1-Swr+gamma/beta))/(1-Swr-Sor);
omega = beta - alpha/(beta-beta*Swr+gamma);

f = @(S) alpha/(beta^2)*(1/(1-Swr+gamma/beta) - 1/(S+gamma/beta));
g = @(S) 1/((beta*S+gamma)^2);

RHS9 = alpha/beta*(beta*S1+gamma) + omega/beta*(beta*S1+gamma)^2;
RHS19 = (1 - f(S1))/g(S1);

BCequivalence = isequal(simplify(RHS9),simplify(RHS19))

% Parameters
Swr = 0.0375; % minimum saturation of water
Sor = 0.15; % minimum saturation of oil
F_visc = 2; % water-to-oil viscosity ratio          
beta = [1,2,4,10]; % relative magnitude of viscous to capillary terms
theta = 1; % theta method weighting parameter (theta = 0, 1/2 or 1)
sigma = 1/2; % weighted average approximation parameter (sigma = 0, 1/2 or 1)
L = 3; % length of domain
N = 200; % number of nodes
x = linspace(0,L,N)'; % node locations
M = 800; % number of time steps
T = 2; % end time
BC = 'zerogradient'; % boundary condition at x = L, 'zeroflux' or 'zerogradient'
jacobian = 'analytical'; % Jacobian generation options 'analytical','finitecolumn','finitetridiag'
ls = 'none'; % line searching options 'simple', 'two-point' or 'three-point'
linear_method = 'sparse'; % method for storing Jacobian 'full', 'sparse' or 'sparsediag'
atol = 1e-10; rtol = 1e-10; maxiters = 100; % tolerence for Newton's method and max iterations

%% ================== 2.0 Analysis ==================
%% 2.1 Comparing solutions at various given times
close all
dt = T/M; % time step
times = [0.01,0.05,0.1,0.3,0.8125,1]; % times to compare solution at
index = times/dt+1; % column index corresponding to previously specified times

% 2.1.1 Solutions for beta = 1

% Semi-analytical solution
Sa_1 = semi_analytical_solution(Swr,Sor,beta(1),F_visc,L,N,M,T);
% Numerical solution
Sn_1 = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
% Plot solutions
plotFunc(length(times),Sn_1,Sa_1,x,index,beta(1))  

% 2.1.2 Solutions for beta = 2

% Semi-analytical solution
Sa_2 = semi_analytical_solution(Swr,Sor,beta(2),F_visc,L,N,M,T);
% Numerical solution
Sn_2 = numerical_solution(Swr,Sor,F_visc,beta(2),theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
% Plot solutions
plotFunc(length(times),Sn_2,Sa_2,x,index,beta(2))  

% 2.1.3 Solutions for beta = 4

% Analytical solution
Sa_4 = semi_analytical_solution(Swr,Sor,beta(3),F_visc,L,N,M,T);
% Numerical solution
Sn_4 = numerical_solution(Swr,Sor,F_visc,beta(3),theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
% Plot solutions
plotFunc(length(times),Sn_4,Sa_4,x,index,beta(3)) 

% 2.1.4 Solutions for beta = 10

% Analytical solution
Sa_10 = semi_analytical_solution(Swr,Sor,beta(4),F_visc,L,N,M,T);
% Numerical solution
Sn_10 = numerical_solution(Swr,Sor,F_visc,beta(4),theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
% Plot solutions
plotFunc(length(times),Sn_10,Sa_10,x,index,beta(4))  

pause; close all

% 2.1.5 Comparing boundary conditions for the numerical solution at x = L

comparingBC(Swr,Sor,beta,F_visc,L,N,M,3,theta,sigma,jacobian,ls,linear_method,atol,rtol,maxiters)

pause; close all

% 2.2 Average reservoir oil saturation

% Comparison for beta = 1
plotAverage(beta(1),Swr,Sor,F_visc,theta,sigma,L,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)

% Comparison for beta = 2
plotAverage(beta(2),Swr,Sor,F_visc,theta,sigma,L,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)
 
% Comparison for beta = 4
plotAverage(beta(3),Swr,Sor,F_visc,theta,sigma,L,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)

% Comparison for beta = 10
plotAverage(beta(4),Swr,Sor,F_visc,theta,sigma,L,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)

pause; close all

%% 2.3 Efficiency & Accuracy Comparison

% Testing efficiency and accuracy of the different numerical strategies
% Default parameters: 
Swr = 0.0375; Sor = 0.15; F_visc = 2; beta = beta(1); 
theta = 1; % backward euler
sigma = 0.5; % averaging
N = 100; M = 400; T = 1; x = linspace(0,L,N)'; 
BC = 'zerogradient'; 
jacobian = 'analytical'; % analytical jacobian
ls = 'none'; % no line searching
linear_method = 'sparse'; % sparse matrix storage for jacobian
atol = 1e-3; rtol = 1e-3; maxiters = 100;
tol = 1e-3; % Acceptable relative tolerance 

%% 2.3.1 Semi-analytical solution timing

Sa_solve = @() semi_analytical_solution(Swr,Sor,beta(1),F_visc,L,N,M,T);
Sa_time = timeit(Sa_solve);
Sa_true = semi_analytical_solution(Swr,Sor,beta(1),F_visc,L,N,M,T); % true solution for comparison

%% 2.3.2 Numerical solution timing

%% 2.3.2.1 Theta method

% Backward Euler: theta = 1
Sn_solve_be = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
Sn_time_be = timeit(Sn_solve_be);
Sn_be = numerical_solution(Swr,Sor,F_visc,beta(1),1,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
% Crank Nicolson: theta = 0.5
Sn_solve_cn = @() numerical_solution(Swr,Sor,F_visc,beta(1),0.5,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
Sn_time_cn = timeit(Sn_solve_cn);
Sn_cn = numerical_solution(Swr,Sor,F_visc,beta(1),0.5,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);

% Compare efficiency and accuracy
methods = {'Semi-analytical','Backward Euler','Crank Nicolson'}; 
methods1 = {'Backward Euler','Crank Nicolson'}; 
times = [Sa_time,Sn_time_be,Sn_time_cn];
solutions = {Sn_be, Sn_cn};
comparison(Sa_true,methods,methods1,times,solutions,tol)

% Crank Nicolson has higher accuracy 
theta = 0.5;

%% 2.3.2.2 Flux Method

% Upwinding: sigma = 0
Sn_solve_up = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,0,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
Sn_time_up = timeit(Sn_solve_up);
Sn_up = numerical_solution(Swr,Sor,F_visc,beta(1),theta,0,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
% Averaging: sigma = 0.5
Sn_solve_av = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,0.5,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
Sn_time_av = timeit(Sn_solve_av);
Sn_av = numerical_solution(Swr,Sor,F_visc,beta(1),theta,0.5,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);

% Compare efficiency and accuracy
methods_flux = {'Semi-analytical','Upwinding','Averaging'}; 
methods_flux1 = {'Upwinding','Averaging'}; 
times_flux = [Sa_time,Sn_time_up,Sn_time_av];
solutions_flux = {Sn_up, Sn_av};
comparison(Sa_true,methods_flux,methods_flux1,times_flux,solutions_flux,tol)

% Averaging is within acceptable relative error and is the default

%% 2.3.2.3. Jacobian matrix generation methods

% Analytical Jacobian 
Sn_solve_aj = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,'analytical',ls,linear_method,atol,rtol,maxiters);
Sn_time_aj = timeit(Sn_solve_aj);
Sn_aj = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,'analytical',ls,linear_method,atol,rtol,maxiters);
% Finite Column Jacobian 
Sn_solve_fcj = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,'finitecolumn',ls,linear_method,atol,rtol,maxiters);
Sn_time_fcj = timeit(Sn_solve_fcj);
Sn_fcj = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,'finitecolumn',ls,linear_method,atol,rtol,maxiters);
% Tridiagonal FD Jacobian 
Sn_solve_tfdj = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,'finitetridiag',ls,linear_method,atol,rtol,maxiters);
Sn_time_tfdj = timeit(Sn_solve_tfdj);
Sn_tfdj = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,'finitetridiag',ls,linear_method,atol,rtol,maxiters);

% Compare efficiency and accuracy
methods_jac = {'Semi-analytical','Analytical Jacobian','Finite Column Jacobian', ...
    'Tridiagonal FD Jacobian'}; 
methods_jac1 = {'Analytical Jacobian','Finite Column Jacobian', ...
    'Tridiagonal FD Jacobian'}; 
times_jac = [Sa_time,Sn_time_aj,Sn_time_fcj,Sn_time_tfdj];
solutions_jac = {Sn_aj,Sn_fcj,Sn_tfdj};
comparison(Sa_true,methods_jac,methods_jac1,times_jac,solutions_jac,tol)

% Analytical jacobian method has highest efficiency and is the default 

%% 2.3.2.4 Line Searching Methods

% Simple line searching
Sn_solve_sls = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,'simple',linear_method,atol,rtol,maxiters);
Sn_time_sls = timeit(Sn_solve_sls);
Sn_sls = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,'simple',linear_method,atol,rtol,maxiters);
% Two-point line searching
Sn_solve_2ls = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,'two-point',linear_method,atol,rtol,maxiters);
Sn_time_2ls = timeit(Sn_solve_2ls);
Sn_2ls = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,'two-point',linear_method,atol,rtol,maxiters);
% Three-point line searching
Sn_solve_3ls = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,'three-point',linear_method,atol,rtol,maxiters);
Sn_time_3ls = timeit(Sn_solve_3ls);
Sn_3ls = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,'three-point',linear_method,atol,rtol,maxiters);

% Compare efficiency and accuracy
methods_ls = {'Semi-analytical','Simple Line-searching','Two-point Line-searching', ...
    'Three-point Line-searching'}; 
methods_ls1 = {'Simple Line-searching','Two-point Line-searching', ...
    'Three-point Line-searching'}; 
times_ls = [Sa_time,Sn_time_sls,Sn_time_2ls,Sn_time_3ls];
solutions_ls = {Sn_sls,Sn_2ls,Sn_3ls};
comparison(Sa_true,methods_ls,methods_ls1,times_ls,solutions_ls,tol)

% Three-point line searching is the most efficient
ls = 'three-point';

%% 2.3.2.5 Jacobian matrix storage methods

% Full matrix
Sn_solve_full = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,'full',atol,rtol,maxiters);
Sn_time_full = timeit(Sn_solve_full);
Sn_full = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,'full',atol,rtol,maxiters);
% Sparse matrix 
Sn_solve_sparse = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,'sparse',atol,rtol,maxiters);
Sn_time_sparse = timeit(Sn_solve_sparse);
Sn_sparse = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,'sparse',atol,rtol,maxiters);
% Sparse tridiagonal matrix
Sn_solve_tri = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,'sparsediag',atol,rtol,maxiters);
Sn_time_tri = timeit(Sn_solve_tri);
Sn_tri = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,'sparsediag',atol,rtol,maxiters);

% Compare efficiency and accuracy
methods_jacstor = {'Semi-analytical','Full Jacobian','Sparse Jacobian','Sparse Tridiagonal Jacobian'}; 
methods_jacstor1 = {'Full Jacobian','Sparse Jacobian','Sparse Tridiagonal Jacobian'};
times_jacstor = [Sa_time,Sn_time_full,Sn_time_sparse,Sn_time_tri];
solutions_jacstor = {Sn_full,Sn_sparse,Sn_tri};
comparison(Sa_true,methods_jacstor,methods_jacstor1,times_jacstor,solutions_jacstor,tol)

% Spares jacobian is best and is the default

%% 2.3.3 Ideal numerical strategy

% The numerical strategy that gives the most efficient and accurate
% solution is one that uses the Crank Nicolson theta method, average flux,
% analytical jacobian, three-point line searching and sparse storage

dt = T/M; % time step
times = [0.01,0.05,0.1,0.3,0.8125,1]; % times to compare solution at
index = times/dt+1; % column index corresponding to previously specified times

Sn_ideal = numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
Sn_solve_ideal = @() numerical_solution(Swr,Sor,F_visc,beta(1),theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
Sn_time_ideal = timeit(Sn_solve_ideal);

% Compare ideal numerical solution visually and computationally
figure;
for n = 1:length(times)
    plot(x,Sn_ideal(:,index(n)),'b-','LineWidth',2.5)
    hold on
    plot(x,Sa_true(:,index(n)),'r--','LineWidth',2.5)
    hold on
    ylim([0,1])
    xlabel('x','FontSize',24,'interpreter','latex')
    ylabel('S(x,t)','FontSize',24,'interpreter','latex')
    set(gca,'FontSize',20)
    legend('Numerical','Semi-Analytical','Location','southeast')
    title(['Semi-analytical vs Ideal Numerical Solution for $\beta$ = 1'], ...
        'Interpreter','latex','FontSize',18);
    text(0.1,0.07,'t = 0.01,0.05,0.1,0.3,0.8125,1 respectively','Interpreter','latex',...
        'FontSize',12)
    drawnow
end
methods_ideal = {'Semi-analytical','Ideal numerical'}; 
methods_ideal1 = {'Ideal numerical'};
times_ideal = [Sa_time,Sn_time_ideal];
solutions_ideal = {Sn_ideal};
comparison(Sa_true,methods_ideal,methods_ideal1,times_ideal,solutions_ideal,tol)

%% 2.4 Comparisons to Buckley-Leverett solution

% For high water injection rates of velocity v, the solution to S should
% resemble a step function, with a sudden change at x=vt (the current x position
% of the shock front)

% Testing for large value of beta
beta = [10 15 20]; T = 1;

BuckleyLeverett_comparison(Swr,Sor,beta(1),F_visc,L,N,M,T,theta,sigma,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)
pause; close;

BuckleyLeverett_comparison(Swr,Sor,beta(2),F_visc,L,N,M,T,theta,sigma,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)
pause; close;
 
BuckleyLeverett_comparison(Swr,Sor,beta(3),F_visc,L,N,M,T,theta,sigma,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)
pause; close;

% Displaying how different beta values look at when compared to the
% quickest breakthrough time for semi-analytical models
S1 = semi_analytical_solution(Swr,Sor,1,F_visc,L,N,M,T);
S10 = semi_analytical_solution(Swr,Sor,10,F_visc,L,N,M,T);

% Defining bore position better in case L or N changes
[~, bore_idx] = min(abs(x - 1)); % position of x=1 in x array
bore_xPosition = x(bore_idx); % sanity check for x value
bt_timing = round(find(S1 < 0.999-Swr,1, 'first') - 1) * dt); % timing for breakthrough with beta = 1

figure
for i = bt_timing
    hold on
    plot(x, S10(:,i),'r-','LineWidth',3)
    plot(x, S1(:,i),'b-','LineWidth',3)
    xline(1, 'k--',LineWidth=2)
    hold off
    ylim([0,1])
    xlim([0,1.5])
    xlabel('x','FontSize',24,'interpreter','latex')
    ylabel('S(x,t)','FontSize',24,'interpreter','latex')
    legend('\beta=10','\beta=1','Location','southeast')
    title({'Impact of beta value on breakthrough time', 'at t = 0.1325'}, ...
        'Interpreter','latex','FontSize',20);
end
pause; close;

%% 2.5 Breakthrough Timings

beta = [1,2,4,10]; % Beta values being examined  
    
bt_tol = 0.999-Swr; % saturation tolerance for a breakthrough
bt_timings = zeros(2,4); % array for tracking breakthrough times
    
% Defining bore position better in case L or N changes
[~, bore_idx] = min(abs(x - 1)); % position of x=1 in x array
bore_xPosition = x(bore_idx); % sanity check for x value
    
% Breakthrough checking
for n = 1:length(beta)
    S_analytical = semi_analytical_solution(Swr,Sor,beta(n),F_visc,L,N,M,T);
    S_numerical = numerical_solution(Swr,Sor,F_visc,beta(n),theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
        
    % Only care about values for the node at x=1
    S_analytical = S_analytical(bore_idx,:);
    S_numerical = S_numerical(bore_idx,:);
        
    % Find earliest change at well-bore then converting that to a time
    % in seconds
    Sa_btTime = (find(S_analytical < bt_tol,1, 'first') - 1) * T/M;
    Sn_btTime = (find(S_numerical < bt_tol,1,  'first') - 1) * T/M;

    bt_timings(:,n) = [Sa_btTime, Sn_btTime]';
end

method_comp = {'Beta value';'Semi-Analytical'; 'Numerical'};
betaVal = {'β = 1', 'β = 2', 'β = 4', 'β = 10'}; 

% Displaying table
comparison_results = table(method_comp, [betaVal; num2cell(bt_timings)], ...
    'VariableNames', {'Method', 'Breakthrough time (seconds)'});
disp(comparison_results)

% Breakthrough timings approaches the expected T_B = 1-Swr-Sor as beta values increase

%% ================== 3.0 Functions ==================
%% 3.1 Plotting Beta Solutions Function
function plotFunc(iters,S_numerical,S_analytical,x,index,beta)
figure;
for n = 1:iters
    plot(x,S_numerical(:,index(n)),'b-','LineWidth',2.5)
    hold on
    plot(x,S_analytical(:,index(n)),'r--','LineWidth',2.5)
    hold on
    ylim([0,1])
    xlabel('x','FontSize',24,'interpreter','latex')
    ylabel('S(x,t)','FontSize',24,'interpreter','latex')
    set(gca,'FontSize',20)
    legend('Numerical','Semi-Analytical','Location','southeast')
    title(['Semi-analytical vs Numerical Solutions for $\beta$ = ', num2str(beta)], ...
        'Interpreter','latex','FontSize',18);
    text(0.1,0.07,'t = 0.01,0.05,0.1,0.3,0.8125,1 respectively','Interpreter','latex',...
        'FontSize',12)
    drawnow
end
end

%% 3.1.2 Plotting Solutions with Different Boundary Conditions
function comparingBC(Swr,Sor,beta,F_visc,L,N,M,T,theta,sigma,jacobian,ls,linear_method,atol,rtol,maxiters)
x = linspace(0,L,N)';
% Semi-analytical solution
Sa = semi_analytical_solution(Swr,Sor,beta(2),F_visc,L,N,M,T);
% Numerical solution: Zero flux
Sn_zeroflux = numerical_solution(Swr,Sor,F_visc,beta(2),theta,sigma,N,M,T,x,'zeroflux',jacobian,ls,linear_method,atol,rtol,maxiters);
% Numerical solution: Zero gradient
Sn_zerogradient = numerical_solution(Swr,Sor,F_visc,beta(2),theta,sigma,N,M,T,x,'zerogradient',jacobian,ls,linear_method,atol,rtol,maxiters);

figure;
for n = M-50:M+1
    plot(x,Sn_zeroflux(:,n),'b-','LineWidth',3)
    hold on
    plot(x,Sn_zerogradient(:,n),'c-','LineWidth',3)
    hold on
    plot(x,Sa(:,n),'r-','LineWidth',3)
    hold off
    ylim([0.1,0.3])
    xlim([1,3])
    xlabel('x','FontSize',24,'interpreter','latex')
    ylabel('S(x,t)','FontSize',24,'interpreter','latex')
    set(gca,'FontSize',18)
    legend('Numerical: Zero Flux BC','Numerical: Zero Gradient BC','Semi-Analytical','Location','northwest')
    title('Semi-analytical vs Numerical Solutions at x = L','Interpreter','latex','FontSize',20)
    drawnow
end
end

%% 3.2 Plotting average oil saturations function
function plotAverage(beta,Swr,Sor,F_visc,theta,sigma,L,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)

dt = T/M;
t = 0:dt:T;

% Analytical average
avg_analyticalfunction = @(t) -t/L + 1 - Swr;
avg_sat_analytical = avg_analyticalfunction(t);

% Numerical Average
S = numerical_solution(Swr,Sor,F_visc,beta,theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);
% Trapezoidal approximation
h = L/(N-1);
avg_sat_numerical = zeros(1,M+1);
for n = 1:M+1
    Sk = S(:,n);
    I = h/2*(Sk(1) + 2*sum(Sk(2:N-1)) + Sk(N)); % integral approximation
    avg_sat_numerical(n) = 1/L * I;
end

figure;
plot(t,avg_sat_numerical,'b-','LineWidth',3)
hold on
plot(t,avg_sat_analytical,'r--','LineWidth',3)
hold off
ylim([0,1])
xlim([0,T])
xlabel('t','FontSize',24,'interpreter','latex')
ylabel('Average Reservoir Oil Saturation $\bar{S}(t)$','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14)
legend('Numerical solution','Analytical solution','Location','southwest','FontSize',18)
title(['Evolution of Average Reservoir Oil Saturation for $\beta$ = ',num2str(beta),''],...
    'Interpreter','latex','FontSize',20)
pause
end

%% 3.3 Efficiency and accuracy function
function comparison(true_Sa,methods,methods1,times,solutions,tol)

efficiency_results = table(methods',times','VariableNames',{'Method','Time'});
disp(efficiency_results)

% Find fastest method
num_methods = methods(2:end);
num_times = times(2:end);
[sort_times,sort_idx] = sort(num_times);
fastest_method = num_methods(sort_idx);
for i = 1:length(num_times)
    fprintf('%d. %s - %.4f seconds\n',i,fastest_method{i},sort_times(i));
end

rel_error = zeros(length(solutions),1);
% Computer relative errors of each solution
for i = 1:length(solutions)
    Sn_i = solutions{i};
    rel_error(i) = norm(true_Sa-Sn_i,2)/norm(true_Sa,2);
end
% Display relative error table
acceptable = rel_error < tol;
accuracy_results = table(methods1(:),rel_error,acceptable,...
    'VariableNames',{'Method','Relative Error','Acceptable'});
disp(accuracy_results)
end

%% 3.4 Comparison of solutions to Buckley-Leverett solution for large beta values
function BuckleyLeverett_comparison(Swr,Sor,beta,F_visc,L,N,M,T,theta,sigma,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters)

S_analytical = semi_analytical_solution(Swr,Sor,beta,F_visc,L,N,M,T);
S_numerical = numerical_solution(Swr,Sor,F_visc,beta,theta,sigma,N,M,T,x,BC,jacobian,ls,linear_method,atol,rtol,maxiters);

dt = T/M;
t = 0:dt:M;

v = 1/(1-Swr-Sor); % Water injection rate v for large beta values
shockFront_position = v * t; % Position of shock front at any given time t

% Prefilling S_buckLev with 1-Swr unless shock front has reached that node
S_buckLev = ones(N,M+1) * (1-Swr);
for i = 1:length(t)
    SF_mask = x < shockFront_position(:,i); 
    S_buckLev(SF_mask==1,i) = Sor; % adjusts any node behind shock front
end

% Relative error of solutions vs Buckley-Leverett solution
relErr_a = zeros(1, M+1); % BL vs Semi-Analytical
relErr_n = zeros(1, M+1); % BL vs Numerical

% Avoids the legend being ugly and blocking the bottom right of the graph
% during the figure animation, will swap position 3/4 through animation
legendSwap_idx = length(t) * (3/4); 

for i = 1:10:M+1
    plot(x, S_numerical(:,i),'b-','LineWidth',2.5)
    hold on
    plot(x, S_analytical(:,i),'r-','LineWidth',2.5)
    plot(x, S_buckLev(:,i),'g--','LineWidth',2.5)
    hold off
    ylim([0,1])
    xlabel('x','FontSize',18,'interpreter','latex')
    ylabel('S(x,t)','FontSize',18,'interpreter','latex')
    if i < legendSwap_idx
        legend('Numerical','Semi-Analytical','Buckley-Leverett','Location','southeast','FontSize',16)
    else
        legend('Numerical','Semi-Analytical','Buckley-Leverett','Location','northwest','FontSize',16)
    end
    title({'Comparison of Semi-analytical, Numerical and'; ...
        ['Buckley-Leverett solutions for $\beta$ = ',num2str(beta),'']}, ...
        'Interpreter','latex', 'FontSize',16)
    drawnow

    relErr_a(i) = norm(S_analytical(:,i) - S_buckLev(:,i), 2) / norm(S_buckLev(:,i), 2);
    relErr_n(i) = norm(S_numerical(:,i) - S_buckLev(:,i), 2) / norm(S_buckLev(:,i), 2);
end

% Average relative errors over time
avg_rel_a = mean(relErr_a);
avg_rel_n = mean(relErr_n);

methods_comp = {'Semi-Analytical'; 'Numerical'};
avg_errors = [avg_rel_a; avg_rel_n];

% Displaying table
comparison_results = table(methods_comp, avg_errors, ...
    'VariableNames', {'Method', 'Average Relative Error'});
disp(comparison_results)
end
