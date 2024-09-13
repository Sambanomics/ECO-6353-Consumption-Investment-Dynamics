
%% Stating the Endogenous Variable: 



var c w n r a k inv y pi_f rk;

%% Stating the Exogenous Variable:  



varexo e;



%% Stating the Parameters:

parameters gamma_c beta_p delta_k phi_u alpha_y psi_n rho_t sigma_t n_target;



%% Calibration of the Parameters:

beta_p = 0.99;

gamma_c = 1;

phi_u = 2;

delta_k = 0.025;

alpha_y = 0.35;

rho_t = 0.92;

sigma_t = 0.04;

n_target = 1/3;



%% Outlining the Model:


model;

% Labour Supply Condition:

c^-gamma_c * w = psi_n * n^phi_u;

% Euler Condition:

c^-gamma_c = beta_p * (1+r) * c(+1)^-gamma_c;

% Capital Return:

rk = delta_k + r;

% Capital Law of Motion:

k = (1 - delta_k) * k(-1) + inv;

% Factor Price of Capital:

rk = alpha_y / k(-1);

% Factor Price of Labour:

w = (1 - alpha_y) * y/n;

% Production Function:

y = a * k(-1)^alpha_y * n^(1 - alpha_y);

% Output:
y = c + inv;

% Profit Level:
pi_f = y - w * n - r * k;

% Technology shock process:

log(a) = rho_t * log(a(-1)) + e;

end;



%% Initial Steady State Values:


initval;

% Productivity at steady state following a log AR(1) process:   
a = 1;  

% Labour; 38 hours per week; assuming 8hours per day approximately = 8/24
n = 1/3;  

% Guess for Capital Stock:
k = 10;       
    
% Guess for Output:
y = 1;   
    
% Guess for Consumption:
c = 0.8;     
    
% Guess for Investment:
inv = 0.1;

% Guess for Rental Rate of Capital:
rk = 0.05; 

% Guess for Wage:
w = 1;  

% Guess for Interest Rate:
r = 0.04;  

% Guess for Profit:
pi_f = 0.01; 

    end;


shocks;

var e; stderr 0.04;

end;

psi_n = (c^-gamma_c * w) / n^phi_u ;



% Stochastic Stimulation Process:

stoch_simul c y inv w r a n ;







