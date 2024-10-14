
%% Stating the Endogenous Variable: 

var c n rk r u q inv k w y a pi i mc m;

%% Stating the Exogenous Variable:  

varexo e_a e_m;

%% Stating the Parameters:

parameters betta lamda gama phhi pssi deltta alph rho_a rho_m sigma_a sigma_m zeta_pi zeta_y xi_p xi_w rho_i eta_w adj_ii kappa_p kappa_w omega;

%% Calibration of the Parameters:

betta = 0.99;           
lamda = 0.75;           
gama = 1;               
phhi = 1.25;            
deltta = 0.025;         
alph = 0.35;            
rho_a = 0.92;           
rho_m = 0.4;            
sigma_a = 0.04;         
sigma_m = 0.1;          
zeta_pi = 1.5;          
zeta_y = 0.1;          
xi_p = 0.6;             
xi_w = 0.9;             
rho_i = 0.78;           
eta_w = 2.1;            
adj_ii = 5;             
omega = 1.6;            
kappa_p = (1 - xi_p) * (1 - betta * xi_p) / xi_p;  
kappa_w = (1 - xi_w) * (1 - betta * xi_w) / xi_w;  

%% Outlining the Model:

model(linear);

% Log-linearized Euler equation

c = lamda * c(-1) + (c(+1) - lamda * c) - (1 /gama) * r;
    
% Labor Supply equation

w = gama * (c - lamda * c(-1)) + phhi * n;
    
% Capital accumulation equation (log-linearized)

k = (1 - deltta) * k(-1) + deltta * inv;

% Investment (Tobin's q)

q = adj_ii*(inv - inv(-1)) - betta*adj_ii*(inv(+1) - inv);

% Wage equation

w = y - n;

% Capital rental rate

rk = y - k(-1);

% Output function

y  = a + alph*u + alph*k(-1) + (1-alph)*n;

% Interest rule

r = rk + deltta; 

% Goods market clearing condition

y = c + inv + u + k(-1);

% Monetary policy rule

i = rho_i * (i(-1) - i) + (1-rho_i) * (zeta_pi * pi + zeta_y *( y - y(-1))) + m;

% New Keynesian Philips curve (NKPC)

pi = kappa_p * mc + betta * pi(+1);

% Marginal cost

mc = (w - a) - (1-alph) * y;

% Wage setting

w = betta*w(+1) + kappa_w * (w - phhi * n);

% Fisher equation

i = r + pi;

% Technology shock process

 a = rho_a * a(-1) + e_a;

% Monetary shock process

 m = rho_m * m(-1) + e_m;

end;


% Initial values

initval;
    c = 0;
    k = 0;
    inv = 0;
    n = 0;
    pi = 0;
    w = 0;
    r = 0;
    y = 0;
    mc = 0;
    a = 0;
    m = 0;

end;

% Outlining the Shocks

shocks;

    var e_a; stderr sigma_a;
    
    var e_m; stderr sigma_m;

end;


% Stochastic Stimulation Process:

stoch_simul c y inv w r a n;