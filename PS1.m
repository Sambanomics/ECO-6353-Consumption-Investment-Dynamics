
stoch_simul c y inv w r a n;

% Outlining the IRFs:
irf_y = oo_.irfs.y_e;
irf_c = oo_.irfs.c_e;
irf_inv = oo_.irfs.inv_e;
irf_n = oo_.irfs.n_e;
irf_r = oo_.irfs.r_e;
irf_w = oo_.irfs.w_e;
irf_a = oo_.irfs.a_e;

% Setting the plot to 10 periods:
nPeriods = 10;


figure;
tiledlayout(3,3)

% Ploting for Output(y):
nexttile;
plot(1:nPeriods, irf_y, 'LineWidth', 2);
title('IRF of Output (y)');
xlabel('Periods');
ylabel('Deviation');

% Ploting for Consumption(c):
nexttile;
plot(1:nPeriods, irf_c, 'LineWidth', 2);
title('IRF of Consumption (c)');
xlabel('Periods');
ylabel('Deviation');

% Ploting for Investment(inv):
nexttile;
plot(1:nPeriods, irf_inv, 'LineWidth', 2)
title('IRF of Investment (inv)');
xlabel('Periods');
ylabel('Deviation');

% Ploting for Labour(n):
nexttile;
plot(1:nPeriods, irf_n, 'LineWidth', 2)
title('IRF of Labour (n)');
xlabel('Periods');
ylabel('Deviation');

% Ploting for Interest(r):
nexttile;
plot(1:nPeriods, irf_r, 'LineWidth', 2)
title('IRF of Interest (r)');
xlabel('Periods');
ylabel('Deviation');

% Ploting for Wage(w):
nexttile;
plot(1:nPeriods, irf_c, 'LineWidth', 2)
title('IRF of Wage (w)');
xlabel('Periods');
ylabel('Deviation');

% Ploting for Technology(a):
nexttile;
plot(1:nPeriods, irf_c, 'LineWidth', 2)
title('IRF of Technology (a)');
xlabel('Periods');
ylabel('Deviation');

% Outlining the IRFs:
y_sequence = oo.endo_simul(strmatch('y', M_.endo_names, 'exact'), :);
c_sequence = oo.endo_simul(strmatch('c', M_.endo_names, 'exact'), :);
inv_sequence = oo.endo_simul(strmatch('inv', M_.endo_names, 'exact'), :);
n_sequence = oo.endo_simul(strmatch('n', M_.endo_names, 'exact'), :);
r_sequence = oo.endo_simul(strmatch('r', M_.endo_names, 'exact'), :);
w_sequence = oo.endo_simul(strmatch('w', M_.endo_names, 'exact'), :);
a_sequence = oo.endo_simul(strmatch('a', M_.endo_names, 'exact'), :);

% Setting number of Lags;
nLags = 6;

figure;
tiledlayout(3,3); 

% Ploting for Output(y):
nexttile;
autocorr(y_sequence, nLags);
title('ACF of Output (y)');

% Ploting for Consumption(c):
nexttile;
autocorr(c_sequence, nLags);
title('ACF of Consumption (c)');

% Ploting for Investment(inv):
nexttile;
autocorr(inv_sequence, nLags);
title('ACF of Investment (inv)');

% Ploting for Labour(n):
nexttile;
autocorr(n_sequence, nLags);
title('ACF of Labour (n)');

% Ploting for Interest(r):
nexttile;
autocorr(r_sequence, nLags);
title('ACF of Interest (r)');

% Ploting for Wage(w):
nexttile;
autocorr(w_sequence, nLags);
title('ACF of Wage (w)');

% Ploting for Technology(a):
nexttile;
autocorr(a_sequence, nLags);
title('ACF of Technology (a)');
