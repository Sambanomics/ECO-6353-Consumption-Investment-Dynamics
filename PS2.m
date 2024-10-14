
%% Outlining the IRFs (Technology Shock Process)

irf_y = oo_.irfs.y_e_a;
irf_c = oo_.irfs.c_e_a;
irf_inv = oo_.irfs.inv_e_a;
irf_n = oo_.irfs.n_e_a;
irf_r = oo_.irfs.r_e_a;
irf_w = oo_.irfs.w_e_a;
irf_a = oo_.irfs.a_e_a;

% Setting the plot to 40 periods:
nPeriods = 40;


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





%% Outlining the IRFs (Monetary Shock Process)

irf_y = oo_.irfs.y_e_m;
irf_c = oo_.irfs.c_e_m;
irf_inv = oo_.irfs.inv_e_m;
irf_n = oo_.irfs.n_e_m;
irf_r = oo_.irfs.r_e_m;
irf_w = oo_.irfs.w_e_m;
irf_a = oo_.irfs.a_e_m;

% Setting the plot to 40 periods:
nPeriods = 40;


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




