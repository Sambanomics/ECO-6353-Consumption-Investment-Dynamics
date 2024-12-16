%Extracting irfs for all observables at the posterior means, medians, and
%Mode

%Technology shock

%Posterior mean
dy_to_ea = oo_.PosteriorIRF.dsge.Mean.dy_ea;
dc_to_ea = oo_.PosteriorIRF.dsge.Mean.dc_ea;
dinve_to_ea = oo_.PosteriorIRF.dsge.Mean.dinve_ea;
dw_to_ea = oo_.PosteriorIRF.dsge.Mean.dw_ea;
pinfobs_to_ea = oo_.PosteriorIRF.dsge.Mean.pinfobs_ea;
robs_to_ea = oo_.PosteriorIRF.dsge.Mean.robs_ea;
labobs_to_ea = oo_.PosteriorIRF.dsge.Mean.labobs_ea;

%Posterior Median
dy_to_ea_median = oo_.PosteriorIRF.dsge.Median.dy_ea;
dc_to_ea_median = oo_.PosteriorIRF.dsge.Median.dc_ea;
dinve_to_ea_median = oo_.PosteriorIRF.dsge.Median.dinve_ea;
dw_to_ea_median = oo_.PosteriorIRF.dsge.Median.dw_ea;
pinfobs_to_ea_median = oo_.PosteriorIRF.dsge.Median.pinfobs_ea;
robs_to_ea_median = oo_.PosteriorIRF.dsge.Median.robs_ea;
labobs_to_ea_median = oo_.PosteriorIRF.dsge.Median.labobs_ea;




% Number of periods to plot
nPeriods = 40;

% Creating a figure using tiledlayout
figure;
tiledlayout(2, 7);  % 4x2 grid for 7 variables

% Ploting IRF for output (y)
nexttile;
plot(1:nPeriods, dy_to_ea, 'LineWidth', 2);
title('IRF of Output (dy)'); 
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for consumption
nexttile;
plot(1:nPeriods, dc_to_ea, 'LineWidth', 2);
title('IRF of consumption (dc)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for Investment
nexttile;
plot(1:nPeriods, dinve_to_ea, 'LineWidth', 2);
title('IRF of Investment (dinve)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for wage
nexttile;
plot(1:nPeriods, dw_to_ea, 'LineWidth', 2);
title('IRF of wage (dw)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Inflation rate
nexttile;
plot(1:nPeriods, pinfobs_to_ea, 'LineWidth', 2);
title('IRF of Inflation rate (pinfobs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Federal funds rate
nexttile;
plot(1:nPeriods, dc_to_ea, 'LineWidth', 2);
title('IRF of federal funds rate (robs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for labor hours
nexttile;
plot(1:nPeriods, dc_to_ea, 'LineWidth', 2);
title('IRF of labor hours (labobs)');
xlabel('Periods');
ylabel('Deviation');



% Ploting IRF for output (y)
nexttile;
plot(1:nPeriods, dy_to_ea_median, 'LineWidth', 2);
title('IRF of Output (dy)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for consumption
nexttile;
plot(1:nPeriods, dc_to_ea_median, 'LineWidth', 2);
title('IRF of consumption (dc)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for Investment
nexttile;
plot(1:nPeriods, dinve_to_ea_median, 'LineWidth', 2);
title('IRF of Investment (dinve)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for wage
nexttile;
plot(1:nPeriods, dw_to_ea_median, 'LineWidth', 2);
title('IRF of wage (dw)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Inflation rate
nexttile;
plot(1:nPeriods, pinfobs_to_ea_median, 'LineWidth', 2);
title('IRF of Inflation rate (pinfobs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Federal funds rate
nexttile;
plot(1:nPeriods, dc_to_ea_median, 'LineWidth', 2);
title('IRF of federal funds rate (robs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for labor hours
nexttile;
plot(1:nPeriods, dc_to_ea_median, 'LineWidth', 2);
title('IRF of labor hours (labobs)');
xlabel('Periods');
ylabel('Deviation');




% Adjusting the layout for better spacing
sgtitle('Impulse Response Functions to a Technology Shock');








 



%Extracting irfs for all observables at the posterior means, mediand and
%modes
%Investment Specific Technology

%Posterior mean
dy_to_eqs = oo_.PosteriorIRF.dsge.Mean.dy_eqs;
dc_to_eqs = oo_.PosteriorIRF.dsge.Mean.dc_eqs;
dinve_to_eqs = oo_.PosteriorIRF.dsge.Mean.dinve_eqs;
dw_to_eqs = oo_.PosteriorIRF.dsge.Mean.dw_eqs;
pinfobs_to_eqs = oo_.PosteriorIRF.dsge.Mean.pinfobs_eqs;
robs_to_eqs = oo_.PosteriorIRF.dsge.Mean.robs_eqs;
labobs_to_eqs = oo_.PosteriorIRF.dsge.Mean.labobs_eqs;


%Posterior median
dy_to_eqs_median = oo_.PosteriorIRF.dsge.Median.dy_eqs;
dc_to_eqs_median = oo_.PosteriorIRF.dsge.Median.dc_eqs;
dinve_to_eqs_median = oo_.PosteriorIRF.dsge.Median.dinve_eqs;
dw_to_eqs_median = oo_.PosteriorIRF.dsge.Median.dw_eqs;
pinfobs_to_eqs_median = oo_.PosteriorIRF.dsge.Median.pinfobs_eqs;
robs_to_eqs_median = oo_.PosteriorIRF.dsge.Median.robs_eqs;
labobs_to_eqs_median = oo_.PosteriorIRF.dsge.Median.labobs_eqs;




% Number of periods to plot
nPeriods = 40;

% Creating a figure using tiledlayout
figure;
tiledlayout(2, 7);  % 4x2 grid for 7 variables

% Ploting IRF for output (y)
nexttile;
plot(1:nPeriods, dy_to_eqs, 'LineWidth', 2);
title('IRF of Output to investment specific technology shock (dy)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for consumption
nexttile;
plot(1:nPeriods, dc_to_eqs, 'LineWidth', 2);
title('IRF of consumption to investment specific technology shock (dc)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for Investment
nexttile;
plot(1:nPeriods, dinve_to_eqs, 'LineWidth', 2);
title('IRF of Investment to investment specific technology shock (dinve)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for wage
nexttile;
plot(1:nPeriods, dw_to_eqs, 'LineWidth', 2);
title('IRF of wage to investment specific technology shock (dw)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Inflation rate
nexttile;
plot(1:nPeriods, pinfobs_to_eqs, 'LineWidth', 2);
title('IRF of Inflation rate to investment specific technology shock (pinfobs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Federal funds rate
nexttile;
plot(1:nPeriods, dc_to_eqs, 'LineWidth', 2);
title('IRF of federal funds rate to investment specific technology shock (robs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for labor hours
nexttile;
plot(1:nPeriods, dc_to_eqs, 'LineWidth', 2);
title('IRF of labor hours to investment specific technology shock (labobs)');
xlabel('Periods');
ylabel('Deviation');


% Ploting IRF for output (y)
nexttile;
plot(1:nPeriods, dy_to_eqs_median, 'LineWidth', 2);
title('IRF of Output to investment specific technology shock (dy)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for consumption
nexttile;
plot(1:nPeriods, dc_to_eqs_median, 'LineWidth', 2);
title('IRF of consumption to investment specific technology shock (dc)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for Investment
nexttile;
plot(1:nPeriods, dinve_to_eqs_median, 'LineWidth', 2);
title('IRF of Investment to investment specific technology shock (dinve)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for wage
nexttile;
plot(1:nPeriods, dw_to_eqs_median, 'LineWidth', 2);
title('IRF of wage to investment specific technology shock (dw)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Inflation rate
nexttile;
plot(1:nPeriods, pinfobs_to_eqs_median, 'LineWidth', 2);
title('IRF of Inflation rate to investment specific technology shock (pinfobs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Federal funds rate
nexttile;
plot(1:nPeriods, dc_to_eqs_median, 'LineWidth', 2);
title('IRF of federal funds rate to investment specific technology shock (robs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for labor hours
nexttile;
plot(1:nPeriods, dc_to_eqs_median, 'LineWidth', 2);
title('IRF of labor hours to investment specific technology shock (labobs)');
xlabel('Periods');
ylabel('Deviation');




% Adjusting the layout for better spacing
sgtitle('Impulse Response Functions to a Investment Specific Technology Shock');



%Extracting irfs for all observables at the posterior means, mrdians, and
%mode
%Investment Specific Technology

%Posterior mean
dy_to_em = oo_.PosteriorIRF.dsge.Mean.dy_em;
dc_to_em = oo_.PosteriorIRF.dsge.Mean.dc_em;
dinve_to_em = oo_.PosteriorIRF.dsge.Mean.dinve_em;
dw_to_em = oo_.PosteriorIRF.dsge.Mean.dw_em;
pinfobs_to_em = oo_.PosteriorIRF.dsge.Mean.pinfobs_em;
robs_to_em = oo_.PosteriorIRF.dsge.Mean.robs_em;
labobs_to_em = oo_.PosteriorIRF.dsge.Mean.labobs_em;

%Posterior median
dy_to_em_median = oo_.PosteriorIRF.dsge.Median.dy_em;
dc_to_em_median = oo_.PosteriorIRF.dsge.Median.dc_em;
dinve_to_em_median = oo_.PosteriorIRF.dsge.Median.dinve_em;
dw_to_em_median = oo_.PosteriorIRF.dsge.Median.dw_em;
pinfobs_to_em_median = oo_.PosteriorIRF.dsge.Median.pinfobs_em;
robs_to_em_median = oo_.PosteriorIRF.dsge.Median.robs_em;
labobs_to_em_median = oo_.PosteriorIRF.dsge.Median.labobs_em;





% Number of periods to plot
nPeriods = 40;

% Creating a figure using tiledlayout
figure;
tiledlayout(2, 7);  % 4x2 grid for 7 variables

% Ploting IRF for output (y)
nexttile;
plot(1:nPeriods, dy_to_em, 'LineWidth', 2);
title('IRF of Output to monetary policy shock (dy)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for consumption
nexttile;
plot(1:nPeriods, dc_to_em, 'LineWidth', 2);
title('IRF of consumption to monetary policy shock (dc)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for Investment
nexttile;
plot(1:nPeriods, dinve_to_em, 'LineWidth', 2);
title('IRF of Investment to monetary policy shock (dinve)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for wage
nexttile;
plot(1:nPeriods, dw_to_em, 'LineWidth', 2);
title('IRF of wage to monetary policy shock (dw)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Inflation rate
nexttile;
plot(1:nPeriods, pinfobs_to_em, 'LineWidth', 2);
title('IRF of Inflation rate to monetary policy shock (pinfobs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Federal funds rate
nexttile;
plot(1:nPeriods, dc_to_em, 'LineWidth', 2);
title('IRF of federal funds rate to monetary policy shock (robs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for labor hours
nexttile;
plot(1:nPeriods, dc_to_em, 'LineWidth', 2);
title('IRF of labor hours to monetary policy shock (labobs)');
xlabel('Periods');
ylabel('Deviation');



% Ploting IRF for output (y)
nexttile;
plot(1:nPeriods, dy_to_em_median, 'LineWidth', 2);
title('IRF of Output to monetary policy shock (dy)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for consumption
nexttile;
plot(1:nPeriods, dc_to_em_median, 'LineWidth', 2);
title('IRF of consumption to monetary policy shock (dc)');
xlabel('Periods');
ylabel('Deviation');

%Plotting IRF for Investment
nexttile;
plot(1:nPeriods, dinve_to_em_median, 'LineWidth', 2);
title('IRF of Investment to monetary policy shock (dinve)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for wage
nexttile;
plot(1:nPeriods, dw_to_em_median, 'LineWidth', 2);
title('IRF of wage to monetary policy shock (dw)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Inflation rate
nexttile;
plot(1:nPeriods, pinfobs_to_em_median, 'LineWidth', 2);
title('IRF of Inflation rate to monetary policy shock (pinfobs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for Federal funds rate
nexttile;
plot(1:nPeriods, dc_to_em_median, 'LineWidth', 2);
title('IRF of federal funds rate to monetary policy shock (robs)');
xlabel('Periods');
ylabel('Deviation');


%Plotting IRF for labor hours
nexttile;
plot(1:nPeriods, dc_to_em_median, 'LineWidth', 2);
title('IRF of labor hours to monetary policy shock (labobs)');
xlabel('Periods');
ylabel('Deviation');




% Adjusting the layout for better spacing
sgtitle('Impulse Response Functions to a monetary policy Shock');




