

% Loading the Data (Using "Import Data" tool)

Data = Data_PS3;

disp(Data(1:5,:));

% Renaming the columns

dates = Data{:,1};  
cons = Data{:,2}; 
inv = Data{:,3};
inc = Data{:,4};
hrwrked = Data{:,5};
inf = Data{:,6};
wage = Data{:,7};
ffrate = Data{:,8};
recession = Data{:,9};
hrwrked1 = Data{:,10};
inflation_cpi = Data{:,11};

% Using HP Filter Parameters

lambda = 1600;

% Applying HP filter for Consumption

[trend, cycle] = hpfilter(cons);

% Plotting the results to depict the detrended data

figure;
subplot(2,1,1);
plot(dates, cons, 'LineWidth', 2); hold on;
plot(dates, trend, '--r', 'LineWidth', 2);
title('consumption and HP Filtered Trend');
legend('consumption', 'Trend');
hold on;
grid on;

subplot(2,1,2);
plot(dates, cycle, 'LineWidth', 2);
title('HP Filtered Cycle (Detrended Consumption)');
hold on;
grid on;


in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Applying the HP filter for Investment

[trend1, cycle1] = hpfilter(inv);

% Plotting the results to depict the detrended data
figure;
subplot(2,1,1);
plot(dates, inv, 'LineWidth', 2); hold on;
plot(dates, trend1, '--r', 'LineWidth', 2);
title('Investment and HP Filtered Trend');
legend('Investment', 'Trend');
hold on;
grid on;

subplot(2,1,2);
plot(dates, cycle1, 'LineWidth', 2);
title('HP Filtered Cycle (Detrended Investment)');
hold on;
grid on;


in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Applying the HP filter for Income

[trend2, cycle2] = hpfilter(inc);

% Plotting the results to depict the detrended data
figure;
subplot(2,1,1);
plot(dates, inc, 'LineWidth', 2); hold on;
plot(dates, trend2, '--r', 'LineWidth', 2);
title('Income and HP Filtered Trend');
legend('Income', 'Trend');
hold on;
grid on;

subplot(2,1,2);
plot(dates, cycle2, 'LineWidth', 2);
title('HP Filtered Cycle (Detrended Income)');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Applying the HP filter for Hours-Worked

[trend3, cycle3] = hpfilter(hrwrked);

% Plotting the results to depict the detrended data
figure;
subplot(2,1,1);
plot(dates, hrwrked, 'LineWidth', 2); hold on;
plot(dates, trend3, '--r', 'LineWidth', 2);
title('Hours worked and HP Filtered Trend');
legend('Hours worked', 'Trend');
hold on;
grid on;

subplot(2,1,2);
plot(dates, cycle3, 'LineWidth', 2);
title('HP Filtered Cycle (Detrended Hours worked)');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;



% Applying the HP filter for Inflation

[trend4, cycle4] = hpfilter(inf);

% Plotting the results to depict the detrended data
figure;
subplot(2,1,1);
plot(dates, inf, 'LineWidth', 2); hold on;
plot(dates, trend4, '--r', 'LineWidth', 2);
title('Inflation rate and HP Filtered Trend');
legend('Inflation rate', 'Trend');
hold on;
grid on;

subplot(2,1,2);
plot(dates, cycle4, 'LineWidth', 2);
title('HP Filtered Cycle (Detrended Inflation rate)');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Applying the HP filter for Wage

[trend5, cycle5] = hpfilter(wage);

% Plotting the results to depict the detrended data
figure;
subplot(2,1,1);
plot(dates, wage, 'LineWidth', 2); hold on;
plot(dates, trend5, '--r', 'LineWidth', 2);
title('wage and HP Filtered Trend');
legend('Wage', 'Trend');
hold on;
grid on;

subplot(2,1,2);
plot(dates, cycle5, 'LineWidth', 2);
title('HP Filtered Cycle (Detrended Wage)');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;



% Applying the HP filter for Federal Funds Rate

[trend6, cycle6] = hpfilter(ffrate);

% Plotting the results to depict the detrended data
figure;
subplot(2,1,1);
plot(dates, ffrate, 'LineWidth', 2); hold on;
plot(dates, trend6, '--r', 'LineWidth', 2);
title('Federal Funds rate and HP Filtered Trend');
legend('Federal Funds rate', 'Trend');
hold on;
grid on;

subplot(2,1,2);
plot(dates, cycle6, 'LineWidth', 2);
title('HP Filtered Cycle (Detrended Federal Funds rate)');
hold on;
grid on;


in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;




%% Demeaning the Data

% Calculating the mean of each series

mean_cons = mean(cons);
mean_inc = mean(inc);
mean_inv = mean(inv);
mean_inf = mean(inf);
mean_ffrate = mean(ffrate);
mean_hrwrked = mean(hrwrked);
mean_wage = mean(wage);

% Demeaning the variables

demeaned_cons = cons - mean_cons;
demeaned_inc = inc - mean_inc;
demeaned_inv = inv - mean_inv;
demeaned_inf = inf - mean_inf;
demeaned_ffrate = ffrate - mean_ffrate;
demeaned_hrwrked = hrwrked - mean_hrwrked;
demeaned_wage = wage - mean_wage;

% Plotting the Graphs for the demeaned variables

% Plotting for Consumption
figure;
plot(dates, demeaned_cons, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned consumption');
title('Deemeaned Consumption with NBER recession regions');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Plotting for Income
figure;
plot(dates, demeaned_inc, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned real gross domestic product');
title('Deemeaned RGDP with NBER recession regions');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Plotting for Investment
figure;
plot(dates, demeaned_inv, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned Inveestment');
title('Deemeaned Investment with NBER recession regions');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;

%Plotting demmeaned_cons
figure;
plot(dates, demeaned_inf, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned Inflation rate');
title('Deemeaned Inflation ratte with NBER recession regions');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        %start of a recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        %endd of aa recession
        r_end = dates(i);
        %add shaading for a recession period
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end

%Bring the Inflation rate plot to the front
uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Plotting for Wage
figure;
plot(dates, demeaned_wage, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned wage');
title('Deemeaned wage with NBER recession regions');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Plotting for Federal Funds Rate
figure;
plot(dates, demeaned_ffrate, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned federal funds rate');
title('Deemeaned federal with NBER recession regions');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;


% Plotting for Hours-Worked
figure;
plot(dates, demeaned_hrwrked, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned Hours worked');
title('Deemeaned hours worked with NBER recession regions');
hold on;
grid on;

in_recession = false;
for i = 1:length(recession)
    if recession(i)== 1 && ~in_recession
        r_start = dates(i);
        in_recession = true;
    elseif recession(i)== 0 && in_recession
        r_end = dates(i);
        y_limits = ylim;
        fill([r_start r_start r_end r_end], [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
        in_recession = false;
    end
end


uistack(findobj(gca, 'Type', 'line'), 'top');
hold off;



