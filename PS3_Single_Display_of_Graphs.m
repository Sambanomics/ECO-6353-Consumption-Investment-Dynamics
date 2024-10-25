
%% Combining the Graphs

% Number of periods to plot
nPeriods = 14;

% Creating a Figure using Tiledlayout
figure;
tiledlayout(7, 2);  


% Generating and ploting the 14 graphs

for i =  1:n_graphs
    
    subplot(n_rows, n_cols, i);


% Plotting for Demmeaned Consumption
nexttile;
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

% Plotting for Demmeaned RGDP
nexttile;
plot(dates, demeaned_inc, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned RGDP');
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


% Plotting for Demmeaned Investment
nexttile;
plot(dates, demeaned_inv, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned investment');
title('Deemeaned investment with NBER recession regions');
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


% Plotting for Demmeaned Inflation Rate
nexttile;
plot(dates, demeaned_cons, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned inflation rate');
title('Deemeaned inflation rate with NBER recession regions');
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


% Plotting for Demmeaned Federal Funds Rate
nexttile;
plot(dates, demeaned_ffrate, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned federal funds rate');
title('Deemeaned federal funds rate with NBER recession regions');
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


% Plotting for Demmeaned Wage
nexttile;
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


% Plotting for Demmeaned Hours-Worked
nexttile;
plot(dates, demeaned_hrwrked, 'LineWidth', 2);
xlabel('Time');
ylabel('demeaned hours worked');
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


% Plotting for Detrended Consumption
nexttile;
plot(dates, cycle, 'LineWidth', 2);
xlabel('Time');
ylabel('consumption');
title('Consumption with NBER recession regions');
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


% Plotting for Detrended Investment
nexttile;
plot(dates, cycle1, 'LineWidth', 2);
xlabel('Time');
ylabel('Investment');
title('Investment with NBER recession regions');
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


% Plotting for Detrended RGDP
nexttile;
plot(dates, cycle2, 'LineWidth', 2);
xlabel('Time');
ylabel('RGDP');
title('RGDP with NBER recession regions');
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


% Plotting for Detrended Hours-Worked
nexttile;
plot(dates, cycle3, 'LineWidth', 2);
xlabel('Time');
ylabel('Hours worked');
title('Hours worked with NBER recession regions');
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


% Plotting for Detrended Inflation rate
nexttile;
plot(dates, cycle4, 'LineWidth', 2);
xlabel('Time');
ylabel('Inflation rate');
title('Inflation rate with NBER recession regions');
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



% Plotting for Detrended Wage
nexttile;
plot(dates, cycle5, 'LineWidth', 2);
xlabel('Time');
ylabel('wage');
title('Wage with NBER recession regions');
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


% Plotting for Detrended Federal Funds Rate
nexttile;
plot(dates, cycle6, 'LineWidth', 2);
xlabel('Time');
ylabel('federal funds rate');
title('federal funds rate with NBER recession regions');
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

% Adjusting the Layout for better spacing
sgtitle('Graph of Demeaned and Detrended Variables');

end



