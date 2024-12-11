%***** Temperature gradient of Borehole D *******************************

% input data, D = depth, T = temperature
D = [51,150,249,351,450,550,649,749,850,950,1049,1150,1253,1352,1449,1552,1651]/1000;
T = [6.5,11.4,15.7,21.1,26,30.3,33.9,38.4,43,47.5,52.6,56.1,60.6,64.8,69.8,74.8,79];

% plot data
figure();
plot(D, T, 'bo');
xlabel('Depth [km]', 'FontSize', 15, 'FontName','Times New Roman');
ylabel('Temperature [\circC]', 'FontSize', 15, 'FontName','Times New Roman');
title('Borehole D Temperature With Depth', 'FontSize', 15,'FontName','Times New Roman');
hold on;

% find slope of line
p = polyfit(D,T,1);
xFit = linspace(min(D), max(D), 1000);
yFit = polyval(p, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 1);
annotation('textbox', [.15 .8 .1 .1], 'String', ['Gradient = ', ...
    num2str(round(p(1,1))), '\circC/km'], 'FontSize', 12, 'FontName','Times New Roman');

