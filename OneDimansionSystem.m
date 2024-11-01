%% Phase portrait of the system (S=0). 
% Setting the initial values, n has to be >> 1, s = 0 and b = 2 
% according with the instructions of the assignment. 
n = 30;
s = 0;
b = 2; 

% Setting the nondimensionalized equation of the system. 
dadt = @(x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;
figure('name', 'Nondimensionalized system'); % Naming and plotting the equation.
plot1 = fplot(dadt, [-0.5, 3]);
hold on
yline(0)
grid on;

% Locating, setting and plotting the fixed points of the system, 
% according to the plot of the equation we can determine the stability of the fixed points. 
fp1 = [0, 0];
fp2 = [0, 1]; 
fp3 = [0, 2]; 
plot2 = plot(fp1(1), fp1(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot3 = plot(fp2(2), fp2(1), 'ro', 'MarkerSize', 10, 'LineWidth',2);
plot4 = plot(fp3(2), fp3(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x'); ylabel('da/dt');
title('Nondimensionalized system with S = 0');

% Here we are creating the vector field for the phase portrait with an
% anonymous function including the time.
dadtA = @(t,x) dadt(x);
dadtA = @(t,x) dadt(x);
xRange = -0.80:.2:3;
plot5 = quiver(xRange, zeros(size(xRange)), dadt(xRange), zeros(size(xRange)), 1);
legend([plot1, plot2, plot3, plot5], {'s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x', 'Stable Point', 'Unstable Point', 'Vector Field'}, 'Location', 'best');
hold off;

%% Anonymous function with time for the vector slope. 
functionVectorSlope = @(t,x) dadt(x);

%Set the values of the fuction over the time in a grid. 
tStep = 1;
tRange = 0:tStep:10;
xRange = 0:.1:3;
[tValues, xValues] = meshgrid(tRange, xRange);
dx = functionVectorSlope(tValues, xValues);
dt = ones(size(dx)) * tStep;

% The quiver plot here will represent the vector solpe of the fuction
% according with the fixed points.
quiver(tValues, xValues, dt, dx);
axis tight;
xlabel('t'); ylabel('x(t)');
hold on;
title('Nondimensionalized slope field S = 0.');

% plot solution curves (ode45 to solve diff. equation)
xStart = 0; xStep = 0.25; xEnd = 1.5;
for x0 = xStart:xStep:xEnd 
    [tSolution, xSolution] = ode45(functionVectorSlope, [tRange(1), tRange(end)], x0); 
  plot(tSolution, xSolution, 'lineWidth', 2);
end
hold off;

%%  Different phase portraits of the system with a constant stimulus (S > 0) and Bifurcation diagram.
% Setting the initial values, n has to be >> 1, s = 1 (S > 0) and b = 2 
% according with the instructions of the assignment. 
n = 30;
s = 1;
b = 2; 

% Setting the nondimensionalized equation of the system. 
dadt = @(x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;
figure('name', 'Nondimensionalized system'); % Naming and plotting the equation.
plot1 = fplot(dadt, [-0.5, 3]);
hold on
yline(0)
grid on;

% Locating, setting and plotting the fixed point of the system, 
% according to the plot of the equation we can determine the stability of the fixed point. 
fp3 = [0, 2]; 
plot4 = plot(fp3(2), fp3(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x'); ylabel('da/dt');

title('Nondimensionalized system with S > 0 (S = 1)');

% Here we are creating the vector field for the phase portrait with an
% anonymous function including the time.
dadtA = @(t,x) dadt(x);
dadtA = @(t,x) dadt(x);
xRange = -0.5:0.2:3;
plot5 = quiver(xRange, zeros(size(xRange)), dadt(xRange), zeros(size(xRange)), 1);
legend([plot1, plot4, plot5], {'s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x', 'Stable Point', 'Vector Field'}, 'Location', 'best');
hold off;

%% Anonymous function with time for the vector slope. 
functionVectorSlope = @(t,x) dadt(x);

%Set the values of the fuction over the time in a grid. 
tStep = 1;
tRange = 0:tStep:10;
xRange = 0:.1:3;
[tValues, xValues] = meshgrid(tRange, xRange);
dx = functionVectorSlope(tValues, xValues);
dt = ones(size(dx)) * tStep;

% The quiver plot here will represent the vector solpe of the fuction
% according with the fixed points.
quiver(tValues, xValues, dt, dx);
axis tight;
xlabel('t'); ylabel('x(t)');
hold on;
title('Nondimensionalized slope field S = 1.');

% plot solution curves (ode45 to solve diff. equation)
xStart = 0; xStep = 0.25; xEnd = 1.5;
for x0 = xStart:xStep:xEnd 
    [tSolution, xSolution] = ode45(functionVectorSlope, [tRange(1), tRange(end)], x0); 
  plot(tSolution, xSolution, 'lineWidth', 2);
end
hold off;

%% Biurcation diagram.
n = 30;
s = 1;
b = 2; 
dadt = @(s, x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x; 
figure('name', 'Saddle node bifurcation');
% Plot the ODE function for a value of S. 
fplot(@(x) dadt(s, x), [-0.5, 3]);
xlabel('x'); ylabel('da/dt');
xlim([-0.5, 3]); ylim([-2, 2]);
line([-0.5, 3], [0, 0], 'linewidth', 2, 'lineStyle', ':');
grid on;
title('s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x');

%% Plot the location of the fixed point and all the different values for S.
fimplicit(dadt, [-5, 5]);
xlabel('S'); ylabel('x^*')
grid on;
title('Bifurcation diagram');
hold on;

%Plot the direction of the vector field indicating the stability of the 
% fixed points of the system. 
sStep = 0.5;
sRange = -5:sStep:5;
xRange = -5:.5:5;
[sValues, xValues] = meshgrid(sRange, xRange);
dx = dadt(sValues, xValues);
dr = zeros(size(dx));
scaling = sqrt(dx.^2 + dr.^2);
quiver(sValues, xValues, dr./scaling, dx./scaling);
hold off; % It is a Saddle node bifurcation. 

%% Bifurcation diagram of the system increasing stimulus (S).
% Stting the parameters, here s has a range from 1 to 2.5.
n = 30; 
s = 1;
b = 2;
dadt = @(s, x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;

% Plot the bifurcation diagram according with the anonymous function and
% the vector field. 
fimplicit(dadt, [0, 10]);
xlabel('S'); ylabel('x^*')
grid on;
axis([0, 2.5, 0, 3]);
title('Bifurcation diagram in dadt'); 
hold on;
sStep = 0.1;
sRange = -0.75:sStep:2.85;
xRange = -0.75:0.1:2.85;
[sValues, xValues] = meshgrid(sRange, xRange);
dx = dadt(sValues, xValues);
ds = zeros(size(dx));
scaling = sqrt(dx.^2 + ds.^2) / 1;
quiver(sValues, xValues, ds./scaling, dx./scaling);
hold off;

%% Stting the parameters, here s has a range from 1 to 2.5.
n = 10; 
s = 0.515; 
b = 1.2;
dadt = @(s, x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;

% Plot the bifurcation diagram according with the anonymous function and
% the vector field. 
fimplicit(dadt, [0, 3]);
xlabel('S'); ylabel('x^*')
grid on;

title('Bifurcation diagram in dadt'); 

%% Proving that the system is irreversible. 
% Plot the location of the fixed point and all the different values for S.
subplot(1, 2, 1);
fimplicit(dadt, [-5, 5]);
xlabel('S'); ylabel('x^*')
grid on;
title('Bifurcation diagram');
hold on;

%Plot the direction of the vector field indicating the stability of the 
% fixed points of the system. 
sStep = 0.5;
sRange = -5:sStep:5;
xRange = -5:.5:5;
[sValues, xValues] = meshgrid(sRange, xRange);
dx = dadt(sValues, xValues);
dr = zeros(size(dx));
scaling = sqrt(dx.^2 + dr.^2);
quiver(sValues, xValues, dr./scaling, dx./scaling);
hold off; % It is a Saddle node bifurcation.
% Stting the parameters, here s has a range from 1 to 2.5.
n = 30; 
s = 1;
b = 2;
dadt = @(s, x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;

% Plot the bifurcation diagram according with the anonymous function and
% the vector field. 
subplot(1, 2, 2);
fimplicit(dadt, [0, 10]);
xlabel('S'); ylabel('x^*')
grid on;
axis([0, 2.5, 0, 3]);
title('Bifurcation diagram in dadt'); 
hold on;
sStep = 0.1;
sRange = -0.75:sStep:2.85;
xRange = -0.75:0.1:2.85;
[sValues, xValues] = meshgrid(sRange, xRange);
dx = dadt(sValues, xValues);
ds = zeros(size(dx));
scaling = sqrt(dx.^2 + ds.^2) / 1;
quiver(sValues, xValues, ds./scaling, dx./scaling);
hold off;

%% Repeat the analysis above, but this time for a value of β<<kdAT.
%% Phase portrait of the system (S=0). 
% Setting the initial values, n has to be >> 1, s = 0 and b = 0 
% according with the instructions of the assignment and to the new 
% nondimensionalized equation the lower value of beta has to be equal to zero. 
n = 30;
s = 0;
b = 0; 

% Setting the nondimensionalized equation of the system. 
dadt = @(x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;
figure('name', 'Nondimensionalized system'); % Naming and plotting the equation.
plot1 = fplot(dadt, [-0.5, 3]);
hold on
yline(0)
grid on;

% Locating, setting and plotting the fixed points of the system, 
% according to the plot of the equation we can determine the stability of the fixed points. 
fp1 = [0, 0];
plot2 = plot(fp1(1), fp1(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x'); ylabel('da/dt');
title('Nondimensionalized system with S = 0 and b = 0');

% Here we are creating the vector field for the phase portrait with an
% anonymous function including the time.
dadtA = @(t,x) dadt(x);
dadtA = @(t,x) dadt(x);
xRange = -0.80:.2:3;
plot5 = quiver(xRange, zeros(size(xRange)), dadt(xRange), zeros(size(xRange)), 1);
legend([plot1, plot2,plot5], {'s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x', 'Stable Point', 'Vector Field'}, 'Location', 'best');
hold off;

%% Anonymous function with time for the vector slope. 
functionVectorSlope = @(t,x) dadt(x);

%Set the values of the fuction over the time in a grid. 
tStep = 1;
tRange = 0:tStep:10;
xRange = 0:.1:3;
[tValues, xValues] = meshgrid(tRange, xRange);
dx = functionVectorSlope(tValues, xValues);
dt = ones(size(dx)) * tStep;

% The quiver plot here will represent the vector solpe of the fuction
% according with the fixed points.
quiver(tValues, xValues, dt, dx);
axis tight;
xlabel('t'); ylabel('x(t)');
hold on;
title('Nondimensionalized slope field with S = 0 and b = 0');

% plot solution curves (ode45 to solve diff. equation)
xStart = 0; xStep = 0.25; xEnd = 1.5;
for x0 = xStart:xStep:xEnd 
    [tSolution, xSolution] = ode45(functionVectorSlope, [tRange(1), tRange(end)], x0); 
  plot(tSolution, xSolution, 'lineWidth', 2);
end
hold off;

%% Different phase portraits of the system with a constant stimulus (S > 0) 
% and Bifurcation diagram.
% Setting the initial values, n has to be >> 1, s = 1 (S > 0) and b = 0 
% according with the instructions of the assignment and and to the new 
% nondimensionalized equation the lower value of beta has to be equal to zero. 
n = 30;
s = 1;
b = 0; 

% Setting the nondimensionalized equation of the system. 
dadt = @(x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;
figure('name', 'Nondimensionalized system'); % Naming and plotting the equation.
plot1 = fplot(dadt, [-0.5, 3]);
hold on
yline(0)
grid on;

% Locating, setting and plotting the fixed point of the system, 
% according to the plot of the equation we can determine the stability of the fixed point. 
fp3 = [1, 0]; 
plot4 = plot(fp3(1), fp3(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x'); ylabel('da/dt');

title('Nondimensionalized system with S > 0 (S = 1) and b = 0');

% Here we are creating the vector field for the phase portrait with an
% anonymous function including the time.
dadtA = @(t,x) dadt(x);
dadtA = @(t,x) dadt(x);
xRange = -0.5:0.2:3;
plot5 = quiver(xRange, zeros(size(xRange)), dadt(xRange), zeros(size(xRange)), 1);
legend([plot1, plot4, plot5], {'s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x', 'Stable Point', 'Vector Field'}, 'Location', 'best');
hold off;

%% Anonymous function with time for the vector slope. 
functionVectorSlope = @(t,x) dadt(x);

%Set the values of the fuction over the time in a grid. 
tStep = 1;
tRange = 0:tStep:10;
xRange = 0:.1:3;
[tValues, xValues] = meshgrid(tRange, xRange);
dx = functionVectorSlope(tValues, xValues);
dt = ones(size(dx)) * tStep;

% The quiver plot here will represent the vector solpe of the fuction
% according with the fixed points.
quiver(tValues, xValues, dt, dx);
axis tight;
xlabel('t'); ylabel('x(t)');
hold on;
title('Nondimensionalized slope field with S = 1 and b = 0');

% plot solution curves (ode45 to solve diff. equation)
xStart = 0; xStep = 0.25; xEnd = 1.5;
for x0 = xStart:xStep:xEnd 
    [tSolution, xSolution] = ode45(functionVectorSlope, [tRange(1), tRange(end)], x0); 
  plot(tSolution, xSolution, 'lineWidth', 2);
end
hold off;

%% Bifurcation diagram
n = 30;
s = 1;
b = 0; 
dadt = @(s, x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x; 
figure('name', 'Saddle node bifurcation');
% Plot the ODE function for a value of S. 
fplot(@(x) dadt(s, x), [-0.5, 3]);
xlabel('x'); ylabel('da/dt');
xlim([-0.5, 3]); ylim([-2, 2]);
line([-0.5, 3], [0, 0], 'linewidth', 2, 'lineStyle', ':');
grid on;
title('s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x');

%%% Plot the location of the fixed point and all the different values for S.
fimplicit(dadt, [-5, 5]);
xlabel('S'); ylabel('x^*')
grid on;
title('Bifurcation diagram b = 0');
hold on;

%Plot the direction of the vector field indicating the stability of the 
% fixed points of the system. 
sStep = 0.5;
sRange = -5:sStep:5;
xRange = -5:.5:5;
[sValues, xValues] = meshgrid(sRange, xRange);
dx = dadt(sValues, xValues);
dr = zeros(size(dx));
scaling = sqrt(dx.^2 + dr.^2);
quiver(sValues, xValues, dr./scaling, dx./scaling);
hold off; % It is a Saddle node bifurcation.

%% Fixed points in the bifurcation diagram. 
% Stting the parameters, here s has a range from 1 to 2.5.
n = 30; 
s = 1;
b = 0;
dadt = @(s, x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;

% Plot the bifurcation diagram according with the anonymous function and
% the vector field. 
fimplicit(dadt, [0, 10]);
xlabel('S'); ylabel('x^*')
grid on;
axis([0, 2.5, 0, 3]);
title('Bifurcation diagram in dadt b = 0'); 
hold on;
sStep = 0.1;
sRange = -0.75:sStep:2.85;
xRange = -0.75:0.1:2.85;
[sValues, xValues] = meshgrid(sRange, xRange);
dx = dadt(sValues, xValues);
ds = zeros(size(dx));
scaling = sqrt(dx.^2 + ds.^2) / 1;
quiver(sValues, xValues, ds./scaling, dx./scaling);
hold off;

%% Proving that the system is reversible.
% Plot the location of the fixed point and all the different values for S.
subplot(1, 2, 1);
fimplicit(dadt, [-5, 5]);
xlabel('S'); ylabel('x^*')
grid on;
title('Bifurcation diagram b = 0');
hold on;

%Plot the direction of the vector field indicating the stability of the 
% fixed points of the system. 
sStep = 0.5;
sRange = -5:sStep:5;
xRange = -5:.5:5;
[sValues, xValues] = meshgrid(sRange, xRange);
dx = dadt(sValues, xValues);
dr = zeros(size(dx));
scaling = sqrt(dx.^2 + dr.^2);
quiver(sValues, xValues, dr./scaling, dx./scaling);
hold off; % It is a Saddle node bifurcation.
% Stting the parameters, here s has a range from 1 to 2.5.
n = 30; 
s = 1;
b = 0;
dadt = @(s, x) s.*(2-x) + b.*((x.^n)./(1+x.^n)) - x;

% Plot the bifurcation diagram according with the anonymous function and
% the vector field. 
subplot(1, 2, 2);
fimplicit(dadt, [0, 10]);
xlabel('S'); ylabel('x^*')
grid on;
axis([0, 2.5, 0, 3]);
title('Bifurcation diagram in dadt b = 0'); 
hold on;
sStep = 0.1;
sRange = -0.75:sStep:2.85;
xRange = -0.75:0.1:2.85;
[sValues, xValues] = meshgrid(sRange, xRange);
dx = dadt(sValues, xValues);
ds = zeros(size(dx));
scaling = sqrt(dx.^2 + ds.^2) / 1;
quiver(sValues, xValues, ds./scaling, dx./scaling);
hold off;