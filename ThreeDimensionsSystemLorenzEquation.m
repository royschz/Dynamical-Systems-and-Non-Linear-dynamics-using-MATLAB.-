%% Illustrating the deterministic chaos of the systems, with the parameters a1 = 5, b1 = 3, a2 = 0.1, b2 = 2, d1 = 0.4 and d2 = 0.01. 
% Setting the parameters values. 
a1 = 5;
b1 = 3; 
a2 = 0.1;
b2 = 2; 
d1 = 0.4; 
d2 = 0.01; 

% Setting the functions. 
xHandle = @(x,y,z) x.*(1 - x) - (a1.*x.*y)./(1 + b1.*x); 
yHandle = @(x,y,z) (a1.*x.*y)./(1 + b1.*x) - d1.*y - (a2.*y.*z)./(1 + b2.*y); 
zHandle = @(x,y,z) (a2.*y.*z)./(1 + b2.*y) - d2.*z;

HasPow = @(t, x) [...
    xHandle(x(1), x(2), x(3));...
    yHandle(x(1), x(2), x(3));...
    zHandle(x(1), x(2), x(3))];

% Simulation 
x0 = [2, 2, 2];
[tSolution, xSolution] = ode45(HasPow, [0, 10000], x0);

% Plot the first equation (Plants). 
figure; 
subplot(2, 2, 1);
plot(tSolution, xSolution(:,1)); 
title('Plants');
xlabel('Time');
ylabel('Plants');

% Plot the second equation (Herbivores). 
subplot(2, 2, 2);
plot(tSolution, xSolution(:,2)); 
title('Herbivores');
xlabel('Time');
ylabel('Herbivores');

% Plot the third equation (Carnivores). 
subplot(2, 2, 3);
plot(tSolution, xSolution(:,3)); 
title('Carnivores');
xlabel('Time');
ylabel('Carnivores');

% Plot the three populations together. 
subplot(2, 2, 4);
plot(tSolution, xSolution(:, 1:3));
title('Hastings and Powell Model');
xlabel('Time');
ylabel('Populations');
legend({'Plants', 'Herbivores', 'Carnivores'});

%% Plot the dynamics of the system. 
figure; 
plot3(xSolution(:,1), xSolution(:,2), xSolution(:,3));
title('Hastings and Powell Model');
xlabel('Plants'); 
ylabel('Herbivores'); 
zlabel('Carnivores'); 
rotate3d on; 

%% New simulation starting close to the strange attractor. 
% Setting the parameters values. 
a1 = 5;
b1 = 3; 
a2 = 0.1;
b2 = 2; 
d1 = 0.4; 
d2 = 0.01; 

% Setting the functions. 
xHandle = @(x,y,z) x.*(1 - x) - (a1.*x.*y)./(1 + b1.*x); 
yHandle = @(x,y,z) (a1.*x.*y)./(1 + b1.*x) - d1.*y - (a2.*y.*z)./(1 + b2.*y); 
zHandle = @(x,y,z) (a2.*y.*z)./(1 + b2.*y) - d2.*z;

HasPow = @(t, x) [...
    xHandle(x(1), x(2), x(3));...
    yHandle(x(1), x(2), x(3));...
    zHandle(x(1), x(2), x(3))];

% New simulation, taking the last points if the xSolutions. 
%x0New = [0.933212439725246, 0.045097418644185, 10.078896369518654];
x0New = xSolution(end,:);
[tSolutionNew, xSolutionNew] = ode45(HasPow, [0, 10000], x0New);

% Plot the first equation (Plants). 
figure;
subplot(2, 2, 1);
plot(tSolutionNew, xSolutionNew(:,1)); 
title('Plants');
xlabel('Time');
ylabel('Plants');

% Plot the second equation (Herbivores). 
subplot(2, 2, 2);
plot(tSolutionNew, xSolutionNew(:,2)); 
title('Herbivores');
xlabel('Time');
ylabel('Herbivores');

% Plot the third equation (Carnivores). 
subplot(2, 2, 3);
plot(tSolutionNew, xSolutionNew(:,3)); 
title('Carnivores');
xlabel('Time');
ylabel('Carnivores');

% Plot the three populations together. 
subplot(2, 2, 4);
plot(tSolutionNew, xSolutionNew(:, 1:3));
title('Hastings and Powell Model');
xlabel('Time');
ylabel('Populations');
legend({'Plants', 'Herbivores', 'Carnivores'});

%% Plot the solutions of the system. 
figure;
plot3(xSolutionNew(:,1), xSolutionNew(:,2), xSolutionNew(:,3));
title('Hastings and Powell Model');
xlabel('Plants'); 
ylabel('Herbivores'); 
zlabel('Carnivores'); 
rotate3d on;

%% Liapunov exponent of the model. 
% I. Simulate the trajectory until the strange attractor and taking the end
% of the trajectory as the new initial condition. 
tEnd = 5000; 
tStep = 0.01; 
[t, trajectory] = ode45(HasPow, 0:tStep:tEnd, x0); 
attractorXO = trajectory(end, :); 

% II. Simulate the trajectory with a fixed time step and generate a second initial 
% condition. 
icChange = 1e-15; 
X0close = attractorXO + icChange; 
[t, trajectory] = ode45(HasPow, 0:tStep:tEnd, attractorXO); 
[t, trajectoryClose] = ode45(HasPow, 0:tStep:tEnd, X0close); 

% III. Compute the distance between the two trajectories. 
distance = sqrt(sum((trajectory - trajectoryClose).^2, 2)); 

% IV. Plot the natural logarithm of the distance (ln) as a function of
% time.
figure; 
subplot(2, 1, 1); 
plot(t, trajectory(:,1), 'DisplayName', 'I.C.'); 
hold on; 
legend;
plot(t, trajectoryClose(:,1), 'DisplayName', 'close');
hold off
subplot(2, 1, 2); 
plot(t, log(distance), 'DisplayName', 'distance');
title('distance'); 
xlabel('t'); 
ylabel('ln(distance)'); 

% V. Measure the initial slope of this curve, straight line fit. 
validTime = t > 1 & t < 2000;
mdl = fitlm(t(validTime), log(distance(validTime)))

%% Vary the value of the parameter b1 in the interval 2 < b1 < 3. 
close all 
clear all
% b1 = 2, Stable fixed point. 
% Setting the parameters values. 
a1 = 5;
b1 = 2;  
a2 = 0.1;
b2 = 2; 
d1 = 0.4; 
d2 = 0.01; 

% Setting the functions. 
xHandle = @(x,y,z) x.*(1 - x) - (a1.*x.*y)./(1 + b1.*x); 
yHandle = @(x,y,z) (a1.*x.*y)./(1 + b1.*x) - d1.*y - (a2.*y.*z)./(1 + b2.*y); 
zHandle = @(x,y,z) (a2.*y.*z)./(1 + b2.*y) - d2.*z;

HasPow = @(t, x) [...
    xHandle(x(1), x(2), x(3));...
    yHandle(x(1), x(2), x(3));...
    zHandle(x(1), x(2), x(3))];

% Simulation 
x0 = [0.4; 0.2; 9];
[tSolution, xSolution] = ode45(HasPow, [0, 10000], x0);

% Plot the first equation (Plants). 
figure; 
subplot(1, 3, 1);
plot(tSolution, xSolution(:,1)); 
title('Plants');
xlabel('Time');
ylabel('Plants');

% Plot the second equation (Herbivores). 
subplot(1, 3, 2);
plot(tSolution, xSolution(:,2)); 
title('Herbivores');
xlabel('Time');
ylabel('Herbivores');

% Plot the third equation (Carnivores). 
subplot(1, 3, 3);
plot(tSolution, xSolution(:,3)); 
title('Carnivores');
xlabel('Time');
ylabel('Carnivores');

%% Plot the three populations together. 
figure
subplot(1, 2, 1);
plot(tSolution, xSolution(:, 1:3));
title('Hastings and Powell Model');
xlabel('Time');
ylabel('Populations');
legend({'Plants', 'Herbivores', 'Carnivores'});
% Plot the dynamics of the system. 
subplot(1, 2, 2) 
plot3(xSolution(:,1), xSolution(:,2), xSolution(:,3));
title('Hastings and Powell Model');
xlabel('Plants'); 
ylabel('Herbivores'); 
zlabel('Carnivores'); 
rotate3d on;

%% 
close all 
clear all
% b1 = 2.15, Hopf bifurcation. 
% Setting the parameters values. 
a1 = 5;
b1 = 2.15;  
a2 = 0.1;
b2 = 2; 
d1 = 0.4; 
d2 = 0.01; 

% Setting the functions. 
xHandle = @(x,y,z) x.*(1 - x) - (a1.*x.*y)./(1 + b1.*x); 
yHandle = @(x,y,z) (a1.*x.*y)./(1 + b1.*x) - d1.*y - (a2.*y.*z)./(1 + b2.*y); 
zHandle = @(x,y,z) (a2.*y.*z)./(1 + b2.*y) - d2.*z;

HasPow = @(t, x) [...
    xHandle(x(1), x(2), x(3));...
    yHandle(x(1), x(2), x(3));...
    zHandle(x(1), x(2), x(3))];

% Simulation 
x0 = [0.4; 0.2; 9];
[tSolution, xSolution] = ode45(HasPow, [0, 10000], x0);

% Plot the first equation (Plants). 
figure; 
subplot(1, 3, 1);
plot(tSolution, xSolution(:,1)); 
title('Plants');
xlabel('Time');
ylabel('Plants');

% Plot the second equation (Herbivores). 
subplot(1, 3, 2);
plot(tSolution, xSolution(:,2)); 
title('Herbivores');
xlabel('Time');
ylabel('Herbivores');

% Plot the third equation (Carnivores). 
subplot(1, 3, 3);
plot(tSolution, xSolution(:,3)); 
title('Carnivores');
xlabel('Time');
ylabel('Carnivores');

%% Plot the three populations together. 
figure;
subplot(1, 2, 1);
plot(tSolution, xSolution(:, 1:3));
title('Hastings and Powell Model');
xlabel('Time');
ylabel('Populations');
legend({'Plants', 'Herbivores', 'Carnivores'});
% Plot the dynamics of the system. 
subplot(1, 2, 2); 
plot3(xSolution(:,1), xSolution(:,2), xSolution(:,3));
title('Hastings and Powell Model');
xlabel('Plants'); 
ylabel('Herbivores'); 
zlabel('Carnivores'); 
rotate3d on; 

%% 
close all 
clear all
% b1 = 2.369, period-2 rythm. 
% Setting the parameters values. 
a1 = 5;
b1 = 2.369;  
a2 = 0.1;
b2 = 2; 
d1 = 0.4; 
d2 = 0.01; 

% Setting the functions. 
xHandle = @(x,y,z) x.*(1 - x) - (a1.*x.*y)./(1 + b1.*x); 
yHandle = @(x,y,z) (a1.*x.*y)./(1 + b1.*x) - d1.*y - (a2.*y.*z)./(1 + b2.*y); 
zHandle = @(x,y,z) (a2.*y.*z)./(1 + b2.*y) - d2.*z;

HasPow = @(t, x) [...
    xHandle(x(1), x(2), x(3));...
    yHandle(x(1), x(2), x(3));...
    zHandle(x(1), x(2), x(3))];

% Simulation 
x0 = [0.4; 0.2; 9];
[tSolution, xSolution] = ode45(HasPow, [0, 10000], x0);

% Plot the first equation (Plants). 
figure; 
subplot(1, 3, 1);
plot(tSolution, xSolution(:,1)); 
title('Plants');
xlabel('Time');
ylabel('Plants');

% Plot the second equation (Herbivores). 
subplot(1, 3, 2);
plot(tSolution, xSolution(:,2)); 
title('Herbivores');
xlabel('Time');
ylabel('Herbivores');

% Plot the third equation (Carnivores). 
subplot(1, 3, 3);
plot(tSolution, xSolution(:,3)); 
title('Carnivores');
xlabel('Time');
ylabel('Carnivores');

%% Plot the three populations together. 
figure; 
subplot(1, 2, 1);
plot(tSolution, xSolution(:, 1:3));
title('Hastings and Powell Model');
xlabel('Time');
ylabel('Populations');
legend({'Plants', 'Herbivores', 'Carnivores'});
% Plot the dynamics of the system. 
subplot(1, 2, 2);
plot3(xSolution(:,1), xSolution(:,2), xSolution(:,3));
title('Hastings and Powell Model');
xlabel('Plants'); 
ylabel('Herbivores'); 
zlabel('Carnivores'); 
rotate3d on; 

%% 
close all 
clear all
% b1 = 2.38, period-4 rythm. 
% Setting the parameters values. 
a1 = 5;
b1 = 2.38;  
a2 = 0.1;
b2 = 2; 
d1 = 0.4; 
d2 = 0.01; 

% Setting the functions. 
xHandle = @(x,y,z) x.*(1 - x) - (a1.*x.*y)./(1 + b1.*x); 
yHandle = @(x,y,z) (a1.*x.*y)./(1 + b1.*x) - d1.*y - (a2.*y.*z)./(1 + b2.*y); 
zHandle = @(x,y,z) (a2.*y.*z)./(1 + b2.*y) - d2.*z;

HasPow = @(t, x) [...
    xHandle(x(1), x(2), x(3));...
    yHandle(x(1), x(2), x(3));...
    zHandle(x(1), x(2), x(3))];

% Simulation 
x0 = [0.4; 0.2; 9];
[tSolution, xSolution] = ode45(HasPow, [0, 10000], x0);

% Plot the first equation (Plants). 
figure; 
subplot(1, 3, 1);
plot(tSolution, xSolution(:,1)); 
title('Plants');
xlabel('Time');
ylabel('Plants');

% Plot the second equation (Herbivores). 
subplot(1, 3, 2);
plot(tSolution, xSolution(:,2)); 
title('Herbivores');
xlabel('Time');
ylabel('Herbivores');

% Plot the third equation (Carnivores). 
subplot(1, 3, 3);
plot(tSolution, xSolution(:,3)); 
title('Carnivores');
xlabel('Time');
ylabel('Carnivores');

%% Plot the three populations together. 
figure; 
subplot(1 ,2, 1);
plot(tSolution, xSolution(:, 1:3));
title('Hastings and Powell Model');
xlabel('Time');
ylabel('Populations');
legend({'Plants', 'Herbivores', 'Carnivores'});
% Plot the dynamics of the system. 
subplot(1, 2, 2);
plot3(xSolution(:,1), xSolution(:,2), xSolution(:,3));
title('Hastings and Powell Model');
xlabel('Plants'); 
ylabel('Herbivores'); 
zlabel('Carnivores'); 
rotate3d on; 