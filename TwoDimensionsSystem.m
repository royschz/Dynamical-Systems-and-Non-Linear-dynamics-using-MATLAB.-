%% Phase portrait of the system with the parameters r1 = 1, r2 = 0.1, K = 7, d = 1, j = 1 and w = 0.3.
% Set the parameter values for the system. 
r1 = 1;
r2 = 0.1; 
K = 7; 
d = 1;
j = 1; 
w = 0.3; 

% Define the differential equations. 
dNdtprey = @(N, P) r1*N.*(1 - (N./K)) - (w*P.*N./(d + N));
dPdtpredator = @(N, P) r2*P.*(1 - (j*P./N));

% Calulate the nullclines.
figure
int = [-0.01,10];
fimplicit(dNdtprey, int, 'b', 'displayname', 'dNdtprey = 0');
hold on;
fimplicit(dPdtpredator, int, 'r', 'displayname', 'dPdtpredator = 0');
legend;
[N, P] = meshgrid(-1:.5:10, -1:.5:10);
DN = dNdtprey(N, P);
DP = dPdtpredator(N, P);

% Setting the trajectories. 
scaling = sqrt(DN.^2 + DP.^2);
quiver(N, P, DN./scaling, DP./scaling, 0.3, 'k', 'DisplayName', 'Trajectories');
functionh = @(t,x) [dNdtprey(x(1), x(2)); dPdtpredator(x(1), x(2))];
N0 = 0.1;
P0 = 0.1;
[~, trajectory] = ode45(functionh, [0,100], [N0;P0]);
plot(trajectory(:, 1), trajectory(:, 2), 'LineWidth', 2, 'DisplayName', 'trajectory');

% Fiexed Points 
plot(7,0, "ko", 'MarkerSize', 7, 'DisplayName','Saddle point');
plot(5.2367, 5.2367, "k.", "MarkerSize", 20, 'DisplayName','Stable point');
title('Portrait fase');
xlabel('N');
ylabel('P');
xlim([-1, 10]);
ylim([-1, 10]);

%% Fixed points of the system and their Jacobian Matrix. 
% Jacobian matrix. 
% Parameters 
r1 = 1;
r2 = 0.1; 
K = 7; 
d = 1;
j = 1; 
w = 0.3; 
%% For the fixed point N = 7 and  P = 0. 
 N = 7; 
 P = 0; 

JacoMax1 = [r1-d*w*(P/((d+N)^2))-2*r1*(N/K), -(N*w)/(d+N);
    (j*(P^2)*r2)/(N^2), r2*(1-(2*j*P)/N)];  
[EigVec1 EgiVal1] = eig(JacoMax1)

% Saddle node. 

%% For the fixed point N = 39/20 + (sqrt(4321)/20) and P = 39/20 + (sqrt(4321)/20);
N = (1/20)*(39 + sqrt(4321)); 
P = (1/20)*(39 + sqrt(4321));
JacoMax2 = [r1-d*w*(P/((d+N)^2))-2*r1*(N/K), -(N*w)/(d+N);
    (j*(P^2)*r2)/(N^2), r2*(1-(2*j*P)/N)];  
[EigVec2 EgiVal2] = eig(JacoMax2)

% Stable node.
%% For the fixed point N = 7 and  P = 0. 
 N = 0; 
 P = 0; 

JacoMax3 = [r1-d*w*(P/((d+N)^2))-2*r1*(N/K), -(N*w)/(d+N);
    (j*(P^2)*r2)/(N^2), r2*(1-(2*j*P)/N)];  
[EigVec3 EgiVal3] = eig(JacoMax3)

%% Phase portrait with w = 1. 
% Set the parameter values for the system. 
r1 = 1;
r2 = 0.1; 
K = 7; 
d = 1;
j = 1; 
w = 1; 
% Define the differential equations. 
dNdtprey = @(N, P) r1*N.*(1 - (N./K)) - (w*P.*N./(d + N));
dPdtpredator = @(N, P) r2*P.*(1 - (j*P./N));

% Calulate the nullclines.
figure
int = [-0.01,10];
fimplicit(dNdtprey, int, 'b', 'displayname', 'dNdtprey = 0');
hold on;
fimplicit(dPdtpredator, int, 'r', 'displayname', 'dPdtpredator = 0');
legend;
[N, P] = meshgrid(-1:.5:10, -1:.5:10);
DN = dNdtprey(N, P);
DP = dPdtpredator(N, P);

% Setting the trajectories.
scaling = sqrt(DN.^2 + DP.^2);
quiver(N, P, DN./scaling, DP./scaling, 0.3, 'k', 'DisplayName', 'Trajectories');
functionh = @(t,x) [dNdtprey(x(1), x(2)); dPdtpredator(x(1), x(2))];
N0 = 0.1;
P0 = 0.1;
[~, trajectory] = ode45(functionh, [0,100], [N0;P0]);
plot(trajectory(:, 1), trajectory(:, 2), 'LineWidth', 2, 'DisplayName', 'trajectory');

% Fiexed Points 
plot(7,0, "ko", "MarkerSize", 7,'DisplayName','Saddle node');
plot(2.1926, 2.1926,"k.", "MarkerSize", 20, 'DisplayName','Stable point');
title('Portrait fase');
xlabel('N');
ylabel('P');
xlim([-1, 10]);
ylim([-1, 10]);

%% Jacobian matrix. 
% Parameters 
r1 = 1;
r2 = 0.1; 
K = 7; 
d = 1;
j = 1; 
w = 1; 

%% For the fixed point N = 7 and  P = 0. 
 N = 7; 
 P = 0; 

JacoMax4 = [r1-d*w*(P/((d+N)^2))-2*r1*(N/K), -(N*w)/(d+N);
    (j*(P^2)*r2)/(N^2), r2*(1-(2*j*P)/N)];  
[EigVec4 EgiVal4] = eig(JacoMax4)

% Saddle node. 
 
%% For the fixed point N = 39/20 + (sqrt(4321)/20) and P = 39/20 + (sqrt(4321)/20);
N = (1/2)*(sqrt(29)-1); 
P = (1/2)*(sqrt(29)-1);

JacoMax5 = [r1-d*w*(P/((d+N)^2))-2*r1*(N/K), -(N*w)/(d+N);
    (j*(P^2)*r2)/(N^2), r2*(1-(2*j*P)/N)];  
[EigVec5 EgiVal5] = eig(JacoMax5)

% Oscillatory behaviour