clc;
clear;

% Create grid
[x, y] = meshgrid(0:0.5:10, 0:0.5:10);

% Define a time-dependent velocity field
u = @(x, y, t) -y + 0.1*t;  % x-component of velocity
v = @(x, y, t)  x + 0.1*t;  % y-component of velocity


% Time range for simulation
tspan = linspace(0, 10, 100);

% Starting point
x0 = 5; 
y0 = 5;

% --- PATHLINE ---
% Trajectory of a single fluid particle over time
path = zeros(length(tspan), 2);
path(1,:) = [x0, y0];
for i = 2:length(tspan)
    dt = tspan(i) - tspan(i-1);
    x_prev = path(i-1,1);
    y_prev = path(i-1,2);
    path(i,1) = x_prev + u(x_prev, y_prev, tspan(i-1)) * dt;
    path(i,2) = y_prev + v(x_prev, y_prev, tspan(i-1)) * dt;
end

% --- STREAMLINE ---
% Instantaneous lines tangent to the velocity vectors (at t = 0)
[xs, ys] = meshgrid(0:1:10, 0:1:10);
us = u(xs, ys, 0);   % Velocity field at time t = 0
vs = v(xs, ys, 0);

% --- STREAKLINE ---
% Locus of particles that have passed through a specific point (x0,y0)
streak = zeros(length(tspan), 2);
for i = 1:length(tspan)
    % Time when the particle is injected
    ti = tspan(i);
    x_temp = x0;
    y_temp = y0;
    for j = i:length(tspan)-1
        dt = tspan(j+1) - tspan(j);
        x_temp = x_temp + u(x_temp, y_temp, tspan(j)) * dt;
        y_temp = y_temp + v(x_temp, y_temp, tspan(j)) * dt;
    end
    streak(i,:) = [x_temp, y_temp];
end

% --- PLOT RESULTS ---
figure;
hold on;
quiver(xs, ys, us, vs, 'k');                           % Velocity field vectors
plot(path(:,1), path(:,2), 'r', 'LineWidth', 2);       % Pathline in red
streamline(xs, ys, us, vs, x0, y0);                    % Streamline (at t=0)
plot(streak(:,1), streak(:,2), 'b--', 'LineWidth', 2); % Streakline in blue dashed

% Annotations
legend('Velocity Field','Pathline','Streamline','Streakline');
title('Streamline vs Pathline vs Streakline');
xlabel('x');
ylabel('y');
axis equal;
grid on;
