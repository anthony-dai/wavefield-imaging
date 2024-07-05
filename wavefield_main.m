%%%
% Author: Anthony Dai
% Course: Wavefield imaging EE4595
% Inverse scattering problem (medium unkown, source, wavefield known at receiver locations)
%%%

close all;
clear all;
clc;

%% 
% Parameters
kb = 1;  % Given kb = 1
lambda = 2 * pi / kb;  % Wavelength of the background field
M = 20;  % Number of receivers
% Source location
rho_s = [(lambda / 2), 10 * lambda];

% Create figure
figure;
hold on;
axis equal;
grid on;

% Draw the rectangle
rectangle('Position', [0, 0, lambda, lambda], 'EdgeColor', 'b', 'LineWidth', 1.5);

% Define the end points of the line segment
x1 = -lambda;  % x-coordinate of the first end point
y1 = 1.5*lambda; % y-coordinate of the first end point
x2 = 2*lambda;   % x-coordinate of the second end point
y2 = 1.5*lambda; % y-coordinate of the second end point

% Plot the receiver line segment
plot([x1 x2], [y1 y2], 'r-', 'LineWidth', 1);

% Create a vector of uniformly spaced grid points along the line segment
t = linspace(0, 1, M); % Parameter t ranges from 0 to 1
x_receivers = x1 + t * (x2 - x1); % Interpolate x-coordinates
y_receivers = y1 + t * (y2 - y1); % Interpolate y-coordinates
scatter(x_receivers, y_receivers, 50, 'g', 'x', 'LineWidth', 1.5); % Grid points

% Mark the source location
plot(rho_s(1), rho_s(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Add labels
xlabel('x');
ylabel('y');
title('Sketch of the Configuration');

% Add text annotations for clarity
text(rho_s(1), rho_s(2), ' \rho_s', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Set axis limits
xlim([-2*lambda, 3*lambda]);
ylim([0, 11*lambda]);

legend('Receiver line');


% Invert y-axis to have positive y going down
set(gca, 'YDir', 'reverse');



hold off;

%%
% Define the step size
h = lambda / 20;

% Compute the number of grid points in the x-direction
num_x_points = lambda / h;

% Compute the number of grid points in the y-direction
num_y_points = lambda / h;

% Compute the total number of grid points
N = num_x_points * num_y_points;

% Display the total number of grid points
disp(['Total number of grid points: ', num2str(N)]);

%% Grid
% Define the dimensions of the rectangle
width = lambda;  % Width of the rectangle
height = lambda;  % Height of the rectangle

% Generate the grid of points within the rectangle
x = h/2:h:width-h/2;  % X-coordinates of the points
y = h/2:h:height-h/2; % Y-coordinates of the points
[X, Y] = meshgrid(x, y);  % Create a 2D grid of points

% Flatten the grid matrices to create a list of points
X = X(:);
Y = Y(:);

% Create a figure
figure;

% Plot the points
scatter(X, Y, 'filled','black');

% Draw the rectangle
hold on;  % Hold on to add the rectangle to the same plot
rectangle('Position', [0 0 width height], 'EdgeColor', 'b', 'LineWidth', 2);

% Set the axis properties
axis equal;
axis([0 width 0 height]);  % Set the axis limits to match the rectangle dimensions
xlabel('X-axis');
ylabel('Y-axis');
title('Points Spread with Step Size h inside object D');
% Invert y-axis to have positive y going down
set(gca, 'YDir', 'reverse');
% Add grid for better visualization
grid on;

%%
% Compute the Hankel function of the second kind and order zero
u_inc = zeros(num_x_points, num_y_points);
for x = 1:num_x_points
    for y = 1:num_y_points
        u_inc(y,x) = -1j/4 * besselh(0, 2, kb * norm( [h/2 + (x-1)*h,  h/2 + (y-1)*h] - rho_s ));
    end
end

% Compute the real, imaginary, and absolute value of u_inc
real_part = real(u_inc);
imaginary_part = imag(u_inc);
absolute_value = abs(u_inc);

% Plot the real part of u_inc
figure;
imagesc(real_part);
colorbar;
title('Real Part of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

% Plot the imaginary part of u_inc
figure;
imagesc(imaginary_part);
colorbar;
title('Imaginary Part of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

% Plot the absolute value of u_inc
figure;
imagesc(absolute_value);
colorbar;
title('Absolute Value of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

%% Move source closer to object
% Source location
rho_s = [(lambda / 2), 10 * lambda];

% Move closer to object
rho_s_closer = [(lambda / 2), lambda];

% Compute the Hankel function of the second kind and order zero
u_inc = zeros(num_x_points, num_y_points);
for x = 1:num_x_points
    for y = 1:num_y_points
        u_inc(y,x) = -1j/4 * besselh(0, 2, kb * norm( [h/2 + (x-1)*h,  h/2 + (y-1)*h] - rho_s_closer ));
    end
end

% Compute the real, imaginary, and absolute value of u_inc
real_part = real(u_inc);
imaginary_part = imag(u_inc);
absolute_value = abs(u_inc);

% Plot the real part of u_inc
figure;
imagesc(real_part);
colorbar;
title('Real Part of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

% Plot the imaginary part of u_inc
figure;
imagesc(imaginary_part);
colorbar;
title('Imaginary Part of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

% Plot the absolute value of u_inc
figure;
imagesc(absolute_value);
colorbar;
title('Absolute Value of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

%% Increase kb by factor

% Compute the Hankel function of the second kind and order zero
u_inc = zeros(num_x_points, num_y_points);
for x = 1:num_x_points
    for y = 1:num_y_points
        u_inc(y,x) = -1j/4 * besselh(0, 2, 2 * kb * norm( [h/2 + (x-1)*h,  h/2 + (y-1)*h] - rho_s_closer ));
    end
end

% Compute the real, imaginary, and absolute value of u_inc
real_part = real(u_inc);
imaginary_part = imag(u_inc);
absolute_value = abs(u_inc);

% Plot the real part of u_inc
figure;
imagesc(real_part);
colorbar;
title('Real Part of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

% Plot the imaginary part of u_inc
figure;
imagesc(imaginary_part);
colorbar;
title('Imaginary Part of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

% Plot the absolute value of u_inc
figure;
imagesc(absolute_value);
colorbar;
title('Absolute Value of u\_inc');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');

%% contrast function
rho_s = [(lambda / 2), 10 * lambda];
kb = 1;

% % Initialize a 20x20 matrix of zeros
% chi = zeros(20, 20);
% 
% % Define the center of the grid
% center = 10;
% 
% % Create the vertical part of the plus sign
% chi(center-5:center+5, center) = 1;
% 
% % Create the horizontal part of the plus sign
% chi(center, center-5:center+5) = 1;

load('chi_1bigdot_1smalldot.mat')
% Display the grid using imagesc
figure;
imagesc(chi);
title('Contrast function')
% Add a colorbar for better visualization
colorbar;

% Ensure the axes are set to equal for correct aspect ratio
axis equal;
axis tight;

%%
delta_x = h; % width of one pixel
delta_y = h; % height of one pixel

A = zeros(M, num_x_points*num_y_points); % M x N system matrix, M : number of receivers, N = total samples

% Data equation using the midpoint rule
for m = 1:M
    n = 1;
    for i = 1:num_x_points
        for j = 1:num_y_points
            % x,y are the midpoints of each pixel 
            x = delta_x/2 + (i-1)*delta_x; 
            y = delta_y/2 + (j-1)*delta_y;
            u_tmp = -1j/4 * besselh(0, 2, kb * norm( [x,y] - rho_s ));
            G_tmp = -1j/4 * besselh(0, 2, kb * norm( [x_receivers(m), y_receivers(m)] - [x,y] ));
            
%             a = 0.5 * min([delta_x, delta_y]);
%             G_tmp = -1j/(2 * kb * a) * besselh(1, 1, kb * a) * besselh(0, 2, kb * norm( [x_receivers(m), y_receivers(m)] - [x,y] )); %weak green
            A(m,n) = kb^2 * delta_x * delta_y * G_tmp * u_tmp;
            n = n + 1;
        end
    end
end

% reshape contrast function
chi_vec = reshape(chi, [], 1);
%%
% Data equation using the midpoint rule
% G_matrix = zeros(num_x_points*num_y_points);

% for q = 1:num_y_points
%     for j = 1:num_y_points
%         G_block = zeros(num_x_points);
%         for p = 1:num_x_points
%             for i = 1:num_x_points
%                 % x,y are the midpoints of each pixel 
%                 x = delta_x/2 + (i-1)*delta_x; 
%                 y = delta_y/2 + (j-1)*delta_y;
%                 % x_p, x_q are the midpoints of each pixel given by p and q
%                 x_p = delta_x/2 + (p-1)*delta_x; 
%                 y_q = delta_y/2 + (q-1)*delta_y;
%                 
%                 a = 0.5 * min([delta_x, delta_y]);
%                 
%                 distance = norm( [x_p, y_q] - [x,y] );
%                 
%                 if distance > 0
%                     G_block(p,i) = -1j/(2 * kb * a) * besselh(1, 1, kb * a) * besselh(0, 2, kb * distance); %weak green
%                 else
%                     G_block(p,i) = -1j/(2 * kb * a) * (besselh(1, 2, kb * a) - 2j/(pi * kb * a) ); %weak green when distance == 0
%                 end
%             end
%         end
%         G_matrix((q-1)*num_x_points + 1:q*num_x_points, (j-1)*num_x_points + 1:j*num_x_points) = kb^2 * delta_x * delta_y * G_block;
%     end
% end
% identity_matrix = eye(num_x_points*num_y_points);
% chi_diag = diag(chi(:));
% A = G*C;
%% svd
singular_values = svd(A);
figure;
stem(singular_values, 'filled');
title('Singular Values of Matrix A');
xlabel('Index');
ylabel('Singular Value');

%% Compute u_sc
u_sc = A * chi_vec;

%% min norm solution for x
% Add noise
% Specify the noise level
noise_level = 1*10^(-7);
% noise_level = 0;

% Generate random noise
noise = noise_level * rand(size(u_sc));

% Load fixed noise at e-7
% load('noise_eminus7.mat');

% Add the noise to the vector
u_sc_noisy = u_sc + noise;

% using pseudo inverse
x_min_norm = pinv(A) * u_sc_noisy;

% % svd
% [U, S, V] = svd(A);
% S = S(1:M, 1:M);
% U = U(:,1:M);
% V = V(:,1:M);
% x_min_norm_svd = real(V*inv(S)*U'*u_sc_noisy);

% truncated svd
[U, S, V] = svd(A);
S = S(1:M, 1:M);
U = U(:,1:M);
V = V(:,1:M);
x_min_norm_tsvd = 0;
r_threshold = 3e-5;
singular_values_retained = 0;
for r = 1:size(S,1)
    if S(r,r) > r_threshold
        x_min_norm_tsvd = x_min_norm_tsvd + 1 / S(r,r) * (U(:,r)' * u_sc_noisy) * V(:,r);
        singular_values_retained = singular_values_retained + 1;
    end
end

% Reshape vector into a m-by-n matrix
X_min_norm = reshape(abs(x_min_norm_tsvd), [num_x_points, num_y_points]);

% Plot X_min_norm
figure;
imagesc(X_min_norm);
colorbar;
title('Image reconstruction');
axis equal tight;
xlabel('X-axis');
ylabel('Y-axis');
