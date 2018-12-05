%-------------------------------------------------------------------------
% This function visualizes convolutional population of tuning 
% curves (and their corresponding fisher information curves)in 2 dimensions
% 
%-------------------------------------------------------------------------

function displacement_pop2D

% define parameters
std        = 1.0;
stim_range = [-6 6 -6 6];
N_x        = 8; 
N_y        = 8;

% define set of axes for plotting
fig = figure; hold on; 
set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig, 'name', 'Displacement Population 1-D');
ax1 = subplot(3,2, 1); hold on; 
ax2 = subplot(3,2, 2); hold on; 
ax3 = subplot(3,2, 3); hold on; 
ax4 = subplot(3,2, 4); hold on; 
ax5 = subplot(3,2, 5); hold on; 
ax6 = subplot(3,2, 6); hold on; 

% define 2-D stimulus space. Pad so that when we visualize the displaced domain
% we are sure to see all of the relevant behavior
range_x = stim_range(2) - stim_range(1);
range_y = stim_range(4) - stim_range(3);
s_x     = linspace(stim_range(1) - range_x / 2, stim_range(2) + range_x / 2, 100);
s_y     = linspace(stim_range(3) - range_y / 2, stim_range(4) + range_y / 2, 100);
[X, Y]  = meshgrid(s_x, s_y);
    
% define the centers of the convolutional population
m_x    = linspace(stim_range(1), stim_range(2), N_x);
m_y    = linspace(stim_range(3), stim_range(4), N_y);

fisher   = meshgrid(zeros(size(s_x)), zeros(size(s_y)));
fisher_w = meshgrid(zeros(size(s_x)), zeros(size(s_y)));

% Define 2D Gaussian Tuning Curves and derivatives in each direction
syms x y
h(x, y)   = exp(-0.5 * (x^2 + y^2)/(std^2)) / (N_x * N_y);
h_x(x, y) = diff(h, x);
h_y(x, y) = diff(h, y);

% Define displacement field and all derivatives
syms x y
f_x(x, y)  = sqrt(abs(x * y));
f_y(x, y)  = x;
f_xx(x, y) = diff(f_x, x);
f_xy(x, y) = diff(f_x, y); 
f_yx(x, y) = diff(f_y, x);
f_yy(x, y) = diff(f_y, y); 

  

for i = 1:N_x
    for j = 1:N_y
        % Generate curve centered at appropriate point
        curve_n = double(h(X - m_x(i), Y - m_y(j)));
        dx_n    = double(h_x(X - m_x(i), Y - m_y(j)));
        dy_n    = double(h_y(X - m_x(i), Y - m_y(j)));
        
        % elements of fisher info matrix
        fish_xx = dx_n.^2 ./ curve_n;
        fish_xy = dx_n.*dy_n ./ curve_n;
        fish_yx = dx_n.*dy_n ./ curve_n;
        fish_yy = dy_n.^2 ./ curve_n;
        % determinant of said matrix
        fish_n  = fish_xx.*fish_yy - fish_yx.*fish_xy;
        
        % Identical steps over warped domain (_w <-> warped). Note that we
        % use the total derivative both arguments now depend on X and Y. 
        curve_n_w = double(h(X - f_x(X,Y)- m_x(i), Y - f_y(X,Y) - m_y(j)));
        
        dx_n_w    = (1 - double(f_xx(X,Y))) .* double(h_x(X-f_x(X,Y)- m_x(i), Y - f_y(X,Y) - m_y(j))) ...
            - double(f_yx(X,Y)) .* double(h_y(X-f_x(X,Y)- m_x(i), Y - f_y(X,Y) - m_y(j)));
        
        dy_n_w    = (1 - double(f_yy(X,Y))) .* double(h_y(X-f_x(X,Y)- m_x(i), Y - f_y(X,Y) - m_y(j))) ...
            - double(f_xy(X,Y)) .* double(h_x(X-f_x(X,Y)- m_x(i), Y - f_y(X,Y) - m_y(j)));
        
        % elements of fisher info matrix (warped)
        fish_xx_w = dx_n_w.^2    ./ curve_n_w;
        fish_xy_w = dx_n_w.*dy_n ./ curve_n_w;
        fish_yx_w = dx_n_w.*dy_n ./ curve_n_w;
        fish_yy_w = dy_n_w.^2    ./ curve_n_w;
        % determinant of said matrix
        fish_n_w  = fish_xx_w.*fish_yy_w - fish_yx_w.*fish_xy_w;
        
        % fisher info of each curve sums. 
        fisher   = fisher + fish_n;
        fisher_w = fisher_w + fish_n_w; 
        % plot both sets of curves. 
        surf(ax3, X, Y, curve_n);
        surf(ax4, X, Y, curve_n_w);
             
    end
end


% less dense grid for better visualizatio of displacement field
x_q        = linspace(stim_range(1), stim_range(2), 20);
y_q        = linspace(stim_range(3), stim_range(4), 20);
[X_q, Y_q] = meshgrid(x_q, y_q);
    
quiver(ax2, X_q, Y_q, meshgrid(zeros(size(X)), zeros(size(Y))));
quiver(ax2, X_q, Y_q, double(f_x(X_q, Y_q)), double(f_y(X_q, Y_q)));

surf(ax5, X, Y, fisher);
surf(ax6, X, Y, fiser_w);
    
end

