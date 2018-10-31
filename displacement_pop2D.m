%-------------------------------------------------------------------------
% This function visualizes a uniform convolutional population of tuning 
% curves in 2 dimensions
% 
% Parameters: 
% -----------
%  curve_shape : code for which tuning curve to use; currently only
%                gaussian tuning curves are implemented and this must be
%                equal to zero
%  stim_range  : specifies the range of stimulus values. is of the form 
%                [s_x_min, s_x_max, s_y_min, s_y_max]
%  N_x         : Number of neurons that span the x-axis
%  N_y         : Numver of neurons that span the y-axis
%  f_x         : A function of 2 variables that defines how the domain will
%                be displaced in the x-direction. 
%  f_y         : A function of 2 variables that defines how the domain will
%                be displaced in the y-direction. 
%  
% Returns: 
% ---------
% Returns None
%-------------------------------------------------------------------------

function displacement_pop2D(curve_shape, stim_range, N_x, N_y, f_x, f_y)


% gaussian tuning curves
if curve_shape == 0
    
    
    
    % define 2-D axis. Pad so that when we visualize the displaced domain
    % we are sure to see all of the relevant behavior
    range_x = stim_range(2) - stim_range(1);
    range_y = stim_range(4) - stim_range(3);
    s_x     = linspace(stim_range(1) - range_x / 2, stim_range(2) + range_x / 2, 100);
    s_y     = linspace(stim_range(3) - range_y / 2, stim_range(4) + range_y / 2, 100);
    [X, Y]  = meshgrid(s_x, s_y);
    
    % define the centers of the convolutional population
    m_x    = linspace(stim_range(1), stim_range(2), N_x);
    m_y    = linspace(stim_range(3), stim_range(4), N_y);
    
    res    = zeros(size(X));
    uni    = zeros(size(X));
    for i = 1:N_x
        for j = 1:N_y
            % Generate curve centered at appropriate point
            h_n     = gaussian2D(m_x(i), m_y(j), N_x * N_y);
            % The curve of interest is h_n as a function over the displaced
            % domain. 
            curve_n = h_n(X + f_x(X, Y), Y + f_y(X, Y));
            curve_u = h_n(X, Y);
            res     = res + curve_n; 
            uni     = uni + curve_u;
             
        end
    end
    
    fig = figure; hold on; 
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig, 'name', 'Displacement Population 2-D');
    
    subplot(1, 3, 1); hold on; 
    imagesc(uni); 
    title('Uniform Tuning Curves'); 
    axis image off; colorbar
    
    subplot(1,3,2); hold on; 
    x_q        = linspace(stim_range(1), stim_range(2), 20);
    y_q        = linspace(stim_range(3), stim_range(4), 20);
    [X_q, Y_q] = meshgrid(x_q, y_q);
    
    title('Negative Displacement Field');
    xlabel('stimulus x');
    ylabel('stimulus y');
    quiver(X_q, Y_q, -1.*f_x(X_q, Y_q), -1.*f_y(X_q, Y_q));
    axis('square'); 
    
    subplot(1, 3, 3); hold on;
    title('Warped Tuning Curves'); 
    imagesc(res); colorbar 
    axis image off; 
    
    
end


end