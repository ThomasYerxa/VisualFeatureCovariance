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
%  
% Returns: 
% ---------
% Returns None
%-------------------------------------------------------------------------


function uniform_pop2D(curve_shape, stim_range, N_x, N_y)


% gaussian tuning curves
if curve_shape == 0
    
    fig = figure; hold on; 
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig, 'name', 'Uniform Population 2-D');
    
    % Define 2D axes
    s_x    = linspace(stim_range(1), stim_range(2), 100);
    s_y    = linspace(stim_range(3), stim_range(4), 100);
    [X, Y] = meshgrid(s_x, s_y);
    
    % Define the centers convolutional population 
    m_x    = linspace(stim_range(1), stim_range(2), N_x);
    m_y    = linspace(stim_range(3), stim_range(4), N_y);
    
    % Visualisation is achieved by adding the 2-D tuning curves together.
    res    = zeros(size(X));
    for i = 1:N_x
        for j = 1:N_y
            h_n     = gaussian2D(m_x(i), m_y(j), N_x * N_y);
            curve_n = h_n(X, Y);
            res     = res + curve_n;     
        end
    end
    imagesc(res); colorbar; 
    axis image off;
end


end