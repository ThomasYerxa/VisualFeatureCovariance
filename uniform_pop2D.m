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
    
    % Visualize by adding the 2-D tuning curves together. (And fisher info
    % sums).
    res    = zeros(size(X));
    fisher = zeros(size(X));
    for i = 1:N_x
        for j = 1:N_y
            % print i,j for progress check
            [i, j]
            % define variables and 2D gaussian tuning curve
            syms x y
            f(x, y) = exp(-0.5 * ((x - m_x(i))^2 + (y-m_y(j))^2)/(0.55^2)) / (N_x*N_y);
            % compute each partial derivative
            dfx     = diff(f, x);
            dfy     = diff(f, y);
            % calculate nth tuning curve and its fisher information
            c_n     = double(f(X, Y));
            c_n_dx  = double(dfx(X, Y));
            c_n_dy  = double(dfy(X, Y));
            fish_n  = (c_n_dx.*c_n_dy).^2 ./ (c_n.^2);
            % sum results. 
            res     = res + c_n;
            fisher  = fisher + fish_n; 
        end
    end
    imagesc(res); colorbar; 
    axis image off;
end

fig = figure; 
set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig, 'name', 'Uniform Population 2-D Fisher Info');
imagesc(fisher); axis image off; 

end