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
%                [s_min, s_max]
%  n_neurins   : Number of neurons that span the s-axis
%  f           : A function of 1 variable that defines how the domain will
%                be displaced. 
%  
% Returns: 
% ---------
% Returns None
%-------------------------------------------------------------------------



function displacement_pop1D(curve_shape, stim_range, n_neurons, f)


% gaussian tuning curve
if curve_shape == 0
    std = 1; 
    
    fig = figure; hold on; 
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig, 'name', 'Displacement Function Population 1-D');
    
    ylabel('firing rate');
    xlabel('s');
    
    % define axis and neuron positions
    s      = linspace(stim_range(1), stim_range(2), 100);
    means  = linspace(stim_range(1), stim_range(2), n_neurons);
    for n = 1:n_neurons
        h_n     = gaussian1D(means(n), std, n_neurons);
        % visualize tuning cuver over displaced domain
        curve_n = h_n(s + f(s));
        plot(s, curve_n);
    end
    
end

end

