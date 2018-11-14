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
%  n_neuron    : Number of neurons that span the s-axis
%  
% Returns: 
% ---------
% Returns None
%-------------------------------------------------------------------------


function uniform_pop1D(curve_shape, stim_range, n_neurons)


% gaussian tuning curves
if curve_shape == 0
    std = 1; 
    
    fig = figure; hold on; 
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig, 'name', 'Uniform Population 1-D');
    ax1 = subplot(2,1,1); hold on; 
    ylabel(ax1, 'firing rate');
    xlabel(ax1, 's');
    
    % Define axis and neuron positions
    s      = linspace(stim_range(1), stim_range(2), 200);
    means  = linspace(stim_range(1), stim_range(2), n_neurons);
    for n = 1:n_neurons
        h_n     = gaussian1D(means(n), std, n_neurons);
        curve_n = h_n(s);
        plot(ax1, s, curve_n); 
        axis([-6 6 -inf inf]);
    end
    
    % separately plot the fisher information for uniform population

    res = zeros(size(s));
    for n = 1:n_neurons
        h         = gaussian1D(0, std, n_neurons);
        curve_n   = h(s - means(n));
        curve_n_p = (means(n)-s) .* curve_n; 
        res       = res + (curve_n_p .^ 2 ./ curve_n); 
        
    end
    
    ax2 = subplot(2,1, 2);
    res_trunc = res(25:175); 
    plot(ax2, s, ones(size(s)).*mean(res_trunc)); hold on; 
    plot(ax2, s, res); hold on;
    legend(ax2, 'Approximate Fisher Information', 'Actural Fisher Information', 'Location', 'southwest');
    ylabel('Fisher Information');
    xlabel('s');
    axis([-6 6 -inf inf]);
    
    
    
   
     
    
end




end