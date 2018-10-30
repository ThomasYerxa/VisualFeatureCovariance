% this function visualizes a convolutional population parameterized by a 
% density function. See equation 2.10 in Ganguli and Simoncelli, 2014

function func_density_pop1D(curve_shape, stim_range, n_neurons, D)


% gaussian tuning curve
if curve_shape == 0
    std = 1; 
    
    fig = figure; hold on; 
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig, 'name', 'Density Function Population 1-D');
    
    ylabel('firing rate');
    xlabel('s');
    
    
    s = linspace(stim_range(1), stim_range(2), 100);
    for n = 1:n_neurons
        h_n     = gaussian1D(0, std, n_neurons);
        curve_n = h_n(D(s) - n);
        plot(s, curve_n);
    end
    
end

end