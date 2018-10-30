% this function visualizes a convolutional population parameterized by a 
% scalar density. See equation 2.9 in Ganguli and Simoncelli, 2014

function scalar_density_pop1D(curve_shape, stim_range, n_neurons, d)


% gaussian tuning curves
if curve_shape == 0
    std = 1; 
    
    fig = figure; hold on; 
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig, 'name', 'Uniform Population 1-D');
    ylabel('firing rate');
    xlabel('s');
    
    s = linspace(stim_range(1), stim_range(2), 100);
    for n = 1:n_neurons
        h_n = gaussian1D(0, std, n_neurons);
        curve_n = h_n(d .* (s - d/n));
        plot(s, curve_n); 
    end
    
end

end