% this function visualizes a convolutional population parameterized by a 
% density function. See equation 2.10 in Ganguli and Simoncelli, 2014

function func_density_pop1D(curve_shape, stim_range, n_neurons)


% gaussian tuning curve
if curve_shape == 0
    sigma    = 0.55;  
    fig = figure; hold on; 
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig, 'name', 'Density Function Population 1-D');
    
    s = linspace(stim_range(1), stim_range(2), 200);
    
    syms x 
    D(x) = (n_neurons/2)*(1+erf(x/(0.5*sqrt(2))));
    d(x) = diff(D, x);
    subplot(4,1,1); hold on; 
    plot(s, D(s)); 
    ylabel('Cumulative Density');
    xlabel('s');
    
    
    for n = 1:n_neurons
        h_n     = gaussian1D(0, sigma, n_neurons);
        curve_n = h_n(D(s) - n);
        subplot(4,1,2); hold on; 
        xlabel('s');
        ylabel('firing rate')
        plot(s, curve_n);
    end
    
    
    
    res = zeros(size(s));
    for n = 1:n_neurons
        h_n       = gaussian1D(0, sigma, n_neurons);
        curve_n   = h_n(D(s) - n);
        curve_n_p = (n-D(s)) .* d(s) .* curve_n;
        res       = res + (curve_n_p .^ 2 ./ curve_n);
        ax4       = subplot(4,1,4); 
        plot(ax4, s, (curve_n_p .^ 2 ./ curve_n)); hold on; 
        
    end
    
    ax2       = subplot(4,1, 3);
    
    
    plot(ax2, s, res); hold on;
   
    ylabel('Fisher Information');
    xlabel('s');
    
end

end