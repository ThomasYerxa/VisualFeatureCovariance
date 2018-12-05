% figure for 1-D population parameterized by a displacement field 

% Define Parameters
std        = 1.0; 
stim_range = [-8, 8];
n_neurons  = 9;

% Define set of axes
fig = figure; hold on; 
set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig, 'name', 'Displacement Population 1-D');
ax1 = subplot(3,2, 1); hold on; 
ax2 = subplot(3,2, 2); hold on; 
ax3 = subplot(3,2, 3); hold on; 
ax4 = subplot(3,2, 4); hold on; 
ax5 = subplot(3,2, 5); hold on; 
ax6 = subplot(3,2, 6); hold on; 
    
% Define axis and neuron positions. Initialize Fisher Info
s             = linspace(stim_range(1), stim_range(2), 200);
m             = linspace(stim_range(1), stim_range(2), n_neurons);
fisher        = zeros(size(s));
fisher_warped = zeros(size(s));


% Define Gaussian tuning curve and its derivative
syms x
h(x)   = exp(-0.5 * (x^2)/(std^2)) / (n_neurons);
h_p(x) = diff(h, x);

% Define displacement field
slope = 0.8;
syms x
f(x)  = slope * sqrt(abs(x));
df(x) = diff(f, x); 


for n = 1:n_neurons
    % nth tuning curve and corresponding fisher information
    curve_n  = double(h(s-m(n)));
    deriv_n  = double(h_p(s-m(n)));
    fisher_n = deriv_n.^2 ./ curve_n; 
    fisher   = fisher + fisher_n; 
    
    % warped curves and corresponding fisher information
    curve_warped_n = double(h(s - f(s) - m(n)));
    deriv_warped_n = (1 - double(df(s))) .* double(h_p(s - f(s) - m(n)));
    
    fisher_warped_n = deriv_warped_n.^2 ./ curve_warped_n; 
    fisher_warped   = fisher_warped + fisher_warped_n; 
    
    % plot both sets of tuning curves
    plot(ax3, s, curve_n); 
    plot(ax4, s, curve_warped_n);
    
end
xlabel(ax3, 's');
ylabel(ax3, 'Firing Rate');
xlabel(ax4, 's');
ylabel(ax4, 'Firing Rate');


% plot each displacement field
plot(ax1,s, zeros(size(s))); 
xlabel(ax1, 's');
ylabel(ax1, 'Displacement Field');

plot(ax2, s, double(f(s))); 
xlabel(ax2, 's');
ylabel(ax2, 'Displacement Field');

% plot both sets of fisher information.
I_conv = mean(fisher(30:170)); 
plot(ax5, s, fisher);
plot(ax5, s, ones(size(s)).*I_conv);
legend(ax5, {'Measured Fisher Information','Approximate Fisher Information'},'Location','southwest');
xlabel(ax5, 's');
ylabel(ax5, 'Fisher Information');
plot(ax6, s, fisher_warped); 
plot(ax6, s, ones(size(s)) .* I_conv.*double((1-df(s))).^2)  
legend(ax6, {'Measured Fisher Information','Approximate Fisher Information'},'Location','northwest');
xlabel(ax6, 's');
ylabel(ax6, 'Fisher Information');







