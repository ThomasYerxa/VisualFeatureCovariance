% reproduce figure 1 from Ganguli and Simoncelli

fig = figure; hold on; 
set(fig, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig, 'name', 'Density Population 1-D');
ax1 = subplot(3,2, 1); hold on; 
ax2 = subplot(3,2, 2); hold on; 
ax3 = subplot(3,2, 3); hold on; 
ax4 = subplot(3,2, 4); hold on; 
ax5 = subplot(3,2, 5); hold on; 
ax6 = subplot(3,2, 6); hold on; 

% define parameters 
std        = 0.55; 
stim_range = [-4, 4];
n_neurons  = 9;
   
% Define axis and neuron positions (for uniform population)
s      = linspace(stim_range(1), stim_range(2), 200);
m      = linspace(stim_range(1), stim_range(2), n_neurons);

% Define Gaussian Tuning curve and its derivative
syms x
h(x)   = exp( (x^2)/(-2*std^2)) / n_neurons;
h_p(x) = diff(h, x);

% Define density function and its cumulative. Note: Ther are alternative
% implementations of D and d commented out, I'm not sure why they aren't
% all identical. 
std_d = 1.0;
syms x 
%D(x) = (n_neurons/2)*( 1+erf( x/(std_d*sqrt(2)) ) );
%d    = diff(D, x);
d(x)  = n_neurons * exp(-x^2/(2*std_d^2))/sqrt(2*pi*std_d^2);
D(x)  = (n_neurons/2)*( 1+erf( x/(std_d*sqrt(2)) ) );
%d(x)  = n_neurons * exp(x^2/(2*std_d^2))/sqrt(2*pi*std_d^2)
%D     = int(d, x)

% Define gain function
syms x
g(x) = 1.0 + abs(x)/3; 

% Find perferred values of warped tuning curves
s_n = zeros(1, n_neurons);
for n = 1:n_neurons
    for x = 1:size(s, 2)
        if double(D(s(x))) - n < .01
            s_n(n) = s(x);
        end
    end
end
            


fisher_uniform = zeros(size(s));
fisher_warped  = zeros(size(s));
for n = 1:n_neurons
    % Uniform tuning curves and corresponding fisher information
    h_n      = double(h(s - m(n)));
    h_n_p    = double(h_p(s- m(n)));
    fisher_n = h_n_p.^2 ./ h_n;
    
    % Tuning curves warped by D and corresponding fisher information
    hw_n      = double(g(s_n(n))) .* double(h(D(s)- n)); 
    hw_n_p    = double(g(s_n(n))) .* d(s) .* double(h_p(D(s)- n)); 
    fisherw_n = hw_n_p.^2 ./ hw_n; 
    
    % Sum respective Fisher infos
    fisher_uniform = fisher_uniform + fisher_n; 
    fisher_warped  = fisher_warped + fisherw_n; 
    
    % plot both sets of tuning curves
    plot(ax3, s, h_n); 
    axis(ax3, [-2.5 2.5 -inf inf])
    plot(ax4, s, hw_n); 
    axis(ax4, [-2.5 2.5 -inf inf])
end

xlabel(ax3, 's');
ylabel(ax3, 'Firing Rate');
xlabel(ax4, 's');
ylabel(ax4, 'Firing Rate');


% plot D(s) for each case
n_s = size(s, 2);
linear_CDF = linspace(0,n_neurons, n_s);
plot(ax1, s, linear_CDF);
xlabel(ax1, 's');
ylabel(ax1, 'Cumulative Density');
plot(ax2, s, double(D(s)));
xlabel(ax2, 's');
ylabel(ax2, 'Cumulative Density');

% plot both fisher information curves
plot(ax5, s, fisher_uniform); 
xlabel(ax5, 's');
ylabel(ax5, 'Fisher Information');
plot(ax6, s, fisher_warped); 
xlabel(ax6, 's');
ylabel(ax6, 'Fisher Information');

