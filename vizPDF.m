%-------------------------------------------------------------------------
% Function for visualizing probability distributions.
% 
% Parameters: 
% -----------
%  PDF_j : the measured joint distribution. size: (n_orientation, n_f)
%  PDF_t : the measured distribution over theta. size: (n_orientation)
%  PDF_f : the measured distribution over f. size:  (n_f)
%  PDF_s : the 'separable 'joint distribution. size: (n_orientation, n_f)
%  f : the spatial frequencies used. size: (n_f), units: 1/pixels
%  orientation : the orientations used. size: (n_orientation), units:
%  radians
%
% Returns: 
% ---------
% Returns 0. (Not for use)
%-------------------------------------------------------------------------

function a = vizPDF(PDF_j, PDF_t, PDF_f, PDF_s, f,orientation)

% set return value (arbitrary)
a = 0;
% convert orientation to degrees
orientation = orientation * 360 /(2*pi); 
% truncate f to 3 decimal places for display
f = f - rem(f, 0.001);
% set logScale
logScale = false; 

% generate tick positions for 2-D plots. One tick per sampled value. 
nX     = size(f, 2);
nY     = size(orientation, 2);
xticks = [0.5];
yticks = [0.5];
for x = 2:nX
    xticks(x) = xticks(x-1) + 1.0;
end
for y = 2:nY
    yticks(y) = xticks(y-1) + 1.0;
end
% ---- FIRST FIGURE: Measured Joint PDF and outer product PDF side by side
fig_joint = figure;

set(fig_joint, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig_joint, 'name', 'Joint Distributions');

ax1 = subplot(1,2,1); 

imagesc(ax1, PDF_j); colorbar; 
set(ax1, 'xtick', xticks, 'ytick', yticks);
set(ax1, 'xticklabels', f, 'yticklabels', orientation);
xtickangle(ax1, 45);
title('Measured Joint Distribution');
xlabel('Frequency (cycles/pixels)');
ylabel('Orientation (Degrees)');

ax2 = subplot(1,2,2); 
imagesc(ax2, PDF_s); colorbar;
set(ax2, 'xtick', xticks, 'ytick', yticks);
set(ax2, 'xticklabels', f, 'yticklabels', orientation); 
xtickangle(ax2, 45);
title('Separable Joint Distribution')
xlabel('Frequency (cycles/pixels)');
ylabel('Orientation (Degrees)');
% ---- END FIRST FIGURE 

% ---- SECOND FIGURE: Individual Distributions side by side.
fig_ind = figure; 

set(fig_ind, 'Units','Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig_ind, 'name', 'Individual Distributions');

% First axis contains the PDF over orientation 
ax3 = subplot(1,2,1);
plot(ax3, orientation, PDF_t, '*');
axis([-inf inf 0 1]);
xlabel('Orientation (Degrees)');
ylabel('Probability');
set(ax3, 'xtick', [0.0, 45.0, 90.0, 135.0], 'ytick', [0.0, 0.5, 1.0]);
title('PDF for Orientation');

% Second axis contains PDF over spatial frequency. 
% if logScale is true, plot with log-log scale. 
ax4 = subplot(1,2,2);
plot(ax4, f, PDF_f, '*');
if logScale
    set(ax4, 'XScale', 'log', 'YScale', 'log'); 
end
axis([-inf inf 0 1]);
set(ax4, 'xtick', f, 'ytick', [0.0, 0.5, 1.0]);
title('PDF for Frequency');
xtickangle(ax4, 45);
xlabel('Frequency (cycles/pixels)');
ylabel('Probability');
% ---- END SECOCND FIGURE 

end


