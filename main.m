% Run to measure the joint distribution. The working directory should
% contain .TIF files from the McGill Colour Image Database.
clear all; close all;

tic;

plotting = false;

% Load images
I = loadImages();
sz = size(I);
nImages = sz(3);

if plotting
    
    figure; hold on;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(gcf, 'name', 'Images')
    
    for x = 1:nImages
        subplot(5,5,x); hold on;
        imagesc(I(:,:,x));
        axis image off;
        colormap gray; colorbar
    end
end

lambda_n     = [40, 39, 38, 37, 35, 30, 20, 5]; % Set of wavelengths (in pixels)
lambda_hf    = [30, 20, 10, 5, 4, 3, 2, 1];
lambda_debug = [40, 38];

lambda      = lambda_n; 
n_waves     = size(lambda, 2);
one         = ones(size(lambda));
f           = one ./ lambda; % Spatial frequency is inverse of wavelength

% Orientations in degrees
orientation_f = [0.0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5] * 2 * pi /360;
orientation_debug = [0.0, 22.5] * 2 * pi /360 ;

orientation = orientation_debug;
n_theta     = size(orientation, 2);

% Generate filter bank
nFilters    = n_waves*n_theta;
bw          = 1.0;
ar          = 1.0;
psi         = 0.0;

%
G_re = gabor_fn(bw, ar, psi, lambda(1), orientation(1));
G_im = gabor_fn(bw, ar, psi + pi/2, lambda(1), orientation(1));

% Build G s.t. size(G) = (sz, sz, n_f, n_theta)
for k = 1:n_theta
    for l = 1:n_waves
        if not(and(k == 1, l == 1))
            G_re(:,:, l, k) = gabor_fn(bw, ar, psi,  lambda(l), orientation(k));
            G_im(:,:, l, k) = gabor_fn(bw, ar, psi + pi/2, lambda(l), orientation(k));
        end
    end
end

% draw filter bank
if plotting
    fig_re = figure;
    set(fig_re, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig_re, 'name', 'Filter bank -- real')
    
    fig_im = figure;
    set(fig_im, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig_im, 'name', 'Filter bank -- imag')
    % NEED TO FIX FOR NEW FORMAT OF G
    for x = 1:size(G_re,3)
        figure(fig_re); hold on;
        subplot(8,8,x); hold on;
        title(['F ' num2str(x)]);
        imagesc(G_re(:,:,x));
        axis image off;
        colorbar;
        
        figure(fig_im); hold on;
        subplot(8,8,x); hold on;
        title(['F ' num2str(x)]);
        imagesc(G_im(:,:,x));
        axis image off;
        colorbar;
    end
end

% make complex filter
G = G_re + 1i * G_im;

% Apply filter bank to each image, summing the response values for each
% image.

[xSize, ySize]  = size(conv2(I(:, :, 1), G(:, :, 1, 1),'valid'));
mag             = zeros([xSize, ySize, n_waves, n_theta, nImages]);

for k = 1:n_theta
    for l = 1:n_waves
        
        if plotting
            figure; hold on;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
            set(gcf, 'name', ['F ' num2str(k)])
        end
    
        for j = 1:nImages
            % get responses of filter (l, k) to image j. Save to mag
            conv_tmp           = abs(conv2(I(:, :, j), G(:, :, l, k),'valid'));
            mag(:, :, l, k, j) = conv_tmp;
            % NEED TO FIX FOR NEW FORMAT OF G
            if plotting
                subplot(5,5,j); hold on;
                imagesc(conv_tmp);
                axis image off;
                colorbar;
            end
        end
    end
       
end

% Calculate and visualize joint PDF from filter responses
[PDF_j, PDF_t, PDF_f, PDF_s] = calcPDF(mag, f, orientation);
% Visualize PDF has not been implemented yet.
vizPDF(PDF_j, PDF_t, PDF_f, PDF_s, f,orientation);

toc;