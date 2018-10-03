% Run to measure the joint distribution. The working directory should
% contain .TIF files from the McGill Colour Image Database.
clear all; close all;

tic;

plotting = 1;

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

lambda      = [5, 10, 20]; % Set of wavelengths (in pixels)
n_waves     = size(lambda);
n_waves     = n_waves(2);
one         = ones(size(lambda));
f           = one ./ lambda; % Spatial frequency is inverse of wavelength

% Orientations in degrees
orientation = [0, 45, 90, 135] * 2 * pi /360 ;
n_theta = size(orientation);
n_theta = n_theta(2);

% Generate filter bank
nFilters    = n_waves*n_theta;
bw          = 1.0;
ar          = 1.0;
psi         = 0.0;

%
G_re = gabor_fn(bw, ar, psi, lambda(1), orientation(1));
G_im = gabor_fn(bw, ar, psi + pi/2, lambda(1), orientation(1));

for k = 1:n_theta
    for l = 1:n_waves
        if not(and(k == 1, l == 1))
            G_re(:,:, (k-1)*n_waves + l) = gabor_fn(bw, ar, psi,  lambda(l), orientation(k));
            G_im(:,:, (k-1)*n_waves + l) = gabor_fn(bw, ar, psi + pi/2, lambda(l), orientation(k));
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
    
    for x = 1:size(G_re,3)
        figure(fig_re); hold on;
        subplot(5,5,x); hold on;
        title(['F ' num2str(x)]);
        imagesc(G_re(:,:,x));
        axis image off;
        colorbar;
        
        figure(fig_im); hold on;
        subplot(5,5,x); hold on;
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

[xSize, ySize]  = size(conv2(I(:, :, 1), G(:, :, 1),'valid'));
mag             = zeros([xSize, ySize, nFilters, nImages]);

for k = 1:nFilters
    
    if plotting
        figure; hold on;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
        set(gcf, 'name', ['F ' num2str(k)])
    end
    
    for l = 1:nImages
        
        conv_tmp        = abs(conv2(I(:, :, l), G(:, :, k),'valid'));
        mag(:, :, k, l) = conv_tmp;
        
        if plotting
            subplot(5,5,l); hold on;
            imagesc(conv_tmp);
            axis image off;
            colorbar;
        end
    end
end

% Calculate and visualize joint PDF from filter responses
PDF = calcPDF(mag, f, orientation);
% Visualize PDF has not been implemented yet.
%vizPDF(PDF);

toc;