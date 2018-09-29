% Run to measure the joint distribution. The working directory should 
% contain .TIF files from the McGill Colour Image Database. 

% Load images
I = loadImages();
sz = size(I);
nImages = sz(3);

lambda = [10, 20]; % Set of wavelengths (in pixels)
n_waves = size(lambda);
n_waves = n_waves(2);
one = ones(size(lambda));
f = one ./ lambda; % Spatial frequency is inverse of wavelength

% Orientations in degrees
orientation = [0, 45] * 2 * pi /360 ; 
n_theta = size(orientation);
n_theta = n_theta(2);

% Generate filter bank 
bw = 1.0;
ar = 0.5; 
psi = 0.0;

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


G = G_re + 1i * G_im; 

% Apply filter bank to each image, summing the response values for each 
% image.
mag = conv2(I(:, :, 1),G(:,:, 1));
mag = abs(mag);
nFilters = n_waves*n_theta;

% Because the preliminary implementation calculates PDFs by summing 
% response intensities, the response to each image does not need to be
% stored separately. This will probably need to be fixed in the future. 
for k = 1:nFilters
    mag(:,:,k) = conv2(I(:, :, 1), G(:, :, k));
    for l = 2:nImages
         current = conv2(I(:, :, l), G(:, :, k));
         current = abs(current);
         mag(:, :, k) = mag(:, :, k) + current; 
         (k-1)*nImages + l; 
    end   
end

% Calculate and visualize joint PDF from filter responses 
PDF = calcPDF(abs(mag), f, orientation); 
% Visualize PDF has not been implemented yet. 
%vizPDF(PDF); 