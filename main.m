% Run to measure the joint distribution. The working directory should
% contain files that are appropriate for the loadCode set in loadImages.m
clear all; close all;

tic;

plotting = false;
use_gpu  = false;

% Load images
I = loadImages();
nImages = size(I, 3);


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

% choose frequencies to be linearlly or logarithmically spaced
% choose n_f to be a low value for debugging. 
logSpacing = false;
if logSpacing
    filterSize = 80;
    n_waves    = 8;
    f          = logspace(1/filterSize,0.5, 8);
else
    filterSize = 100;
    n_waves    = 8;
    f          = linspace(3/filterSize,0.5, 8);

end

% wavelength is inverse of spatial frequency. 
lambda = 1./f; 

% Orientations in degrees. Includes a regular sampling of orientations
% and a debug option. 
orientation_f     = [0.0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5] * 2 * pi /360;
orientation_debug = [0.0,45.0, 90.0] * 2 * pi /360 ;

% choose orientation option 
orientation = orientation_f;
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
    
    index = 0; 
    for x = 1:size(G_re,3)
        for y = 1:size(G_re, 4)
            index = index + 1; 
            figure(fig_re); hold on;
            subplot(8,8,index); hold on;
            title(['F ' num2str(index)]);
            imagesc(G_re(:,:,x, y));
            axis image off;
            colorbar;

            figure(fig_im); hold on;
            subplot(8,8,index); hold on;
            title(['F ' num2str(index)]);
            imagesc(G_im(:,:,x, y));
            axis image off;
            colorbar;
        
        end    
    end
end

% make complex filter. clear i to make sure you are using sqrt(-1). 
clear i;
G = G_re + 1i * G_im;

% put G and I on GPU
if use_gpu 
    G_gpu = gpuArray(G);
    I_gpu = gpuArray(I);
end 

% Apply filter bank to each image, summing the response values for each
% image.

[xSize, ySize]  = size(conv2(I(:, :, 1), G(:, :, 1, 1),'valid'));

index = 0;
% When working with a large number of images, we may not be able to store
% evevry convolution in one very large matrix (too much memory is used).
% In this case we partition the PDF calculation into batches of 10 images.

% choose whether or not to partition. 
partition = false;

if partition == false
    mag = zeros([xSize, ySize, n_waves, n_theta, nImages]);
end

for j=1:nImages
    for k = 1:n_theta
        for l = 1:n_waves
            index = index + 1;
            if plotting
                figure; hold on;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
                set(gcf, 'name', ['F ' num2str(index)])
            end
            
            % if partition is true the size of mag nevery gros above
            % [xSize , ySize , n_f , n_theta , 10]
            if partition
                % get responses of filter (l, k) to image j. Save to mag
               if use_gpu
                    conv_tmp                       = abs(conv2(I_gpu(:, :, j), G_gpu(:, :, l, k),'valid'));
                    mag(:, :, l, k, mod(j-1,10)+1) = gather(conv_tmp);
                    
                    % clear mag from GPU Memory, reset I and G
                    gpuDevice(1);
                    I_gpu = gpuArray(I);
                    G_gpu = gpuArray(G); 
                    
               else
                    conv_tmp                       = abs(conv2(I(:, :, j), G(:, :, l, k),'valid'));
                    mag(:, :, l, k, mod(j-1,10)+1) = conv_tmp;
               end
                  
            else
                % note the last index of mag is not bounded modularly. 
                % get responses of filter (l, k) to image j. Save to mag
                conv_tmp           = abs(conv2(I(:, :, j), G(:, :, l, k),'valid'));
                mag(:, :, l, k, j) = conv_tmp;

                if plotting
                    subplot(5,5,j); hold on;
                    imagesc(conv_tmp);
                    axis image off;
                    colorbar;
                end
            end 
        end
    end
    
    % check to see whether this is an iteration at which we should
    % calculate an intermediate PDF. Happens at j%10=0 and last image.
    if partition
         
        % adjust weight if last partition is not length 10
        if mod(j, 10) == 0
            [j_temp, t_temp, f_temp, s_temp] = calcPDF(mag, f, orientation);
            weight = 10;
            PDF_j = PDF_j + j_temp*weight; 
            PDF_t = PDF_j + t_temp*weight; 
            PDF_f = PDF_j + f_temp*weight; 
            PDF_s = PDF_j + s_temp*weight; 
            clear mag; 
        end 
        if j == nImages
            % the last intermediate PDF may have been calculalted using <
            % 10 images. Reset its weight in the sum accordingly. 
            [j_temp, t_temp, f_temp, s_temp] = calcPDF(mag, f, orientation);
            weight = mod(nImages, 10);
            PDF_j = PDF_j + j_temp*weight; 
            PDF_t = PDF_j + t_temp*weight; 
            PDF_f = PDF_j + f_temp*weight; 
            PDF_s = PDF_j + s_temp*weight; 
        end
    end    
end

% Calculate and visualize joint PDF from filter responses

% if partition is on, we take the final PDF to be the average of all 
% of the intermediate PDF values that have already been summed. 
if partition
    PDF_j = PDF_j/nImages; 
    PDF_t = PDF_t/nImages; 
    PDF_f = PDF_f/nImages; 
    PDF_s = PDF_s/nImages; 
    
else
    % else we calculate the total pdf using every filter response. 
    [PDF_j, PDF_t, PDF_f, PDF_s] = calcPDF(mag, f, orientation); 
end


vizPDF(PDF_j, PDF_t, PDF_f, PDF_s, f,orientation);

toc;