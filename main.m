% Run to measure the joint distribution. The working directory should
% contain files that are appropriate for the loadCode set in loadImages.m
clear all; close all;

tic;

plotting = false;
use_gpu  = false;
sv_vars  = false; 

% Load image directories
loadCode = 0; 

% loadCode = 0 --> loading from McGill Tabby
if loadCode == 0
    % load .TIF files into struct
    images = dir('*.TIF');
end

% loadCode = 1 --> loading from van Hateren
if loadCode == 1
    % load .iml files into struct
    images = dir('*.iml'); 
end
nImages = length(images);

% reduce nImages for debugging
debug = true; 
if debug
    nImages = 2; 
end

% NEED TO FIND BETTER WAY TO PLOT IMAGES
%{ 
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
%}
    
filterSize = 100;
n_waves    = 16;
f          = linspace(3/filterSize,0.5, n_waves);

% Orientations in degrees. Includes a regular sampling of orientations
% and a debug option. 
orientation_f     = linspace(0, 180, 16) * 2 * pi /360;
orientation_f     = orientation_f(1:end-1); 

orientation_debug = [0.0,45.0, 90.0] * 2 * pi /360 ;

% choose orientation option 
orientation = orientation_f;
n_theta     = size(orientation, 2);

% Generate filter bank
G           = filterBank(f, orientation_f, filterSize, plotting); 

% Apply filter bank to each image
index = 0;
for j=1:nImages
    tic
    currentName = images(j).name; 
    if loadCode == 0
        % Load image, convert to intensity values
        [LMS] = rgb2lms(currentName); 
        [I] = LMS(:,:,2);
    end
    if loadCode == 1
        w = 1536; h = 1024;
        f = fopen(currentName, 'rb','ieee-be');
        temp = fread(f, [w, h], 'uint16');
        [I] = temp(4:1532, 4:1020);
    end
    
    for k = 1:n_theta
        for l = 1:n_waves
            % Plot the current image
            if plotting
                figure; hold on;
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
                set(gcf, 'name', ['F ' num2str(index)])
            end
            
            % Calculate response of image J to filter LK
            mag(:, :, l, k) = abs(conv2(I, G(:, :, l, k),'valid'));
            mag             = single(mag);
            index           = index + 1;
            
        end
    end
    if loadCode == 0
        save('mag_mcgill_'+string(j), 'mag', '-v7.3');
    end
     if loadCode == 1
        save('mag_vanhateren_'+string(j), 'mag', '-v7.3');
    end
    clear('mag');
    toc
end   

save('f', 'f');
save('orientation','orientation');

toc;