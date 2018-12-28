% Run to measure the joint distribution. The working directory should
% contain files that are appropriate for the loadCode set in loadImages.m
clear all; close all;

tic;

plotting = false;
use_gpu  = false;
sv_vars  = false; 

% Load images
I = loadImages(1);
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

% Orientations in degrees. Includes a regular sampling of orientations
% and a debug option. 
orientation_f     = [0.0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5] * 2 * pi /360;
orientation_debug = [0.0,45.0, 90.0] * 2 * pi /360 ;

% choose orientation option 
orientation = orientation_f;
n_theta     = size(orientation, 2);

% Generate filter bank
G           = filterBank(f, orientation_f, filterSize, plotting); 

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
partition     = true;
partition_num = 0;

if partition == false
    mag = zeros([xSize, ySize, n_waves, n_theta, nImages]);
end

for j=1:nImages
    tic
    
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
        
        if mod(j, 10) == 0 
            mag = single(mag);
            partition_num = partition_num + 1;
            save('mag'+string(partition_num), 'mag', '-v7.3');
            clear('mag');
        end
        
        if j == nImages
            mag = single(mag);
            partition_num = partition_num + 1;
            save('mag'+string(partition_num), 'mag', '-v7.3');
            clear('mag');
        end
        
    end
    
    toc;
    
    % load more images if needed and possible. 
    if j == nImages
       if j >= 100
           continue;
       else
           I       = loadImages(j);
           nImages = nImages + size(I, 3);
       end   
    end
   
end

% Calculate and visualize joint PDF from filter responses

% if partition is on, we take the final PDF to be the average of all 
% of the intermediate PDF values that have already been summed. 
if partition
    % save the variables needed to run calcPDF on the saved "mag#.mat" 
    % files
    save('theta', 'orientation');
    save('f', 'f');
    save('partition_num', 'partition_num');
      
else  
    % else we calculate the total pdf using every filter response. 
    [PDF_j, PDF_t, PDF_f, PDF_s] = calcPDF(mag, f, orientation); 
    vizPDF(PDF_j, PDF_t, PDF_f, PDF_s, f,orientation);

end



toc;