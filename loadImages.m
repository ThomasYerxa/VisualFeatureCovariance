%-------------------------------------------------------------------------
% Function for loading images.
%
% Functionality depends on the value of loadCode and the working directory:
%
% loadCode = 0: load images from the McGill Colour Image Database. 
% Loads all .TIF files in the working directory and returns them as 
% an array of intensity values. when run the working directory must
% contain .TIF files. 
%
% loadCode = 1: load Noise images from NS1,2,3 files. 
% Loads all .mat files in the working directory and returns them as 
% an array of intensity values. Working directory must be NS1, NS2, or NS3
%
% loadCode = 2: load images from the Van Hateren Database. 
% Loads all .iml files in the working directory and returns them as 
% an array of intensity values. Working directory must contain .iml files.
%
% 
% Parameters: 
% -----------
%  NONE
%
% Returns: 
% ---------
% M : array of images. size: (image_x, image_y, nImages) where nImages
%     is either the number of images in the working directory, or 2 if 
%     debug is set to true.
%-------------------------------------------------------------------------

function [M] = loadImages(startIndex)


loadCode = 0;
debug    = false; 

if loadCode == 0
   % load .TIF files into struct
    images = dir('*.TIF');
    nImages = length(images); 
    if debug
        nImages = 10; % reset nImages to 2 to reduce runtime. 
    end
    
    % start index must be below the number of images
    if startIndex > nImages
        return;
    end
    
    % Collect LMS value arrays, extract M channel (which encodes intensity)
    firstName = images(startIndex).name; 
    [LMS] = rgb2lms(firstName); 
    [M] = LMS(:,:,2);
    
    for i=startIndex+1:nImages
        % only return 50 images at a time. 
        if i > startIndex+50
            continue;
        end
        % some files are mishaped. Skip them. 
        currentFileName = images(i).name;
        a = size(rgb2lms(currentFileName));
        b = size(LMS(:, :, :, 1));
        if a(1) ~= b(1)
            continue
        end
        if sum(a(:,:,2), 'all') <= 1
            continue
        end
        
        
        LMS(:, :, :, i) = rgb2lms(currentFileName); 
        M(:, :, i) = LMS(:, :, 2, i);  
        fclose('all');
        
   
    end  
    
    return;
end

if loadCode == 1
    nImages = 10; 
    if debug
        nImages = 2;
    end
    for i=1:nImages
       if i < 10
           temp = load(['noiseim_00' num2str(i)]);
       else 
           temp = load('noiseim_010');
       end
       M(:, :, i) = temp.im; 
    end
    return 
end

if loadCode == 2
   w = 1536; h = 1024;
   % load all .iml files in current directory
   images = dir('*.iml');
   nImages = length(images); 
   if debug
       nImages = 2; % reset nImages to 2 to reduce runtime. 
   end
   firstName = images(1).name; 
   f = fopen(firstName, 'rb','ieee-be');
   temp = fread(f, [w, h], 'uint16');
   M(:, :, 1) = temp(4:1532, 4:1020); 
   
   for i=2:nImages
      name = images(i).name;
      f = fopen(name, 'rb', 'ieee-be');
      temp = fread(f, [w, h], 'uint16'); 
      M(:, :, i) = temp(4:1532, 4:1020);  
   end
   
   return;
end


end

    
    
    
    






