%-------------------------------------------------------------------------
% Function for loading images from the McGill Colour Image Database. 
% Loads all .TIF files in the working directory and returns them as 
% an array of intensity values.
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

function [M] = loadImages()
% load .TIF files into struct
images = dir('*.TIF');
nImages = length(images); 

debug = true;
if debug
    nImages = 2; % reset nImages to 2 to reduce runtime. 
    
% Collect LMS value arrays, extract M channel (which encodes intensity)
firstName = images(1).name; 
[LMS] = rgb2lms(firstName); 
[M] = LMS(:,:,2);

for i=2:nImages
    currentFileName = images(i).name;
    LMS(:, :, :, i) = rgb2lms(currentFileName); 
    M(:, :, i) = LMS(:, :, 2, i);  
end  

end

    
    
    
    






