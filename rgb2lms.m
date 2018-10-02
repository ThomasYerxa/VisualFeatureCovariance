%-----------------------------------------------------------------------
% This program computes the LMS cone responses of the RGB images in the 
% McGill Colour Image Database.
% The LMS cone responses are based on the Smith and Pokorney functions
%
% Usage:          
%       [LMS] = rgb2lms(I)
% Input:
%       I  = Image name
% Output:
%       LMS = LMS image matrix.  LMS is a three dimensional matrix (type double )
%       LMS(:,:,1) = L cone responses
%       LMS(:,:,2) = M cone responses
%       LMS(:,:,3) = S cone responses
%
% Author: 
%       Adriana Olmos Mar/2004 - Mcgill Vision Research
%
% Last modification:
%       Jul/2005 - Faster computation of conversion (Thorsten Hansen, PhD)
%       Feb/2005 - Exposure time stamp added
%       Mar/2004 - Main code
%-----------------------------------------------------------------------

function [LMS] = rgb2lmsnew(imageName)

%Checking the input arguments
if nargin ~= 1
     error('Wrong number of input arguments')
end
%Error message
Message = ['The name of the camera used to take this image could not be found in the header.  Please make sure that the image is TIF format and downloaded from the McGill Colour Calibration Database.'];

% Checking that the image can be read.
fid = fopen(imageName, 'r');
if (fid == -1)
  error(['The file could not be open, please check the image name, path and/or image permisions.']);
end

% Reading image
I = double(imread(imageName));
% Reading Image header
info = imfinfo(imageName);

% PARAMETRES
if (info(1).ImageDescription(1) == ' ')
 error(Message);
 else
  if(all(info.ImageDescription(1:5)=='merry') && all(info(1).Format(1:3)=='tif'))
  % Camera parameters - Merry
     T = [0.428443253	0.495562896	 0.075993851
          0.243026144	0.614128681	 0.142845175
          0.155766424	0.132343175	 0.711890401];     
   a_R = 7.320565961;
   a_G = 12.0579051;
   a_B = 10.6112984;
     b = 1.008634316;     
 else
   if (all(info(1).ImageDescription(1:5)=='pippi') && all(info(1).Format(1:3)=='tif'))
  % Camera parameters - Pippin
     T = [0.431088433	 0.494438389  0.074473178
          0.245488691	 0.614786761  0.139724548
          0.166472303	 0.124487321  0.709040376];
   a_R = 5.562441185;
   a_G = 8.876002262;   
   a_B = 7.233814813;
     b = 1.009696031;  
     else
      error(Message);
    end
  end
end

%Linearising the Image
II = zeros(size(I));
II(:,:,1) = a_R.*(b.^I(:,:,1)-1);
II(:,:,2) = a_G.*(b.^I(:,:,2)-1);
II(:,:,3) = a_B.*(b.^I(:,:,3)-1);

%Exposure Time Stamp
ET = info(1).ImageDescription(8:length(info(1).ImageDescription)); % Exposure Time
%display(['Exposure time = ' ET 'secs']);
II = II./str2double(ET);

rows=size(II,1); cols=size(II,2);
LMS = zeros(size(II));

%Computing the LMS cone responses through the matrix transformation "T"
% This is much faster (about 30 times) than the two for loops...
% At every image position with column vector [r; g; b]
% compute the matrix transformation T*[r; g; b]
%
% LMS = reshape((T*reshape(II, rows*cols, 3)')', rows, cols, 3);
%
% Equivalently, one can compute [r g b]*T' 
% (which runs faster in this case because instead of two transposition of
% the whole image one needs only one for T)
% Contribution: Thorsten Hansen, PhD at psychol.uni-giessen.de

LMS = reshape(reshape(II, rows*cols, 3)*T', rows, cols, 3);
 
 
