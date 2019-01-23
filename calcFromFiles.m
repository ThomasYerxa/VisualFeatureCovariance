% calculate and visualize PDF from saved files.
% assumes that in the working directory there are files specifying: 
% 1) partition_num # of partitions of mag data
% 2) files of mag data in the form "mag_#.mat" where # is a number between 1
%    and partition_num. 
% 3) the frequency and orientation settings used to generate the mag data


% load settings data
frequency    = load('f');
theta        = load('orientation');
nTheta       = size(theta ,2);
nF           = size(frequency,2);
 
% load mag data. 
magNames = dir('*.mat');
nMags    = length(magNames); 

% set range of image data to use, default is to use all
range1 = 1;
range2 = nMags;

% option to take slice of PDF (i.e. view 1D PDF of orientation given 
% specific frequency value)
sliceTheta = false; 
sliceF     = false; 
sliceIndex = 1; 
if sliceTheta
    nF = 1;
end
if sliceF
    nTheta = 1;
end

% initialize distributions
probJoint  = zeros(nTheta,nF);
probF      = zeros(1,nF);
probTheta  = zeros(1, nTheta);
total      = 0.0; % tracker used for normalization

sizeRef = load('mag_1');

for i = range1:range2
    curMag = load('mag_' + string(i)); 
    
    % check that size matches refernce size
    a = size(rgb2lms(curMag));
    b = size(LMS(:, :, :, 1));
    if a(1) ~= b(1)
        continue
    end
    
    % take slice if appropriate
    if sliceTheta
        curMag = curMag(:, :, sliceIndex, :);
    end
    if sliceF
        curMag = curMag(:, :, :, sliceIndex);
    end
    
    % Assign response values. If slice is taken then probJoint = probSliced
    % and probNotSliced will be = 1.0
    for k = 1:n_theta
        for l = 1:n_f 
            filterScore     = sum(curMag(:,:, l, k), 'all');
            probJoint(k, l) = probJoint(k, l) + filterScore; 
            probF(l)        = probF(l) + filterScore; 
            probTheta(k)    = probTheta(k) + filterScore;  
            total           = total + filterScore; 
        end
    end  
end

% normalize the distributions
PDF_joint = prob_joint / total;
PDF_theta = prob_theta / total;
PDF_f     = prob_f / total; 

