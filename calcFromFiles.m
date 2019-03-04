% calculate and visualize PDF from saved files.
% assumes that in the working directory there are files specifying: 
% 1) partition_num # of partitions of mag data
% 2) files of mag data in the form "mag_#.mat" where # is a number between 1
%    and partition_num. 
% 3) the frequency and orientation settings used to generate the mag data

winnerTakeAll = false; 
% load settings data
frequency    = load('f.mat');
frequency    = frequency.f;
theta        = load('orientation.mat');
theta        = [theta.orientation pi];
nTheta       = size(theta ,2);
nF           = size(frequency,2);
nFilters     = nF * nTheta; 
 
% load mag data. 
magNames = dir('mag*.mat');
nMags    = length(magNames); 

% set range of image data to use, default is to use all
range1 = 1;
range2 = nMags;

% initialize distributions
probJoint  = zeros(nTheta,nF);
probF      = zeros(1,nF);
probTheta  = zeros(1, nTheta);
total      = 0.0; % tracker used for normalization

% option to take slice of PDF (i.e. view 1D PDF of orientation given 
% specific frequency value)
sliceTheta = false; 
sliceF     = false; 
sliceIndex = 5; 
if sliceTheta
    nF = 1;
end
if sliceF
    nTheta = 2;
end





for i = range1:range2
    curName = magNames(1).name; 
    curMag  = load(curName);
    curMag  = curMag.mag;
    curMag  = curMag.^2; 
    %{
    % check that size matches refernce size
    a = size(rgb2lms(curName));
    b = size(LMS(:, :, :, 1));
    if a(1) ~= b(1)
        continue
    end
    %}
    
    % take slice if appropriate
    if sliceTheta
        curMag = curMag(:, :, sliceIndex, :);
    end
    if sliceF
        curMag = curMag(:, :, :, sliceIndex);
    end
    % Assign response values. If slice is taken then probJoint = probSliced
    % and probNotSliced will be = 1.0
    
    if winnerTakeAll
        xSize = size(curMag, 1); ySize = size(curMag, 2);
        for x = 1:xSize
            for y = 1:ySize
                
                % reshape unrolls by traversing down columns
                % [[a,b];[c,d]]-->[a, c, b, d]. 
                % This matches the convention of ind2sub, so we recover 
                % the frequency/orientation index of the filter that caused
                % the maximal response at position (x, y) in image n. 
                r                   = squeeze(curMag(x, y, :, :));
                r                   = reshape(r, [1, nFilters]);
                [~, max_filter_ind] = max(r); 
                [f_ind, theta_ind]  = ind2sub([nF, nTheta], max_filter_ind); 
                
                probJoint(theta_ind, f_ind) = probJoint(theta_ind, f_ind) + 1;
                probF(f_ind)                = probF(f_ind) + 1; 
                probTheta(theta_ind)        = probTheta(theta_ind) + 1;
                total                        = total + 1.0;
            end
        end
    else
        for k = 1:nTheta - 1
            for l = 1:nF
                filterScore     = sum(curMag(:,:, l, k), 'all');
                probJoint(k, l) = probJoint(k, l) + filterScore; 
                probF(l)        = probF(l) + filterScore; 
                probTheta(k)    = probTheta(k) + filterScore;  
                total           = total + filterScore; 
                
                if k == 1
                    probJoint(nTheta, l) = probJoint(k, l) + filterScore; 
                    probF(l)             = probF(l) + filterScore; 
                    probTheta(nTheta)    = probTheta(k) + filterScore;  
                    total                = total + filterScore; 
                end
                
                
            end
        end     
    end
    
    
end

% check sizes
size(theta)

% normalize the distributions
PDF_joint = probJoint / total;
PDF_theta = probTheta / total;
PDF_f     = probF / total;
PDF_s     = PDF_theta' * PDF_f;

% save distributions
%{
save('PDF_joint', 'PDF_joint');
save('PDF_theta', 'PDF_theta');
save('PDF_f', 'PDF_f');
save('PDF_s', 'PDF_s');
%}

save('PDF_joint', 'PDF_joint');
save('PDF_theta', 'PDF_theta');
save('PDF_f', 'PDF_f');
save('PDF_s', 'PDF_s');
save('theta', 'theta'); 
save('frequency','frequency');

vizPDF(PDF_joint, PDF_theta, PDF_f, PDF_s, frequency, theta);