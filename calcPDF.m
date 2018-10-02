%-------------------------------------------------------------------------
% Calculates PDF(f, theta), PDF(theta), and PDF(f), by summing the 
% intensity values of all responses, and normalizing by total intensity. 
% Another implementation could use a 'winner take all,' approach 
% by keeping track of which filter causes the maximum response for each 
% image patch, then normalize by the number of image patches. 

% Parameters: 
% -----------
%  responses : The set of filter responses (summed over images)
%              size: (sz_conv1, sz_conv2, nFilters, nImages)

%  frequencies : the list of spatial frequencies used. 
%                size: (1, n_frequencies)

%  orientations : the list of orientations used. 
%                 size: (1, n_theta)
%
% Returns: 
% ---------
% joint_PDF : elements sum to one. size: (n_theta, n_frequencies)
%
% Currently set to return the joint distribution, but prints 
% each distribution. 
%-------------------------------------------------------------------------

function [PDF] = calcPDF(responses, frequencies, orientations)

% initialize all counts to zero
n_theta = size(orientations);
n_theta = n_theta(2);
n_f = size(frequencies);
n_f = n_f(2);
prob_joint = zeros(n_theta,n_f);
prob_f = zeros(1,n_f);
prob_theta = zeros(1, n_theta);
total = 0.0;
winnerTakeAll = true; 
[xSize, ySize, nFilters, nImages] = size(responses);

if winnerTakeAll
    % for each filter response (one per (x,y) coordinate per image)
    % determine which filter caused the maximal response. 
    % Increment the winner's bin. 
    for n = 1:nImages
        for x = 1:xSize
            for y = 1:ySize
                r = squeeze(responses(x, y, :, n));
                [~, max_filter_ind] = max(r); 
                [f_ind, theta_ind] = ind2sub([n_f, n_theta], max_filter_ind); 
                prob_joint(theta_ind, f_ind) = prob_joint(theta_ind, f_ind) + 1;
                prob_f(f_ind) = prob_f(f_ind) + 1; 
                prob_theta(theta_ind) = prob_theta(theta_ind) + 1;
                total = total + 1.0;
            end
        end
    end
            
else
    % sum the intensity accrued by each filter, and save into approriate 
    % variables. 
    for j = 1:nImages
        for k = 1:n_theta
            for l = 1:n_f
                index = (k-1)*n_f + l;
                filterScore = sum(responses(:,:,index, j), 'all');
                prob_joint(k, l) = prob_joint(k, l) + filterScore; 
                prob_f(l) = prob_f(l) + filterScore; 
                prob_theta(k) = prob_theta(k) + filterScore;  
                total = total + filterScore; 
            end
        end
    end
end
    
    
% normalize the distributions (and output results for now)
PDF_joint = prob_joint / total
PDF_theta = prob_theta / total
PDF_f = prob_f / total 

% if the two distributions are separale PDF_joint = PDF_sep
PDF_sep = (PDF_theta' * PDF_f)

[PDF] = PDF_joint;

end