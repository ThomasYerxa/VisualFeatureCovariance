% speed v. depth calculations

% load file directories containing range data
imageNames = dir('rRange*.*');
nImages    = size(imageNames, 1); 


index = 1; 
for i = 1:nImages
    % load range map 
    currentName  = imageNames(i).name;
    currentImage = load(currentName);
    currentRange = currentImage.rangeMap;
    
    % randomly generate displacement vector in x-y plane
    mag   = 1.0; 
    theta = rand * 2 * pi; 
    d     = mag .* [cos(theta), sin(theta)]; 
    
    % sweep pixel values
    for j = 1:size(currentRange, 1)
        for k = 1:size(currentRange, 2)
             % if pixel value is -1, no range data is available
            if currentRange(j, k, 1) == -1 || currentRange(j, k, 3) == -1 || currentRange(j, k, 3) == -1
                continue
            end
            % calculate angle change that results from displacement
            r_sq     = currentRange(j, k, 1)^2 + currentRange(j, k, 3)^2 + currentRange(j, k, 3)^2;
            cosTheta = r_sq - (currentRange(j, k, 1) * d(1) + currentRange(j, k, 2) * d(2));
            cosTheta = cosTheta / (sqrt(r_sq) * mag); 
            theta    = acos(cosTheta); 
            
            depth(index) = currentRange(j, k, 3); 
            angle(index) = theta; 
              
            index = index + 1;    
        end
    end
end

% Kernel Density Estimation to generate probability distributions
joint = [depth(:), angle(:)];

[f_depth, xi_depth] = ksdensity(depth);
[f_angle, xi_angle] = ksdensity(angle);
[f_joint, xi_joint] = ksdensity(joint);
ksdensity(depth);
ksdensity(angle);
ksdensity(joint);
