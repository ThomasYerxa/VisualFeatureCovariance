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
    bad_images = [16, 21, 25, 30, 31, 37, 43, 45, 48, 49, 56, 60, 65, 80, 83, 84, 85, 86, 94 ];
    
    if ismember(i, bad_images)
        continue; 
    end
    
    % randomly generate displacement vector in x-y plane
    mag   = 0.05; 
    theta = rand * 2 * pi; 
    d     = mag .* [cos(theta), sin(theta)]; 
    
    % sweep pixel values
    for j = 1:size(currentRange, 1)
        for k = 1:size(currentRange, 2)
             % if pixel value is -1, no range data is available
            if currentRange(j, k, 1) == -1 || currentRange(j, k, 2) == -1 || currentRange(j, k, 3) == -1
                continue
            end
            % calculate angle change that results from displacement
            x  = currentRange(j, k, 1); y = currentRange(j, k, 2); z = currentRange(j, k, 3); 
            
            r_sq     = x^2 + y^2 + z^2; 
            cosTheta = x*(x - d(1)) + y*(y-d(2)) + z^2;
            mag_rd   = (x - d(1))^2 + (y - d(2))^2 + z^2;
            cosTheta = cosTheta / (sqrt(r_sq) * sqrt(mag_rd));
            
            
            
            if cosTheta > 1 

          
                continue; 
           
            elseif cosTheta < -1 
                
                
                continue; 
            end
            theta    = acos(cosTheta); 
            
            
            
            chance = rand;
            if x < -100 || theta > 0.04
                continue; 
            end
                
            
            if chance < 0.0003
                depth(index) = -currentRange(j, k, 1); 
                angle(index) = theta; 
                index = index + 1;    
            end
            
            
        end
    end
end

% Kernel Density Estimation to generate probability distributions
joint = [depth(:), angle(:)];

[f_depth, xi_depth, bwd] = ksdensity(depth);
[f_angle, xi_angle, bwa] = ksdensity(angle);
[f_joint, xi_joint, bwj] = ksdensity(joint);


figure;
ksdensity(joint, xi_joint)

figure;
plot(xi_depth, f_depth);



figure; 
plot(xi_angle, f_angle); 

depth_j = unique(xi_joint(:, 1));
angle_j = unique(xi_joint(:, 2)); 
PDF_j   = zeros(30);
for i = 1:30
    for j = 1:30
        for k = 1:900
            if xi_joint(k,1)==depth_j(i) && xi_joint(k,2)==angle_j(j)
                PDF_j(i,j) = f_joint(k);
            end
        end
    end
end
save('depth', 'depth_j');
save('angle', 'angle_j');
save('PDF_jda', 'PDF_j');
        


