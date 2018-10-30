% A 2 dimensional gaussian function with Sigma = I. 
%
% This could be generalized to include the case of non-uniform variance and
% non-zero covariance, but for the purpose of this work we are interested 
% in how the displacement of a domain changes the shapes of evenly spaced 
% identical tuning curves, not the shapes of functionally different 
% tuning curves. 
%
% The ability to scale the function is included for convenience

function [func] = gaussian2D(mu_x, mu_y, scale)
func = @TwoDGaussian;

    function y = TwoDGaussian(x, y)
        y = exp(-0.5 * ((x - mu_x).^2 + (y-mu_y).^2)) ./ scale;
    end

end