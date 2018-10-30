% This function defines a 1-D gaussian tuning curve

function [func] = gaussian1D(mean, std, scale)

func = @(x) exp(-(x - mean).^2 ./ (2 * std^2)) ./ scale;

end