% generate a filter bank 

function G = filterBank(frequency, orientation, filterSize, plotting)

% filter varialbes
bw  = 1.0;
ar  = 1.0;
psi = 0.0;

lambda = 1./frequency;
n_t    = size(orientation, 2);
n_f    = size(lambda, 2);

G_re = gabor_fn(bw, ar, psi, lambda(1), orientation(1), filterSize);
G_im = gabor_fn(bw, ar, psi + pi/2, lambda(1), orientation(1), filterSize);
% Build G s.t. size(G) = (sz, sz, n_f, n_theta)
for k = 1:n_t
    for l = 1:n_f
        if not(and(k == 1, l == 1))
            G_re(:,:, l, k) = gabor_fn(bw, ar, psi,  lambda(l), orientation(k), filterSize);
            G_im(:,:, l, k) = gabor_fn(bw, ar, psi + pi/2, lambda(l), orientation(k), filterSize);
        end
    end
end

index = 0;
if plotting
    fig_re = figure;
    set(fig_re, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig_re, 'name', 'Filter bank -- real')
    
    fig_im = figure;
    set(fig_im, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig_im, 'name', 'Filter bank -- imag')
    
    index = 0; 
    for x = 1:size(G_re,3)
        for y = 1:size(G_re, 4)
            index = index + 1; 
            figure(fig_re); hold on;
            subplot(8,8,index); hold on;
            title(['F ' num2str(index)]);
            imagesc(G_re(:,:,x, y));
            axis image off;
            colorbar;

            figure(fig_im); hold on;
            subplot(8,8,index); hold on;
            title(['F ' num2str(index)]);
            imagesc(G_im(:,:,x, y));
            axis image off;
            colorbar;
        
        end    
    end
end

G = G_re + G_im;

end

