% test PDF calculations on synthetic image


plotting = true; 
filterSize = 80;
% largest possible frequency range given filterSize
f      = linspace(1/filterSize,0.5, 8);
lambda = 1./f; 
theta  = [0.0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5] * 2 * pi /360;

n_f = size(f, 2);
n_t = size(theta, 2);

% desgin test stimuli
stim_harmonic = [5, 10];
stim_angle    = [22.5, 112.5]; 
stim1         = mkSine(800, stim_harmonic(1));
stim1         = imrotate(stim1, stim_angle(1), 'bilinear', 'crop');
stim1         = stim1(200:600, 200:600);

stim2         = mkSine(800, stim_harmonic(2));
stim2         = imrotate(stim2, stim_angle(2), 'bilinear', 'crop');
stim2         = stim2(200:600, 200:600);

% plot stimulus 
if plotting
    figure; hold on; 
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(gcf, 'name', 'Stimulus 1');
    imagesc(stim1);
    axis image off;
    colormap gray; colorbar
    figure; hold on; 
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(gcf, 'name', 'Stimulus 2');
    imagesc(stim2);
    axis image off;
    colormap gray; colorbar
end

% filter varialbes
bw          = 1.0;
ar          = 1.0;
psi         = 0.0;


% set figure variables if plotting.
if plotting
    fig_re = figure;
    set(fig_re, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig_re, 'name', 'Filter bank -- real')
    
    fig_im = figure;
    set(fig_im, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(fig_im, 'name', 'Filter bank -- imag')
    
    index = 0; 
    
end

% intitialize all responses to zero
response_f1 = zeros(1, n_f);
response_t1 = zeros(1, n_t);
response_f2 = zeros(1, n_f);
response_t2 = zeros(1, n_t);

for j = 1:n_f
    for k = 1:n_t
        G_re = gabor_fn(bw, ar, psi,  lambda(j), orientation(k));
        G_im = gabor_fn(bw, ar, psi + pi/2,  lambda(j), orientation(k));
        % visualize filters
        if plotting
            index = index + 1; 
            figure(fig_re); hold on;
            subplot(8,8,index); hold on;
            title(['F ' num2str(index)]);
            imagesc(G_re);
            axis image off;
            colorbar;

            figure(fig_im); hold on;
            subplot(8,8,index); hold on;
            title(['F ' num2str(index)]);
            imagesc(G_im);
            axis image off;
            colorbar;
        end
        % calculate response of filter
        r1             = abs(conv2(stim1, (G_re + 1i*G_im), 'valid'));
        r1             = sum(r1, 'all');
        response_f1(j) = response_f1(j) + r1;
        response_t1(k) = response_t1(k) + r1;
        
        r2             = abs(conv2(stim2, (G_re + 1i*G_im), 'valid'));
        r2             = sum(r2, 'all');
        response_f2(j) = response_f2(j) + r2;
        response_t2(k) = response_t2(k) + r2;
    end
end

fig_responses = figure; 
set(fig_responses, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig_responses, 'name', 'Filter Responses Stimulus One: harmonic = 5, theta = 22.5');
subplot(1, 2, 1); hold on; 
title('Response v. Orientation');
plot(theta * 360 / (2*pi), response_t1, '*');
xlabel('Orientation');
ylabel('Filter Response');
subplot(1, 2, 2); hold on; 
title('Response v. Frequency');
plot(f, response_f1, '*');
xlabel('Frequency');
ylabel('Filter Response')

fig_responses = figure; 
set(fig_responses, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig_responses, 'name', 'Filter Responses Stimulus Two: harmonic = 10, theta = 112.5');
subplot(1, 2, 1); hold on; 
title('Response v. Orientation');
plot(theta * 360 / (2*pi), response_t2, '*');
xlabel('Orientation');
ylabel('Filter Response');
subplot(1, 2, 2); hold on; 
title('Response v. Frequency');
plot(f, response_f2, '*');
xlabel('Frequency');
ylabel('Filter Response');

% test for filter dependence on magnitude of contrast cd .. 
 G_re      = gabor_fn(bw, ar, psi,  lambda(3), orientation(1));
 G_im      = gabor_fn(bw, ar, psi + pi/2,  lambda(3), orientation(1));
 G         = G + 1i*G_im; 
 responses = zeros(1,5);
for k = 1:5
    stim = stim1.*k;
    r            = sum(abs(conv2(stim, G, 'valid')), 'all'); 
    responses(k) = responses(k) + r; 
end


fig_responses = figure; hold on; 
set(fig_responses, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
set(fig_responses, 'name', 'Filter Responses as function of contrast'); 
title('Response v. Orientation');
subplot(1,1,1); hold on; 
plot([1,2,3,4,5], responses, '*'); 
xlabel('Stimulus scalar coefficient');
ylabel('Filter Response');
    
        
        
        
