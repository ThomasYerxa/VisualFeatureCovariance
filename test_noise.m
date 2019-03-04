% Test filters using mkFract to generate noise w/ known frequency
% distribution

filterSize  = 80;
plotting    = true; 
remove_bias = false; 
% largest possible frequency range given filterSize 
f      = linspace(1/filterSize,0.5, 9);
f      = f(2:end - 1);
theta  = [0.0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5] * 2 * pi /360;

G      = filterBank(f, theta, filterSize, false); 

% desgin test stimuli
fract_dim = 1;
stim      = mkFract([512, 512], fract_dim); 

% plot stimulus 
if true
    figure; hold on; 
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
    set(gcf, 'name', 'Stimulus 1');
    imagesc(stim);
    axis image off;
    colormap gray; 
end

responses = calcResponses(G, stim);

%Responses squared give power spectrum 
[PDF_j, PDF_t, PDF_f, PDF_s] = calcPDF(responses.^2, f, theta); 

if remove_bias
    PDF_f = PDF_f(1:end - 1);
    PDF_f = PDF_f ./ sum(PDF_f, 'all');
    f = f(1:end-1);
end

% calculate line of best fit (specifically its slope in log-log space)
p_f       = polyfit(log(f), log(PDF_f), 1);
fit_slope = p_f(1);
Xf        = linspace(log(f(1)), log(f(end)), 100);

% plot results
figure; hold on; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);

c1 = plot(log(f), log(PDF_f), '*'); hold on;    l1 = 'Measured Log Probability';
c2 = plot(Xf, p_f(2) + p_f(1).*Xf);             l2 = 'Best Fit, Slope =  ' + string(fit_slope);
legend([c1,c2], [l1,l2]);

% Use Latex to display title
syms f
f  = 1/f^(5-fract_dim*2);

title(['Response to ', '$', latex(f),'$ noise'], 'Interpreter', 'latex');
xlabel('Log Frequency');
ylabel('Log Probability');

figure; 
plot(theta, PDF_t);

% make figure
figure; hold on;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
f      = linspace(1/filterSize,0.5, 9);
f      = f(2:end - 1);
theta  = [0.0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5] * 2 * pi /360;
G      = filterBank(f, theta, filterSize, false); 

% plot gabors
subplot(4,4,1);
imagesc(real(G(:,:,1,1))); axis image off; 
subplot(4,4,2);
imagesc(real(G(:,:,1,2))); axis image off; 
subplot(4,4,5);
imagesc(real(G(:,:,2,1))); axis image off; 
subplot(4,4,6);
imagesc(real(G(:,:,2,2))); axis image off; 


PDF_test_theta = zeros([size(PDF_t), 10]);
PDF_test_f     = zeros([size(PDF_f), 10]);
for i = 1:10
    stim      = mkFract([512, 512], fract_dim); 
    responses = calcResponses(G, stim);

    %Responses squared give power spectrum 
    [PDF_j, PDF_t, PDF_f, PDF_s] = calcPDF(responses.^2, f, theta); 
    PDF_test_theta(:,:, i)       = PDF_t; 
    PDF_test_f(:,:, i)           = PDF_f; 
     
end

theta_means = mean(PDF_test_theta, 3);
f_means     = mean(PDF_test_f, 3);
theta_error = std(PDF_test_theta, 0, 3) ./ sqrt(10); 
f_error     = std(PDF_test_f,0, 3) ./ sqrt(10);

p_f       = polyfit(log(f), log(f_means), 1);
fit_slope = p_f(1);


subplot(4, 4, [9, 10, 13,14]);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
%set(gcf, 'name', 'Stimulus 1');
imagesc(stim);
axis image off;
colormap gray;  

subplot(4,4,[3,4,7,8]); 
c1 = errorbar(log(f), log(f_means), f_error); hold on;    l1 = 'Measured Log Probability';
c2 = plot(Xf, p_f(2) + p_f(1).*Xf);  l2 = 'Best Fit, Slope =  ' + string(fit_slope);
legend([c1,c2], [l1,l2]);
ylabel('Log Probability');
xlabel('Log Frequency');


subplot(4,4,[11,12,15,16]); 
c3 = errorbar(theta, theta_means, theta_error);
xlim([-.01,3]);
ylim([0, 0.5]);
ylabel('Probability');
xlabel('Orientation (Radians)');

    






