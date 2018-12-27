% calculate and visualize PDF from saved files.
% assumes that in the working directory there are files specifying: 
% 1) partition_num # of partitions of mag data
% 2) files of mag data in the form "mag#.mat" where # is a number between 1
%    and partition_num. 
% 3) the frequency and orientation settings used to generate the mag data


% load settings data
partition_num = load('partition_num');
frequency     = load('f');
theta         = load('theta');
 
% load mag data. Sum each summed version of mag.
mag = zeros(size(load('mag1')));
for i = 1:partition_num
   mag = mag + load('mag' + string(i)); 
end

% calc and visualize final PDFs.
[PDF_j, PDF_t, PDF_f, PDF_s] = calcPDF(mag, frequency, theta); 
vizPDF(PDF_j, PDF_t, PDF_f, PDF_s, frequency, theta);
