% Calculate responses of filter bank in response to stimulus set 

function responses = calcResponses(filters, stimuli)

[filterX, filterY, n_f, n_t] = size(filters);
nStims                       = size(stimuli, 3);

mag        = abs(conv2(stimuli(:,:,1), filters(:,:,1,1), 'valid'));
responses  = zeros(size(mag, 1), size(mag, 2), n_f, n_t); 


for n = 1:nStims
    for j = 1:n_f
        for k = 1:n_t
            % calculate response of filter
            mag                = abs(conv2(stimuli(:,:,n), filters(:,:,j,k), 'valid'));
            r                  = sum(mag, 'all');
            responses(:,:,j,k) = mag;     
        end
    end
end


end