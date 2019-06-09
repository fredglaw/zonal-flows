function data_os = manual2Doversample(data_h,N,M)
% Routine to manually oversample 2D data.
% Tried nesting calls to interpft but was unable to get it to work
% Assumes input data is in FFT ordering
% 
% Args: data_h -- 2D array of data to do anti-aliasing on, both dimensions.
%                 Assumes input is in Fourier space, FFT ordering
%       N -- size of the square array data
%       M -- number of nodes to use to anti-alias
% Output: data_os -- oversampled 2D array, in real space
% 
if M <= N
    data_os = ifft2(data_h);
    if M < N
        disp('Error, M < N');
    end
else
    data_h = fftshift(data_h); %shift to a natural ordering
    is_even = rem(N,2) == 0; %figure out if N is even, if so manually put in conjugates
    
    %in the case of even, manually replicate the unmatched modes (row +
    %column) by conjugating, also split by 2 or 4
    if is_even
        bigger_array = zeros(N+1);
        bigger_array(1:N,1:N) = data_h;
        
        corner_term = real(data_h(1,1))/4; %replicate in corners
        first_row = data_h(1,2:end)/2; %repeat this term
        first_col = data_h(2:end,1)/2; %repeat this term
        
        %fill corner terms
%         bigger_array(1,end) = conj(corner_term);
%         bigger_array(end,1) = conj(corner_term);
        bigger_array(1,end) = corner_term;
        bigger_array(end,1) = corner_term;
        bigger_array(end,end) = corner_term;
        bigger_array(1,1) = corner_term; %divide original by 4
        
        %fill last row and column
        bigger_array(2:end-1,end) = flip(conj(first_col));
        bigger_array(end,2:end-1) = flip(conj(first_row));
%         bigger_array(2:end-1,end) = first_col;
%         bigger_array(end,2:end-1) = first_row;
        bigger_array(2:end-1,1) = first_col; %divide original by 2
        bigger_array(1,2:end-1) = first_row; %divide original by 2
        
        data_h = bigger_array;
    end
    
    %now zero pad
    num_to_add = floor((M-N-1)/2);
    data_h = padarray(data_h,[num_to_add,num_to_add],0,'both');
    
    if rem(M,2)==0 %if M is even, pad in one extra layer of zeros
        data_h = padarray(data_h,[1,1],0,'pre');
    end
    data_h = ifftshift(data_h); %inverse shift
    data_os = (((M/N))^2)*ifft2(data_h);
end


end