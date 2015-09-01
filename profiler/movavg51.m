function[smooth_array] = movavg51(array,pixel_size,wind)

% movavg51.m performs a moving average of the elevation array, within a
% window of user-specified width. It takes as input the raw data array, 
% and the pixel size.  The code returns an array of the same size.  Edges 
% of the array are averaged on a progressively smaller window. 
%
% Code written by Eric Kirby (date unknown, but sometime in late 1998)
%
% USAGE:
%     smooth_array = movavg51(array,pixel_size,wind)
%
% INPUT:
%     array - input vector of data
%     pixel_size - distance between data points
%     wind - moving window size
%
% OUTPUT:
%      smooth_array - output of smoothed data
%
% Assumes units are in meters according to text, though it doesn't really 
% matter as long as entered window size and pixel_size are internally 
% consistent

    loop=1;
        
    %error catching loop to make sure that wind is at least pixel_size:
    while loop,
        if wind >= pixel_size,
            loop=0;
        else
            disp(sprintf('Error:  smoothing window size (%sm) must be >= 3x DEM posting (ie %sm)',num2str(wind),num2str(pixel_size*3)))
            wind = input('Enter the smoothing window size (m): ');
        end
    end
    
    % convert to pixels and test for an odd number, ie make it odd
	wind_pix = round(wind/pixel_size);
	if rem(wind_pix,2) == 0
		wind_pix = wind_pix - 1;
	else
	end
    
    if wind_pix==1; %if smoothing window size is 1, no smoothing is done. exit fcn
        disp('No smoothing done in movavg3.m; smoothing window too small')
        smooth_array=array;
        return
    end

	half_width = (wind_pix-1)/2;

% set length of raw data array and initialize new array
	max_x = length(array);
	smooth_array(1:max_x) = zeros;

% set running mean for center of smooth array
	for i = half_width+1:max_x-half_width
		smooth_array(i) = mean(array(i-half_width:i+half_width));
	end
    
% set running mean for array ends
	smooth_array(1) = array(1);
	for i = 2:half_width
		num = i+(i-1);
		num_2 = (num-1)/2;
		smooth_array(i) = mean(array(i-num_2:i+num_2));
	end

	smooth_array(max_x) = array(max_x);
	for i = max_x-half_width+1:max_x-1
		num = max_x - i;
		smooth_array(i) = mean(array(i-num:i+num));
	end

