function [cum_dist,all_z,new_x,new_y,new_z] = new_order(filename,flag,num)
%
% new_order
% an algorithm to order a 2D point data 
% written by byron adams november 2014
%
% algorithm was written to organize bivariate data assuming one path
% through the data.
%
% INPUTS:
% filename -- the full name of an acsii file created from gridded elevation
%             data exported by arcgis, including extension. this script  
%             expects a string for this input (e.g. 'hoh_dem.txt'). the
%             ascii must only contain elevations from one complete drainage
%             basin
% num -- figure number
%
% OUTPUTS:
% new_dist -- the cummulative distance along a path through the data
% new_z -- the reordered z comonent of the data
  
%
% REQUIRED SCRIPTS:
% just_asciin -- script that reads data from arc ascii files
% livinondaedge -- a simple edge filter
%
% tashi delek!
%
%-------------------------------------------------------------------------%
%
% read in ascii data
	[dem,~,~,~,~] = just_asciin(filename);
%
% change arcgis no data values to matlab no data values (NaN)
    dem(dem < 0) = NaN;
    no_data = NaN;
%
% find the elevataion values along the bounding ridges of the basin
    [~,X,Y,Z] = livinondaedge(dem,no_data);
%
% concatonate xy data into a matrix
    data = [X Y];
%
% calculate the distance between neighboring points
    dist = pdist2(data,data);
%
% find how many points there are
    N = size(data,1);
%
% re-order point indicies based on shortest distance to next closest point
    result(1) = 1;
    for ii = 2:N
        dist(:,result(ii - 1)) = Inf;
        [~,closest_idx] = min(dist(result(ii - 1),:));
        result(ii) = closest_idx; %#ok<*AGROW>
    end
%
% re-order xyz data according to new sorting
    for j = 1:length(X)
        new_x(j) = X(result(j)); %#ok<*SAGROW>
        new_y(j) = Y(result(j));
        new_z(j) = Z(result(j));
    end
%
% store this array of z values
    all_z = new_z;
%
% calculate the distance between points
    new_dist(1) = 0;
    for i = 2:N
        new_dist(i) = sqrt((new_x(i - 1) - new_x(i))^2 + (new_y(i - 1) - new_y(i))^2);
    end
%
% calculate the cumulative sum distance along the path
    cum_dist(1) = 0;
    for i = 2:length(new_dist)
        cum_dist(i) = new_dist(i) + cum_dist(i - 1);
    end
% 
% find the indicies of the points that were not included
    for i = 2:length(new_dist)
        diff_dist(i) = cum_dist(i) - cum_dist(i - 1);
    end
    %
    problem_index = find(diff_dist > 100);
    end_pt = problem_index(1) - 1;
%
% remove these points
     new_x = new_x(1:end_pt);
     new_y = new_y(1:end_pt);
     new_z = new_z(1:end_pt);
%
% plot the ordered profile
    figure(num)
    plot(cum_dist(1:end_pt),new_z(1:end_pt),'.k')
    xlabel('Distance (m)')
    ylabel('Elevation (m)')
%
% use electdata.m to subtract points from ridge elevation before surface
% calculation
if flag == 'y'
    xpoint = [];
    ypoint = [];
    pts = [];
    truth = 'y';
    %
    while truth == 'y'
        clear pt_list list this_x this_y
        disp('Drag cursor over basin divde points to exclude from surface interpolation')
        [pt_list,~,~] = selectdata('return','selected','selectionmode','rect','identify','on','verify','on');
        list = pt_list;
        this_x = cum_dist(list);
        this_y = new_z(list);
        pts = [pts;list];
        xpoint = [xpoint this_x];
        ypoint = [ypoint this_y];
        %
        truth = input('Do you want to select more points?  ','s');
    end   
%
% plot selected points
    figure(num)
    hold on
    plot(xpoint,ypoint,'.m')
%    
% subtract the selected points from the basin divde points
    new_x(pts) = [];
    new_y(pts) = [];
    new_z(pts) = [];
end
%