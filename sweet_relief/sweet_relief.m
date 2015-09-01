function [basin_relief,new_surface] = sweet_relief(filename,flag,num)
%
% sweet_relief
% a script that calculates geophyiscal relief within a basin
% written by byron adams in decemeber 2014
%
% this script was built to calculate the difference between modern
% elevations within a drainage basin and the elevation of basin ridges. in
% order to make this measurement a 3D surface is calculated using a cubic
% interpolation between ridge elevations within the basin.
%
% INPUTS:
% filename -- the full name of an acsii file created from gridded elevation
%             data exported by arcgis, including extension. this script  
%             expects a string for this input (e.g. 'hoh_dem.txt'). the
%             ascii must only contain elevations from one complete drainage
%             basin
% flag -- y to select data points to exclude from surface interpolation, n
%         for no
%
% OUTPUTS:
% basin_relief -- a matix of geophysical relief values of the input basin
%
% REQUIRED SCRIPTS:
% just_asciin -- script that reads data from arc ascii files
% livinondaedge -- a simple edge filter
% new_order -- plots the drainage basin divide in one y-z plane
% selectdata -- script that allows for interactive selection of data points
%
%
% tashi delek!
%
%-------------------------------------------------------------------------% 
%
% read in ascii data
	[dem,~,~,ncols,nrows] = just_asciin(filename);
%
% change arcgis no data values to matlab no data values (NaN)
    dem(dem < 0) = NaN;
    no_data = NaN;
%
% find the indicies where elevations exist within the dataset
    [y_all,x_all] = find(dem >= 0);
%
% find the elevataion values along the bounding ridges of the basin
    [~,x_ext,y_ext,z_ext] = livinondaedge(dem,no_data);
%
% call new_order to analyze ridge elevations
    [~,~,new_x,new_y,new_z] = new_order(filename,flag,num+3);
%
% plot the dem and ridge elevations
    figure(num)
    hold on
    imagesc(dem-max(max(dem)));
    scatter(new_x,new_y,5,new_z);
    xlabel('Distance (m)')
    ylabel('Distance (m)')
    set(gca,'YDir','reverse');
    axis image
    mycolormap1 = [gray(30);jet(30)];
    colormap(mycolormap1)
    h = colorbar;
	h.Label.String = 'Elevation (m)';
%
% use a 2D spline to interpolate a surface across the basin divide
% elevations
    elev = round(griddata(new_x,new_y,new_z,x_all,y_all,'cubic'));
%
% regrid the surface data
    surface = zeros(nrows,ncols);
    for l = 1:length(x_all)
        pos1 = y_all(l);
        pos2 = x_all(l);
        %
        surface(pos1,pos2) = elev(l);
    end
%
% reset no data values
    surface(surface <= 0) = NaN;
%
% find where the original dem sticks out above the interpolated surface
    [m,n] = find(dem > surface);
%
% find the values at these indicies
    for p = 1:length(m)
        z(p) = dem(m(p),n(p)); %#ok<AGROW>
    end
%
% add the new points to the existing basin divide points
    x_new = [x_ext;n];
    y_new = [y_ext;m];
    z_new = [z_ext z];
%
% use a 2D spline to interpolate a new surface across new points
    new_elev = round(griddata(x_new,y_new,z_new,x_all,y_all,'cubic'));
%
% regrid new surface data
    new_surface = zeros(nrows,ncols);
    for l = 1:length(x_all)
        pos1 = y_all(l);
        pos2 = x_all(l);
        %
        new_surface(pos1,pos2) = new_elev(l);
    end
%
% reset no data values
    new_surface(new_surface <= 0) = NaN;
%
% find where the original dem sticks out above the interpolated surface
    [m,n] = find(dem > new_surface);
%
% find the values at these indicies
    for p = 1:length(m)
        z2(p) = dem(m(p),n(p)); %#ok<AGROW>
    end
%
% add the new points to the existing basin divide points
    x_new2 = [x_new;n];
    y_new2 = [y_new;m];
    z_new2 = [z_new z2];
%
% use a 2D spline to interpolate a new surface across new points
    new_elev = round(griddata(x_new2,y_new2,z_new2,x_all,y_all,'cubic'));
%
% regrid new surface data
    new_surface = zeros(nrows,ncols);
    for l = 1:length(x_all)
        pos1 = y_all(l);
        pos2 = x_all(l);
        %
        new_surface(pos1,pos2) = new_elev(l);
    end
%
% reset no data values
    new_surface(new_surface <= 0) = NaN;
%
%
% calculate basin relief and write it to an ascii file to import into
% arcgis
    basin_relief = new_surface - dem;
    %
    filename2 = 'relief.txt';
    asciirout(filename,filename2,basin_relief);
%
% plot interpolated surface and geophysical relief
    figure(num + 1)
    imagesc(new_surface)
    xlabel('Distance (m)')
    ylabel('Distance (m)')
    set(gca,'YDir','reverse');
    axis image
    mycolormap2 = [ones(1,3);gray(3000)];
    colormap(mycolormap2)
    h = colorbar;
	h.Label.String = 'Elevation (m)';
    %
    figure(num + 2)
    imagesc(basin_relief)
    xlabel('Distance (m)')
    ylabel('Distance (m)')
    set(gca,'YDir','reverse');
    axis image
    caxis([0 max(max(basin_relief))])
    mycolormap3 = [ones(1,3);jet(3000)];
    colormap(mycolormap3)
    h = colorbar;
	h.Label.String = 'Relief (m)';
    caxis([0 1600])
%