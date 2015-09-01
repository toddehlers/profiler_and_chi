function [gridout,x,y,value] = livinondaedge(gridin,no_data)
%
% livinondaedge
% a simple edge filter 
% written by byron adams november 2014
%
% this filter is desinged to find the outlines of shapes in gridded data
% surrounded by no data values. it retains the values of positions around
% the perimeter of the shapes and writes no data values to the interior
% postions. the filter only searches for data in adjacent positions.
%
% INPUTS:
% gridin -- any matrix of data where the shape of interest is surrounded by
%           no data values
% no_data -- value that denotes a matrix element contains no data
%
% OUTPUTS:
% gridout -- matrix of size gridin where only the data values at the shape
%            outline are recorded
% x and y -- matrix indicies where the shape outline exist
% value -- value of the matrix where the shape outline exists (i.e. at 
%          position x,y) 
%
%
% tashi delek!
%
%-------------------------------------------------------------------------%
%
% find the size of the matrix
    [nrows,ncols] = size(gridin);
%
% reset all no data pixels to matlab no data values
    gridin(gridin == no_data) = NaN;
%
% add 1 pixel buffer to ensure the object is surrounded by no data
    dummy1 = NaN(1,ncols+2);
    dummy2 = NaN(nrows,1);
    %
    new_grid = [dummy2 gridin dummy2];
    new_grid = [dummy1;new_grid;dummy1];
%
% intialize output matrix
    gridout = NaN(nrows+2,ncols+2);
%
% find matrix postions where real values exist adjacent to no data values 
% and write them to the output matrix
    for r = 2:nrows
        for c = 2:ncols
            if isnan(new_grid(r-1,c)) || isnan(new_grid(r,c-1)) || isnan(new_grid(r+1,c)) || isnan(new_grid(r,c+1))
                gridout(r,c) = new_grid(r,c);
            else
                gridout(r,c) = NaN;
            end
        end
    end
%
% find the indicies of shape outline
    [y,x] = find(gridout >= 0);
%
% find the matrix values at these indicies
    for i = 1:length(x)
        value(i) = gridout(y(i),x(i)); %#ok<*AGROW>
    end
%
% reset x and y indicies
    x = x - 1;
    y = y - 1;
%
% discard synthetic buffer
    gridout = gridout(2:nrows+1,2:ncols+1);
%