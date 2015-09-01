function M = ReadArcGrid(filename)

% M = ReadArcGrid('filename') 
%
% Reads the ArcInfo ASCII grid filename.asc into a matrix of grid values
% with associated spatial reference information. The grid is stored as a 
% struct array with elements: 
%
% M.grid      (the matrix of grid values)
% M.ncols     (# of columns in grid)
% M.nrows     (# of rows in grid)
% M.x         (x coordinates of centers of pixels)
% M.y         (y coordinates of centers of pixels)
% M.xllcorner (x coordinate of lower-left element)
% M.yllcorner (y coordinate of lower-left element)
% M.dx        (grid spacing in x direction)
% M.dy        (grid spacing in y direction)
% 
%
% Note: assumes that the .asc file header is in the standard Arc format:
%
% NCOLS xxxxxx
% NROWS xxxxxx
% XLLCORNER xxxxxx
% YLLCORNER xxxxxx
% CELLSIZE xxxxxx
% NODATA_VALUE xxxxxx
%
% If the name of the grid is filename.asc, the input argument can be either
% 'filename.asc' or 'filename'. 

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

if (nargin ~= 1), help(mfilename), return, end

[sizefnr sizefnc] = size(filename);
if strcmpi(filename(sizefnc-3:sizefnc),'.asc')
    filename = filename(1:sizefnc-4);    
end

if ~exist([filename '.txt'],'file')
    error(['File ' filename '.txt does not exist.']);
end

% open the file and check the first string. Should read "ncols". If it doesn't,
% tell the user it isn't a valid file & bail out.

fid = fopen([filename '.txt'],'r');
if ~strcmpi(fscanf(fid, ' %s',1),'ncols')
    fclose(fid);
    error(['File ' filename '.asc is not an ArcInfo ASCII grid file.']);
end

% read in the grid description and calculate cellsize etc.
% we are now right after 'ncols' in the file
M.ncols = fscanf(fid, '%g', 1);
M.nrows = fscanf(fid, '%*s%g', 1);
M.xllcorner = fscanf(fid, '%*s%g', 1);
M.yllcorner = fscanf(fid, '%*s%g', 1);
cellsize = fscanf(fid, '%*s%g', 1);
M.dx = cellsize;
M.dy = cellsize;
nodata_value = fscanf(fid, '%*s%g', 1);
M.x = M.xllcorner + M.dx * (0.5 + (0:M.ncols-1));
M.y = flipud(M.yllcorner + M.dy * (0.5 + (0:M.nrows-1)'));

% read the grid values into the variable M
M.grid = fscanf(fid,'%g',[M.ncols M.nrows]);
fclose(fid);
M.grid = M.grid';

% convert integer data to double precision if necessary
if isa(M.grid,'integer')
    M.grid = double(M.grid);
end

% Replace nodata values with NaNs
M.grid(M.grid==nodata_value)=NaN;
