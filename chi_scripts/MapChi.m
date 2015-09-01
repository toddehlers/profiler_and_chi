% MapChi.m
%
% Creates a grid, or map, of chi within a channel network.  Based on the
% method for calculating chi described by Perron and Royden (2012).
%
% Inputs: Flow accumulation and flow direction (D8) grids, in ArcINFO
% gridascii format. These grids must be already projected (e.g., UTM) with
% horizontal units of meters.
% 
% Note: Grid should consist of data where drainage basins (one or more)
% exist.  Calculation of chi is only correct if full, entire basins are
% included in the grids.  In other words, calculations of chi are not
% correct if a channel leaves the grid in one place and re-enters the grid
% in another.  Mask out grid cells outside the basins of interest with
% NODATA values (-9999), at least in the flow accumulation grid.  This
% script requires that there is a buffer at least one cell wide between the
% grid edge and the drainage basin(s).
% 
% Reference: 
%   Perron, J.T., and Royden, L., 2012, An integral approach to
%       bedrock river profile analysis: Earth Surface Processes and Landforms, v.
%       38, n. 6, p. 570-576, doi: 10.1002/esp.3302.
% 
% Calls: ReadArcGrid.m and WriteArcGrid.m (in Perron's TopoTools set of functions,
% http://web.mit.edu/perron/www/files/TopoTools-v1.3.zip)
%
% S. Miller, 5/2013
% Updated to fix edge problems and error in chi calculation by R. DiBiase,
% 3/2014
tic

%% Define parameters and files to be loaded
movern = 0.45;  % m/n or theta.  Ideally, this value is estimated from one or more representative stream profiles from the basin, before being applied here to the whole.
Aref = 10^6;    % Reference drainage area in square meters ("A0" in Perron and Royden).  10^6 is a good value.
Amin = 10^6;    % Minimum drainage area of channel cells in square meters.

dir = 'Q:\data2\badams\cascadia\10m_data\chi';
flowaccfile = 'chi_acc_clp'; % flow accumulation file in ArcInfo ASCII format.  Flow accumulation is in units of grid cells.
flowdirfile = 'chi_dir_clp'; % flow direction file in ArcInfo ASCII format.  Values of flow direction follow convention used by ArcINFO.

%% Load the data
flowacc = ReadArcGrid(strcat(dir,'\',flowaccfile));
flowdir = ReadArcGrid(strcat(dir,'\',flowdirfile));

cellsize = flowacc.dx;  % grid cell size in meters
A = flowacc.grid * cellsize^2;    % Convert flow accumulation to drainage area
fdir = flowdir.grid;


%% Move up channels, calculating the area-function term and inter-pixel distances along the way

% Rank A from largest to smallest.  This process includes a temporary hack,
% because the sort function puts NaN at the high end.
A(isnan(A)) = -Inf;
[B,IX] = sort(A(:),'descend'); % Sort from largest drainage area to smallest. Converts array to a vector first. B is a sorted vector; IX is an index.
A(isinf(A)) = NaN;

% Figure out how many cells are within channels, so we don't have to search through all cells in the grid 
n_chancells = length(B(B>=Amin));   

% Set up some variables
sumAterm = zeros(size(A));      % Create empty matrix for summed area-function term%
m = 1;                          % Counter for vector of between-pixel distances

% Loop through every channel cell
for i=1:n_chancells
    
    % Find indices for the unvisited cell with the largest drainage area
    [j,k]=ind2sub(size(A),IX(i));
        
    % Look to see if downstream cell has already had sumAterm calculated for
    % it.  If not, this cell is the basin outlet and we begin summing
    % drainage area here.  Otherwise, we add the drainage area at this cell
    % do the downstream summed area.
    if fdir(j,k)==1         % flow direction to the E
        if k == size(A,2);
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        elseif sumAterm(j,k+1)==0
            sumAterm(j,k) = (Aref/A(j,k))^movern;   % Calculation if current cell is the basin outlet
        else
            dx=cellsize;    % Calculate distance current cell and downstream cell
            sumAterm(j,k) = sumAterm(j,k+1) + dx*(Aref/A(j,k))^movern; % Calculation if current cell is NOT the basin outlet            
        end
    elseif fdir(j,k)==2     % SE
        if j == size(A,1) || k == size(A,2);
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        elseif sumAterm(j+1,k+1)==0
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        else
            dx=cellsize*sqrt(2);
            sumAterm(j,k) = sumAterm(j+1,k+1) + dx*(Aref/A(j,k))^movern;
        end
    elseif fdir(j,k)==4     % S
        if j == size(A,1);
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        elseif sumAterm(j+1,k)==0
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        else
            dx=cellsize;
            sumAterm(j,k) = sumAterm(j+1,k) + dx*(Aref/A(j,k))^movern;
        end
    elseif fdir(j,k)==8     % SW
        if j == size(A,1) || k == 0;
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        elseif sumAterm(j+1,k-1)==0
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        else
            dx=cellsize*sqrt(2);
            sumAterm(j,k) = sumAterm(j+1,k-1) + dx*(Aref/A(j,k))^movern;            
        end
    elseif fdir(j,k)==16    % W
        if k == 0;
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        elseif sumAterm(j,k-1)==0
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        else
            dx=cellsize;
            sumAterm(j,k) = sumAterm(j,k-1) + dx*(Aref/A(j,k))^movern;            
        end
    elseif fdir(j,k)==32    % NW
        if j == 0 || k == 0;
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        elseif sumAterm(j-1,k-1)==0
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        else
            dx=cellsize*sqrt(2);
            sumAterm(j,k) = sumAterm(j-1,k-1) + dx*(Aref/A(j,k))^movern;            
        end
    elseif fdir(j,k)==64    % N
        if j == 0;
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        elseif sumAterm(j-1,k)==0
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        else
            dx=cellsize;
            sumAterm(j,k) = sumAterm(j-1,k) + dx*(Aref/A(j,k))^movern;            
        end
    elseif fdir(j,k)==128   % NE
        if j == 0 || k == size(A,2);
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        elseif sumAterm(j-1,k+1)==0
            sumAterm(j,k) = (Aref/A(j,k))^movern;
        else
            dx=cellsize*sqrt(2);
            sumAterm(j,k) = sumAterm(j-1,k+1) + dx*(Aref/A(j,k))^movern;            
        end
    end
end

% Calculate average distance in meters between pixel centers
%avgdx = sum(dx)/(m-1);

% Now calculate chi  >>i left this in just for now, but should replace
% sumAterm with chi everywhere in loop above - RAD<<
chi = sumAterm;

% Change all cells that have values of zero nan
chi(chi==0) = nan;

% Create a structure array for chi, which will be used next to save the
% chi as an ascii file for use in ArcINFO
M.grid = chi;                       % the matrix of grid values
M.ncols = flowacc.ncols;            % # of columns in grid
M.nrows = flowacc.nrows;            % # of rows in grid
M.x = flowacc.x;                    % x coordinates of centers of pixels
M.y = flowacc.y;                    % y coordinates of centers of pixels
M.xllcorner = flowacc.xllcorner;    % x coordinate of lower-left element
M.yllcorner = flowacc.yllcorner;    % y coordinate of lower-left element
M.dx = flowacc.dx;                  % grid spacing in x direction
M.dy = flowacc.dy;                  % grid spacing in y direction

toc/60


WriteArcGrid(M,strcat(dir,'\chimap.asc'))
% % Save chi array as a txt file
% reply = input('Do you want to save "chi" to file in gridascii format? Y/N [Y]: ', 's');
% if isempty(reply)
%     reply = 'Y';
% end
% if reply == 'Y' | reply == 'y'
%     WriteArcGrid(M,strcat(dir,'\chimap.asc'))
% end



