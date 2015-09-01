function [matrix,no_data,cell_size,ncols,nrows] = just_asciin(filename)
%
% just_asciin
% no frills script to read in data from ascii files exported by arcgis
% written by byron adams december 2014
%
% INPUTS:
% filename -- the full name of any acsii file created from gridded acrgis
%             data (including extension). this script expects a string for 
%             this input (e.g. 'hoh_dem.txt")
%
% OUTPUTS:
% matrix -- a matrix containing the values of each pixel in the gridded
%           dataset
% no_data -- value that denotes a matrix element contains no data. this
%            value was determined within arcgis
% cell_size -- the length dimnesion of each matrix element. this value was
%              determined within arcgis
% ncols -- number of columns in the matrix
% nrows -- number of rows in the matrix
%
%
% tashi delek!
%
%-------------------------------------------------------------------------%
%
% read in the header
    header = textread(filename,'%s',12); %#ok<DTXTRD>
%      
% record useful header info
    ncols = str2double(header(2));
    nrows = str2double(header(4));
    cell_size = str2double(header(10));
    no_data = str2double(header(12));
%
% read in main ascii body
    matrix = dlmread(filename,' ',[6 0 nrows+5 ncols-1]);
%