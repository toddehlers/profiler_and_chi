function asciirout(filename1,filename2,matrix)
%
% asciirout
% no frills script to to write an arcgis ascii file
% written by byron adams january 2015
%
% INPUTS:
% filename1 -- the full name of an acsii file created from gridded acrgis
%              data (including extension) that has the same spatial data as 
%              the new fiel to be written. this script expects a string for 
%              this input (e.g. 'hoh_dem.txt")
% filename2 -- the full name of the acsii file to be created (including 
%              extension). this script expects a string for this input
%              (e.g. 'hoh_relief.txt")
% matrix -- any matrix of data that is to be written to the file. this
%           matrix may include NODATA values
%
%
% tashi delek!
%
%-------------------------------------------------------------------------%
%
% read in the header from filename1
    header = textread(filename1,'%s',12); %#ok<DTXTRD>
%      
% record useful header info
    ncols = str2double(header(2));
    nrows = str2double(header(4));
    xllcorner = str2double(header(6));
    yllcorner = str2double(header(8));
    cellsize = str2double(header(10));
    NODATA = str2double(header(12));
%
% set all matlab no data values to arcgis no data values
    matrix(isnan(matrix)) = NODATA;
%
% write a new ascii file with the header data from filename1 with the new
% data within matrix
    fid = fopen(filename2,'w');
    fprintf(fid,'%s %d\n%s %d\n%s %f\n%s %f\n%s %f\n%s %d\n','ncols',ncols,'nrows',nrows,'xllcorner',xllcorner,'yllcorner',yllcorner,'cellsize',cellsize,'NODATA',NODATA);
    fclose(fid);
    % 
    dlmwrite(filename2,matrix,'-append','delimiter',' ','precision','%d');
%