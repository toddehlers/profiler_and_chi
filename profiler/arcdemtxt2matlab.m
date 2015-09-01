function [metadata,x,y,varargout]=arcdemtxt2matlab(matlab_directory,fprefix,varargin)
%SYNTAX
%       [metadata,x,y,dem,accm]=arcdemtxt2matlab('matlab_directory','fprefix');
%       [metadata,x,y,dem,accm]=arcdemtxt2matlab('matlab_directory','fprefix',2);
%       [metadata,x,y,dem]=arcdemtxt2matlab('matlab_directory','fprefix',1);
%       [metadata,x,y]=arcdemtxt2matlab('matlab_directory','fprefix',0);
%       arcdemtxt2matlab('matlab_directory','fprefix',3);
%
%Loads ascii output from arc, writes .mat file for dem, optionally also saves accumulation data.
%Either:  unedited files fprefixacc.txt and fprefixdem.txt must exist in the local directory.
%Or:  enter something bogus for fprefix input, and code will ask for valid filenames.
%INPUT OPTIONS:  
%   If only input variables are 'matlab_directory','fprefix', or ('matlab_directory','fprefix',2), processes and saves metadata, x, y, DEM and ACCM.
%   If input variables are ('matlab_directory','fprefix',1), processes and saves metadata, x, y, and DEM.
%   If input variables are ('matlab_directory','fprefix',0), processes and saves metadata, x, y only.
%        BUT:  still need *dem.txt to exist with headers.
%   If input variables are ('matlab_directory','fprefix',3), processes and saves metadata, x,
%       y, DEM, and ACCM, but no direct output! variables saved in files only.
%       Use this option to save memory.
% fprefix is a string, so put it in single quotes.
%
% Don't cut headers from the text files!
%
%Code looks for filenames ['matlab_directory',fprefix,'dem.txt'] and optionally ['matlab_directory', fprefix,'dem.prj'] and
% ['matlab_directory',fprefix,'acc.txt']

%Saves files fprefixdemm.mat and fprefixaccm.mat, containing variables "dem" and "accum",
%   to be processed using the ArcMap profiler toolbar and profile51.m .
%   Also saves fprefixmeta.mat:  contains variable "metadata" with info from DEM header, and optionally *.prj file,
%   and saves vectors x and y (easting, northing): contain pixel locations for georeferencing.
%
%   See code to set some flags that control what variable names are saved, etc.
%
%Metadata output:  initially 2 columns, matching general format of both 
% header of gridascii command from arc and also the info in the arc *.prj
% file.  But, 12/05, edited to give single column of strings (for maximum
% flexibility since entries *.prj file can vary).
%

%add a backslash to the directory strings if it isn't already there:
slash = '\';
if strcmp(matlab_directory(length(matlab_directory)),slash),
    mat_workdir = matlab_directory;
else
    mat_workdir = strcat(matlab_directory,slash);
end

process_acc=1;  %1 to process accumulation, 0 to not.
%set flag3 to 1 to save both variables "dem" and "fprefixdem" (taking up hard drive space, but being able to load multiple dems into matlab at once w/o renaming them)   
flag3=0;

%flag1 and 2 are changed within code when needed,ie don't change them here:
flag1=0;
flag2=0;
flag4=0;    %flag4 is 0 if .prj exists; 1 if not (will be changed in fcn--don't change here)
flag5=0;    %changed to 1 in code if only to process metadata

disp(sprintf('Function searches for files %s, %s, %s, and loads/saves them',[fprefix,'dem.txt'],[fprefix,'dem.prj'],[fprefix,'acc.txt']));
disp('  to matlab, for use with the ArcMap tool and profile#.m .')
disp('  If files do not exist, user is queried for filenames.')

% keyboard

if length(varargin)==1, %probably case where only metadata will be processed.
    flag6=varargin{:};
    if flag6==0,
        flag5=1;
        disp('Only processing metadata, not exporting dem or accum.')
    elseif flag6==1,
        disp('Only processing and exporting metadata and dem, not accum.')
        process_acc=0;
    elseif flag6==2,
        disp('Processing metadata, dem and accum.')
    elseif flag6==3,
        disp('Processing metadata, dem and accum; memory-saving mode.')    
    end
else
    flag6=2;  %case of processing everything; 
end

filedem=strcat(mat_workdir,fprefix,'dem.txt');      %check for dem file existence
while length(dir(filedem))==0,
    disp(' ');
    disp(sprintf('DEM file %s not found in specified directory!',filedem));
    filedem=input('Enter complete DEM filename (or CTRL-c to quit):  ','s');
    flag1=1;
end

projfile=[mat_workdir,fprefix,'dem.prj'];   %check for prj file existence, OR: give option to not load prj, ie you don't need to have a prj file
if length(dir(projfile))==0,
    disp(' ')
    disp(sprintf('Projection file %s not found.',projfile));
    projfile=input('Enter complete PRJ filename (OR leave blank to read projection info from headers):  ','s');
    if isempty(projfile),
        %case where info from *.prj will not be recorded in *meta.mat; ncols, etc and x,y will still be created
        disp('Extracting projection info from headers of dem.txt file exported from ArcGIS.')
        flag4=1;
    end
    flag2=2;
end

if process_acc & ~flag5,    %check for accumulation existence, provided that acc should be processed, and not just processing metadata
    fileacc=strcat(mat_workdir,fprefix,'acc.txt');
    while length(dir(fileacc))==0,
        disp(' ')
        disp(sprintf('Accumulation file %s not found in specified directory.',fileacc));
        fileacc=input('Enter complete path to accumulation file:  ','s');
        flag1=1;
    end
end

if ~flag4,
    raw_proj=textread(projfile,'%s');
end
%read in entire dem including the 1st 6 lines of header, as string:

demhead=textread(filedem,'%s',12);
%keyboard

%this should always be the same, because its from dem text file:
metadata = [demhead(1); demhead(2); demhead(3); demhead(4); demhead(5); demhead(6); demhead(7); demhead(8);...
            demhead(9); demhead(10); demhead(11); demhead(12)];
%this part will vary, because projection files can have different info:
if exist('raw_proj','var'),
    for i=1:length(raw_proj),
        metadata = [metadata; raw_proj(i)];
    end
end

%header info:
ncols=str2double(demhead(2));
nrows=str2double(demhead(4));
xllcorner=str2double(demhead(6));
yllcorner=str2double(demhead(8));
cellsize=str2double(demhead(10));
NODATA_value=str2double(demhead(12));

%calculate x (easting) and y (northing) vectors:
x=((xllcorner+cellsize/2):cellsize:xllcorner+cellsize*(ncols-1)+cellsize/2)';
y=(yllcorner+cellsize/2):cellsize:(yllcorner+cellsize*(nrows-1)+cellsize/2);
y=(fliplr(y))';

if flag1,   %only interactive if original files weren't found.
    if flag5,
        disp('Data will be saved in file fileprefixmeta.mat:')
        fprefix=input('Enter file_prefix:  ','s');
    elseif process_acc==0,
        disp('Data will be saved in files fileprefixdemm.mat and fileprefixmeta.mat:')
        fprefix=input('Enter file_prefix:  ','s');
    else
        disp('Data will be saved in files fileprefixdemm.mat, fileprefixmeta.mat and fileprefixaccm.mat:')
        fprefix=input('Enter file_prefix:  ','s');
    end
end

%save metadata and x,y vectors separately:
eval(['save ',mat_workdir,fprefix,'meta metadata x y'])

if flag5,  %if only saving metadata, exit now
    disp(sprintf('Exiting; metadata saved to %s.',[mat_workdir,fprefix,'meta.mat']))
    return
end

%keyboard
%load dem, using dlmread:
dem=dlmread(filedem,' ',[6 0 nrows+5 ncols-1]);

%clear demmdata demmout

if flag3,       %case to save both variables dem and fprefixdem    
    eval([mat_workdir,fprefix,'dem=dem;']);
    eval(['save ',mat_workdir,fprefix,'demm ',mat_workdir,fprefix,'dem dem'])
    eval(['clear ',mat_workdir,fprefix,'dem dem'])
else
    eval(['save ',mat_workdir,fprefix,'demm dem'])
    if flag6==3,  %memory-saving mode
        clear dem
        %pack
    end
end

%------------------------------------------------------------------
accm=[];
if process_acc,  %change to zero to not proccess accumulation>
    
    %header info should be identical:
    disp('Not recording header or projection data for Accumulation; should be same as DEM')
    
    %load accum, assuming same number of rows and columns of data as in dem:
    accum=dlmread(fileacc,' ',[6 0 nrows+5 ncols-1]);
    
    %     
    if flag3,  %as above, save accum and fprefixaccum--assumes fprefix is save as for 
        eval([mat_workdir,fprefix,'accum=accum;']);
        eval(['save ',mat_workdir,fprefix,'accm ',mat_workdir,fprefix,'accum accum'])
        eval(['clear ',mat_workdir,fprefix,'accm accum'])
    else    %case for just saving "accum"
        eval(['save ',mat_workdir,fprefix,'accm accum'])
        if flag6==3,
            clear accum
            %pack
        end
    end
    
end

if flag6==3,
else
    varargout(1)={dem};
    if process_acc==1,
        varargout(2)={accum};
    end
end

%keyboard