function [chandata] = profile51_batch(fname,arcmap_directory,matlab_directory,exist_cd,crit_area)
% profile51_batch.m is the automated version of the core stream profile analysis code.
%   It works in concert with movavg51.m, answer_yn.m, slopes.m, closeto.m,
%   closeto_allvalues.m, auto_ks_calc.m, and auto_chn_finder.m
%
% USAGE:
% OPTION 1:
%       profile51_batch('filename','arcmap_directory','matlab_directory','exist_cd',crit_area);
%           exist_cd is a flag that indicates whether a new channel should
%           be extracted, or whether one intends to re-process
%           previously-extracted data; must be 'y' or 'n' (usually 'n')
%           only specify crit_area when running auto channel head selection 
%           (4 arguments means use location_ij.txt for chn pts; 
%           5 arguments mean use the auto channel head selection routine.
% OPTION 2:
%       profile51_batch('filename','arcmap_directory','matlab_directory','exist_cd');
%           this mode assumes you have selected a list of channel head
%           locations in arcmap, written to location_ij.txt.
%
% INPUT:
%       'filename' is the part of the dem or acc file that does NOT contain -accm or -demm.
%           Because filename and directory_name are strings they MUST be written within single quotes.
%           It loads files exported as ASCII from ARC, headers trimmed, and saved as .mat files.
%           (these files are the masked dem and flowaccumulation matrices.  
%           Make these files manually, or use arcdemtxt2matlab.m).
%           filenameaccm.mat and filenamedemm.mat are loaded from the 'matlab_directory'.
%           The files must contain, respectively, variables named "accum" and "dem".
%                                           OR:  
%       'filename' is the prefix of an existing chandata file before _chandata.mat.
%       'arcmap_directory' and 'matlab_directory' are strings that give the complete directory (folder) path (see below).
%            Use '.' as shorthand to specify the local matlab directory.
%       NOTE: In addition to these input parameters, files
%       run_parameters.txt and (in some cases) location_ij.txt must exist in
%       'arcmap_directory'.  These files are created by the profile toolbar
%       in ArcGIS.
% OUTPUT:
%   filename_ksdata.txt: output text data of k_sn data for all channels, for import
%        into ArcMap.
%   chandata:  10 column matrix of DEM points along channel; columns are:
%        chandata = [dfd' pelev' drainarea' smooth_pelev' ptargi' ptargj' dfm' auto_ks_vals' x_coord' y_coord']
%        dfd = distance from divide
%        pelev = elevation
%        drainarea = drainage area
%        smooth_pelv = smoothed elevation
%        ptargi = x-coordinate (matrix value)
%        ptargj = y-coordinate (matrix value)
%        dfm = distance from mouth
%        auto_ks_vals = automatically extracted steepness indices along stream
%        x_coord = geographic coordinate (easting)
%        y_coord = geographic coordinate (northing)
%   Various matlab, postscript, and arcmap-loadable files are also interactively saved.
%
% 'arcmap_directory' is the directory where arcmap writes all stream-related shapefiles.
% IT MUST MATCH THE DIRECTORY SPECIFIED IN ARCMAP ("set parameters to sent to matlab" button).
%
% Matlab writes files for input into arcmap into 'arcmap_directory' (e.g., array_stream.txt, section#.txt).
% 'matlab_directory' is the directory where new matlab files are saved.
%
%


disp('You are using the BATCH profile extraction tool')
disp(' ')
disp('If you specified crit_area this triggers auto channel head location routine')
disp(' ')
disp('      (there is no step removal option in this code)')

% set global varible (tribsection).  This is used to number the
% files associated with the trib sections, which are read in Arcview
global TRIBSECTION

% set tribsection
TRIBSECTION = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%flag to import geographic coordinates, rather than matlab indices:
use_geographic_coords=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist_cd=='y',
    disp('Processing PREVIOUS chandata.')
    %       'exist_cd' is a flag to indicate whether a new chandata file is to be created (this is the default
    %           value: 'n', meaning there is no existing chandata file), or whether a pre-existing chandata file
    %           (with name fname_chandata) is to be processed ('y').
elseif nargin==4, %nargin is the number of inputs to this function
    disp('Processing hand-selected batch, NEW chandata.')
    exist_cd='n';
elseif nargin==5,
    disp('Auto location of channel heads and chandata creation in process.')
    exist_cd='n';
else
    disp('Error--number of input arguments is wrong.  Exiting')
    return
end

%expand directory strings if '.' was used:
if strcmp(arcmap_directory,'.'),
    arcmap_directory=pwd;
end
if strcmp(matlab_directory,'.'),
    matlab_directory=pwd;
end
%add a backslash to the directory strings if it isn't already there:
slash = '\';
if strcmp(arcmap_directory(length(arcmap_directory)),slash),
    arc_workdir = arcmap_directory;
else
    arc_workdir = strcat(arcmap_directory,slash);
end
if strcmp(matlab_directory(length(matlab_directory)),slash),
    mat_workdir = matlab_directory;
else
    mat_workdir = strcat(matlab_directory,slash);
end

%
% *************************************************************************
% Load Data Loop: Dem, Accum for exist_cd = 'n' (no existing chandata),
% Chandata for exist_cd = 'y'
% *************************************************************************
%

%loop to enter correct exist_cd parameter, if its wrong:
while exist_cd ~= 'n' & exist_cd ~= 'y',
    exist_cd=input(sprintf('Use previous chandata file %s?  (y/n)',[fname,'_chandata;']),'s');
end

if exist_cd == 'n'
    % Load files saved to matlab binaries w/ matrices named dem, accum
    try
        eval(['load ',mat_workdir,fname,'demm;']);
    catch
        disp(sprintf('attempt to load dem file %s failed!',[mat_workdir,fname,'demm.mat;']));
        dem_fname = input(sprintf('Enter complete dem filename in %s:  ',mat_workdir),'s');
        eval(['load ',mat_workdir,dem_fname]);
    end

    try
        eval(['load ',mat_workdir,fname,'accm']);
    catch
        disp(sprintf('attempt to load accum file %s failed!',[mat_workdir,fname,'accm.mat;']));
        acc_fname = input(sprintf('Enter complete accum filename in %s:  ',mat_workdir),'s');
        eval(['load ',mat_workdir,acc_fname]);
    end

    if use_geographic_coords,   %load *meta.mat file, that contains x,y--the geographic coordinates corresponding to the dem.
        try
            eval(['load ',mat_workdir,fname,'meta']);
            easting=x;
            northing=y;
            clear x y
        catch
            disp(sprintf('attempt to load metadata file %s failed!',[mat_workdir,fname,'meta.mat;']));
            meta_fname = input(sprintf('Enter complete metadata filename in %s:  ',mat_workdir),'s');
            eval(['load ',mat_workdir,meta_fname]);
            easting=x;
            northing=y;
            clear x y
        end
    end

    %error checking to make sure that variables named dem and accum now exist:
    if isempty(who('accum')),
        disp(sprintf('Error--file %s apparently does not contain variable "accum".',[fname,'accm']));
        %         disp('Manually create variable "accum", or load appropriate accumulation data, or type "dbquit" to exit profile3.m');
        %         keyboard
        disp('Variable "accum" must exist; profile3.m exiting.')
        return
    end
    if isempty(who('dem')),
        disp(sprintf('Error--file %s apparently does not contain variable "dem".',[fname,'demm']));
        %         disp('Manually create variable "dem", or load appropriate DEM data, or type "dbquit" to exit profile3.m')
        %         keyboard
        disp('Variable "dem" must exist; profile3.m exiting.')
        return
    end

    % load trib starting point, cellsize, and reference concavity from
    % arcview/arcgis

    eval(['load ', arc_workdir, 'run_parameters.txt']);
    cellsize = run_parameters(1);
    theta_ref = run_parameters(2);
    rmspike = run_parameters(3);
    no_step = run_parameters(4);
    smooth_prof = run_parameters(5);
    wind = run_parameters(6);
    cont_intv = run_parameters(7);
    ks_window = run_parameters(8);
    Pix_Ran_Dnst = run_parameters(9);
    MinAccum2start = run_parameters(10);
    gridsize = size(accum);
    %
    movernset = -1*theta_ref;
    %
    % Set cellsize for analysis
    pix = round(10*cellsize)/10;
    diag = round(10*(sqrt((pix^2)+(pix^2))))/10;
    ar = round(pix^2);

elseif exist_cd == 'y'
    %
    % Else load existing chandata file; name = fname (for uniform input
    % convenience)
    %
    try
        eval(['load ',mat_workdir,fname,'_chandata;']);
    catch
        disp(sprintf('attempt to load chandata file %s failed!',[mat_workdir,fname,'_chandata.mat;']));
        chan_fname = input('Enter correct chandata filename:  ','s');
        try
            eval(['load ',mat_workdir,chan_fname]);
        catch
            disp(sprintf('Attempt to load chandata file %s failed again, exiting.',[mat_workdir,chan_fname]));
            return
        end
    end
    name = fname;
    name = deblank(name);
    %
    %   read cellsize and reference concavity theta_ref from arcview directory
    %
    eval(['load ', arc_workdir, 'run_parameters.txt'])
    cellsize = run_parameters(1);
    theta_ref = run_parameters(2);
    rmspike = run_parameters(3);
    no_step = run_parameters(4);
    smooth_prof = run_parameters(5);
    wind = run_parameters(6);
    cont_intv = run_parameters(7);
    ks_window = run_parameters(8);
    Pix_Ran_Dnst = run_parameters(9);
    MinAccum2start = run_parameters(10);
    gridsize = size(accum);
    %
    movernset = -1*theta_ref;
    pix = round(10*cellsize)/10;
    %
    ar = round(pix^2);
    %

else
    text = 'No valid option specified. Input parameter exist_cd must be either n or y. Restart'
    return
end
%
% *************************************************************************
% End Load Data Loop
% *************************************************************************
%
crit_pix = crit_area/ar;
ks_table = [];
network = zeros(size(accum));
%
%   if running auto channel head location (crit_area specified), then find
%   channel heads and over-write location_ij.txt
%
%loop to load ALL points in location_ij.txt, FIRST TIME!  after first time, it can be re-read inside loop, see below.

bigloopnum=1;
if nargin==5 %batch-code
    [atargiall,atargjall,atargxall,atargyall,nameall] = auto_chn_finder(fname,arc_workdir,mat_workdir,accum,network,easting,northing,crit_pix);
    numchanpts = length(atargiall); %total number of streams to process
else
    if exist_cd == 'n',
        file = 'location_ij.txt';
        source = strcat(arc_workdir,file);
        %allow reading in and processing of multiple points
            % import format: [row column easting northing  name]
            % note: atargx = northing (rows), atargy = easting (columns)
        [atargiall,atargjall, atargyall,atargxall,nameall] = textread(source,'%f %f %f %f %s',-1);      %-1 flag will load all lines.  Assumes 5 columns
        %         nameall = char(nameall);
        clear file source;
        numchanpts=size(atargiall,1); %number of rows in location_ij.txt, ie number of points to loop through
    end
end

answer1 = 1;
while answer1,
    clear pdist paccum pelev drainarea targi targj ptargi ptargj dfm dfd p_x p_y ks east north area
    %
    TRIBSECTION = 1;
    %
    % *************************************************************************
    % Create Chandata Loop for default, exist_cd = 'n' (no existing chandata)
    % *************************************************************************
    %
    if exist_cd == 'n'
        if nargin==4
            if bigloopnum > numchanpts,   %case to reread location_ij.txt
                bigloopnum=1;
                file = 'location_ij.txt';
                source = strcat(arc_workdir,file);
                %allow reading in and processing of multiple points
                % import format: [row column easting northing  name]
                % note: atargx = northing (rows), atargy = easting (columns)
                [atargiall,atargjall, atargyall,atargxall,nameall] = textread(source,'%f %f %f %f %s',-1);      %-1 flag will load all lines.  Assumes 5 columns
                %                nameall = char(nameall);
                clear file source;
                numchanpts=size(atargiall,1); %number of rows in location_ij.txt, ie number of points to loop through
            end
        end
        atargi=atargiall(bigloopnum);
        atargj=atargjall(bigloopnum);
        atargx=atargxall(bigloopnum);
        atargy=atargyall(bigloopnum);
        %         name=nameall(bigloopnum,:);
        if nargin == 4
            name=char(nameall(bigloopnum));
            disp(sprintf('Processing stream %s out of %i from ArcMap',name,numchanpts));
        elseif nargin == 5
            name = num2str(nameall(bigloopnum));
            disp(sprintf('Processing stream %s out of %i from ArcMap',name,numchanpts));
        end

        if use_geographic_coords,       %convert location_ij just loaded from geographic coords to actually be matlab i,j
            %error check; make sure current point falls inside the bounds of the easting, northing vectors, ie in the bounds of the dem.
            if atargx<min(northing) | atargx>max(northing) | atargy<min(easting) | atargy>max(easting),
                disp('ERROR!  Starting point chosen, geographic coordinates, out of bounds of DEM.  Aborting.')
                warning('hi')
                keyboard
                return
            end
            northindex=closeto(atargx,northing);
            eastindex=closeto(atargy,easting);
            atargj=eastindex;
            atargi=northindex;
        end

        % %Uncomment to see if geographic locations are matching the matlab index locations
        % warning('check atargi,atargj ')
        % keyboard
        %
        targj=round(atargj);
        targi=round(atargi);

        % First, head DOWNSTREAM some # of steps (stepb)to make sure you are in the channel.
        if nargin ~= 5,
            for steps = 0:Pix_Ran_Dnst
                steps = steps+1;
                try
                    patch = accum(targi-1:targi+1,targj-1:targj+1);
                catch
                    warning('Error!  Most likely, channel starting point is not located on the DEM.  Aborting script.')
                    return
                end
                [i_dnst,j_dnst] = find(patch==max(max(patch)));
                targi=targi+(i_dnst-2);
                targj=targj+(j_dnst-2);
            end
        end

        % Now, search back UPSTREAM, following the path of largest upstream drainage areas.
        upst_accum_present = accum(targi(1),targj(1));
        while upst_accum_present > MinAccum2start
            patch = accum(targi-1:targi+1,targj-1:targj+1);
            list = sort(reshape(patch,9,1));
            [i_list,j_list] = find(list==upst_accum_present);
            output_val = list(i_list-1, j_list);
            [i_upst,j_upst] = find(patch==output_val(1));
            targi=targi+(i_upst-2);
            targj=targj+(j_upst-2);
            targi=targi(1);
            targj=targj(1);
            upst_accum_present = patch(i_upst(1),j_upst(1));
        end
        % Now work back DOWNSTREAM from the channel head to the outlet.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        true=1;
        step=1;
        distance = [diag pix diag pix 0 pix diag pix diag];
        % Gather data from the first pixel of the channel
        pdist(step)=0;
        paccum(step)=accum(targi,targj);
        pelev(step)=dem(targi,targj);
        %         plith(step)=rock(targi,targj);
        ptargi(step)=targi;
        ptargj(step)=targj;

        % Begin the channel extraction while loop:
        while true == 1
            step=step+1;
            patch = accum(targi-1:targi+1,targj-1:targj+1);
            patch(2,2)=0;
            if step>2
                patch(2-(targi-ptargi(step-2)),2-(targj-ptargj(step-2)))=0;
            end
            [list,index] = sort(reshape(patch,9,1));
            dnst_accum = list(9);
            [i_dnst,j_dnst] = find(patch==dnst_accum);
            targi=targi+(i_dnst(1)-2);
            targj=targj+(j_dnst(1)-2);

            % now extract data for this new point
            pdist(step)=distance(index(9));
            paccum(step)=dnst_accum;
            pelev(step)=dem(targi,targj);
            %             plith(step)=rock(targi,targj);
            ptargi(step)=targi;
            ptargj(step)=targj;
            p_x(step) = targi;
            p_y(step) = targj;

            % Note: will attempt to step out of bounds and will crash if last
            % point on channel is at edge of dem (i.e., if there is no NODATA
            % buffer)
            if paccum(step) < paccum(step-1)
                true = 0;
            elseif targi >= gridsize(1)
                true = 0;
            elseif targj >= gridsize(2)
                true = 0;
            elseif targi <=1
                true = 0;
            elseif targj <=1
                true = 0;
            end
        end

        % Calculate cumulative distance from divide and convert to distance from the mouth.
        % Flip pelev, dfd  and paccum vectors.
        dfd = cumsum (pdist);
        pelev = fliplr(pelev);
        %         plith = fliplr(plith);
        paccum = fliplr(paccum);
        dfd = fliplr(dfd);
        ptargi = fliplr(ptargi);
        ptargj = fliplr(ptargj);
        p_x = fliplr(p_x);
        p_y = fliplr(p_y);
        dfm = max(dfd)-dfd;
        % Convert drainage area from pix to m^2.  Create dummy data for smooth_pelev
        drainarea = (paccum.*ar);
        smooth_pelev = zeros(1,max(size(pelev)));
        auto_ks_vals = zeros(1,max(size(pelev)));
        x_coord = zeros(1,max(size(pelev)));
        y_coord = zeros(1,max(size(pelev)));

        % make the chandata matrix
        chandata = [dfd' pelev' drainarea' smooth_pelev' ptargi' ptargj' dfm' auto_ks_vals' x_coord' y_coord'];
    end
    %
    % *************************************************************************
    % END Create Chandata Loop for default, exist_cd = 'n' (no existing chandata)
    % *************************************************************************
    % *************************************************************************
    % Chandata now loaded for either option, clear memory and begin analysis
    % *************************************************************************
    %

    % reset variables to ensure all are in the same row vector format.
    dfd = chandata(:,1)';
    pelev = chandata(:,2)';
    drainarea = chandata(:,3)';
    smooth_pelev = chandata(:,4)';
    ptargi = chandata(:,5)';
    ptargj = chandata(:,6)';
    p_x = ptargi;
    p_y = ptargj;
    dfm = chandata(:,7)';
    %
    % if running auto_chn_finder routine, set accum = NaN for pixels
    % already traversed -- this will eliminate over-lapping points.
    % NAH - for many cases will need full chandata files to mouth on many
    % if not all streams - comment this routine out for now.  Below will
    % add routine to trim duplicates out of ks_table.
    %
    if exist_cd == 'n',
        if use_geographic_coords,    %convert j and i locations to x and y locations using the meta.mat file data
            chandata(:,9) = easting(ptargj,1);
            chandata(:,10) = northing(ptargi,1);
        end
    end

    % *************************************************************************
    % Loop to perform analysis starts here
    % *********************************************************************
    %%%%%%%%%%%%
    if rmspike,
        temp_elev = fliplr(pelev);
        q = length(temp_elev)-1;
        for i = 1:q
            if temp_elev(i+1) > temp_elev(i)
                temp_elev(i+1) = temp_elev(i);
            end
        end
        new_pelev = fliplr(temp_elev);
    else
        new_pelev = pelev;
    end
    %
    no_step=0;      %no_step=1 to use a step remover; 0 for no step remover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if no_step,
        smooth_pelev = new_pelev;
        chandata(:,4) = smooth_pelev';
        wind = 0;
    else
        %only allow smoothing if step remover will not be used:  smooth raw
        %dem data
        answer100 = 1;  %1 to do smoothing; 0 no smoothing
        if wind == 0
            answer100 = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if answer100,
            [smooth_pelev] = movavg51(new_pelev,pix,wind);
            % replace with option movavg51a to break into segments, avoid
            % smoothing across tributary junctions
            % [smooth_pelev] = movavg51a(new_pelev,drainarea,pix,wind);
            chandata(:,4) = smooth_pelev';
        else
            smooth_pelev = new_pelev;
            chandata(:,4) = smooth_pelev';
            wind = 0; %no smoothing done
        end
    end
    % end spike removal / step removal / smoothing routine


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % in BATCH profile code, you directly call out to this code
    [chandata] = auto_ks_calc(chandata,movernset,cont_intv,ks_window);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NMG Put in a warning to user if a channel was encountered without enough
    %relief to do contour smoothing
    if chandata(1,8)==-9999;
        disp(sprintf('WARNING Stream %s did not have enough relief to ',name));
        disp(sprintf('calculate ks.  Assigned values of -9999.'));
    end

    %
    ks = chandata(:,8)';
    east = chandata(:,9)';
    north = chandata(:,10)';
    itarg = chandata(:,5)';
    jtarg = chandata(:,6)';
    area = chandata(:,3)';
    count = length(ks);
    namevec = ones(count,1);
    namevec = namevec*bigloopnum;
    if bigloopnum == 1,  % first time through include entire stream
        ks_table=[ks_table;namevec,east',north',area',ks'];
        for c = 1:count
            network(itarg(c),jtarg(c))=1; % mark positions recorded in ks_table
        end
    else
        %
        for c = 1:count
            if network(itarg(c),jtarg(c))==1,
                area(c)=0; %mark points already analysed to be trimmed out
            end
        end
        ks_trim = ks(find(area));
        % this just saves the FIRST one, not the minimum.  edit to achieve
        % that goal
        east_trim = east(find(area));
        north_trim = north(find(area));
        itarg_trim = itarg(find(area));
        jtarg_trim = jtarg(find(area));
        area_trim = area(find(area));
        new_count = length(ks_trim);
        %     if new_count == 0   %debug: print out to id full repeat cases
        %         new_count
        %     end
        namevec = ones(new_count,1);
        namevec = namevec*bigloopnum;
        %
        ks_table=[ks_table;namevec,east_trim',north_trim',area_trim',ks_trim'];
        for d = 1:new_count
            network(itarg_trim(d),jtarg_trim(d))=1;  % mark positions recorded in ks_table
        end
    end % END if bigloopnum = 1, else statement
    %
   
    %
    % Different save filename options depending on whether a pre-existing
    % chandata file is being reprocessed
    %

    if exist_cd == 'n'
        %         t5 = answer_yn('Save this long profile matlab data?');
        t5 = 1;
        if t5,
            %             sa_data = [ida' islope'];
            %             lb_sa_data = [lbarea' lbslope'];
            name = deblank(name);

            eval ([' save ',mat_workdir,name,'_chandata.mat chandata -MAT'])
            %             %        eval ([' save ',mat_workdir,name,'_chandata_ns.mat chandata_ns -MAT'])
            %             eval ([' save ',mat_workdir,name,'_sa_data.mat sa_data -MAT'])
            %             eval ([' save ',mat_workdir,name,'_lb_sa_data.mat lb_sa_data -MAT'])
        end
       
        bigloopnum=bigloopnum+1;
        if bigloopnum <= numchanpts,    
%          
            answer1 = 1;
        else
            disp('No more points chosen for processing, writing auto_ks to text file ...')
       
            answer1 = 0;
        end

    else if exist_cd == 'y'

            disp('somethings wrong: the save routine in code thinks you have existing chandata')

        end
        %
        % End Save Options Loop
        %
    end
    %
    % End While Loop for repeat analysis if answer1 = 1,
    %
end
%
% write out ks_data.txt for import into arc then END profile42_batch.m
%
[v,d]=version;
if str2num(v(1))==6
    dlmwrite([arc_workdir,fname,'ks_data.txt'], ks_table,' ');
else
    dlmwrite([arc_workdir,fname,'ks_data.txt'], ks_table,'delimiter',' ','precision','%1.14g');
end
disp(sprintf('Finished writing %sks_data.txt, ready for import auto_ks values into arcgis',fname))