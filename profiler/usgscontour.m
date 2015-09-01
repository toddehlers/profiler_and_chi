function [ind,cont_intv]=usgscontour(chandata)

% usgscontour.m is specifically tailored to remove step-like changes in 
% elevation from USGS DEMs.  The code extracts elevation data only from 
% the intersections of the stream channel and the contour lines. It takes 
% DEMs that have contour stair steps and take indices of chandata only at 
% points where contours cross.  The output is index numbers for rows in 
% chandata; also the contour interval used
%
% usgscontour.m calls the following functions: closeto.m, slopes.m, and
% closeto_allvalues.m
%
% Originally Made by Joel Johnson for 10m DEMs of Utah based on USGS quad maps.
% 
% Usage:
%   [ind,cont_intv]=stepremover(chandata)
% Output:
%   ind - index rows
%   cont_intv - contour interval
% Input:
%   chandata

%Optional interactive mode, to PLOT, check and change indices:  
interactive=1;
%interactive=0;

warning off MATLAB:divideByZero

usgsinterval=40;     %contour interval in FEET; usually 40, sometimes 20
disp(sprintf('contour interval is set to %s FEET (converted to m).',num2str(usgsinterval)))
disp('(To change contour interval, edit usgscontour.m, line 17)')

%convert feet to meters, multiply feet by .3048
usgsm=usgsinterval.*.3048;
cont_intv=usgsm;
minelev=min(chandata(:,2));
maxelev=max(chandata(:,2));

%assume min contour isn't less than zero (below sea level): cint or contour intervals:
%these are where the contour steps should be; ie these are the elevations I want data at:
cint=0:cont_intv:maxelev;

%only take pts within the range of minelev to maxelev:
i=find(cint<minelev);
i=max(i);
cint=cint(i+1:length(cint));

%find the chandata indices within a range of each cint value, and pick the highest one:
% this assumes that the correct elevation is at the bottom of each stairstep
%pm (for plus or minus) is the range in meters that pts will be considered.
%pm=0.1;
pm=0.5;

ind=[];
%keyboard
for i=1:length(cint),

    indrange=find(chandata(:,2)>cint(i)-pm & chandata(:,2)<cint(i)+pm);
    if isempty(indrange),
%        ind(i)=min(find(chandata(:,2)>cint(i)-pm));
        ind(i)=closeto(cint(i),chandata(:,2));
    else
        ind(i)=max(indrange);
    end
end

%Remove points in ind that are the same:  this happens if contours are cut out at a cliff or knickpt.
df=diff(ind);
df=[df df(length(df))];
notsameind=find(df);
ind=ind(notsameind);



%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%above is unchanged original.  Now, try and make it smarter:
[slope1,junk]=slopes(chandata,1:size(chandata,1));
% "trick" slopes.m into taking 2nd derivative of elevation, ie slope of slope:
[slope2,junk]=slopes([chandata(:,1) slope1 chandata(:,3:7)],1:size(chandata,1));

%find all local maxima in 2nd derivative (curvature), slope2:  these points will be the bottoms of steps in the profiles!
localmaxind=[];
for i=2:size(chandata,1)-2,
    %either solitary local max, or 2 pts local max together, 2nd case needed for central differencing of slope:
    if ( slope2(i)>slope2(i-1) & slope2(i)>slope2(i+1) ) | ( slope2(i)>slope2(i-1) & slope2(i)==slope2(i+1) & slope2(i+1)>slope2(i+2) ),
        %loop to make sure the local max is above zero:
        if slope2(i)>0,
            localmaxind=[localmaxind; i];
        end
    end
end


%I'm just going to move each ind value to the nearest local maxima:  seems like it will usually work well enough.
%maybe move it to the nearest large step? 
%But--put in condition that points can't move more than a certain vertical distance--otherwise they stey where they are, or move
%only part way to local max?  looking at plots says 3-5m should be max distance allowed.
%%%  max vertical change allowed to go to a local maximum:
max_vertmove=5;
%%%
indcorrected=[];    %indcorrected is the indices that have been shifted to coincide with the nearest local maximum of slope2.
for i=1:length(ind),
    %if ind(i) matches a point in localmaxind that index is used; otherwise ind(i) is changed to be the nearest localmaxind value 
    q=find(abs(chandata(localmaxind,2)-chandata(ind(i),2))<=max_vertmove);  
    rows=closeto_allvalues(ind(i),localmaxind(q));
%    keyboard
    %rows=closeto_allvalues(ind(i),localmaxind);
    if length(rows)>1,
        %if multiple matching values, use the index that has the highest slope2 value, ie a larger step, still within max_vertmove:
        [maxval,ind2]=max(slope2(localmaxind(q(rows))));
        %[maxval,ind2]=max(slope2(localmaxind(rows)));
        rows=rows(ind2);
        indcorrected(i)=localmaxind(q(rows));
    elseif isempty(rows),
        %case where no local maxima were within max_vertmove of ind(i); keep location unchanged
        indcorrected(i)=ind(i);
    else
        indcorrected(i)=localmaxind(q(rows));
    end
end

%keyboard
%check for repeated indices:

removeind=[];
for i=1:length(indcorrected)-1,
    if indcorrected(i)==indcorrected(i+1),
        removeind=[removeind; i];
    end
end
indcorrected(removeind)=[];

%keyboard


if interactive,
%    figure(110), clf
    figure(4), clf
    plot(chandata(:,7),chandata(:,2));
    hold on
    plot(chandata(ind,7),chandata(ind,2),'g.-');
    plot(chandata(indcorrected,7),chandata(indcorrected,2),'r.-');
    xlabel('distance from mouth, m');
    ylabel('elevation, m');
    title('Profile with picked usgscontour.m steps.  Red points will be used')
    disp(' ')
    disp('Check to see that usgscontour.m did a good job (figure 4).')
    disp('green line:  uncorrected, original usgscontour picks')
    disp('red line:  "corrected" contour crossings, all should coincide with a local step base')
    disp(' ')
end

%keyboard

ind=indcorrected;

warning on MATLAB:divideByZero