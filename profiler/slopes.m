function [s,a]=slopes(chandata,ind)

% slopes.m calculates channel slopes based on elevation and distance from 
% mouth.  It returns slope and area, using a central differencing based 
% on a 3-pixel window. 
% 
% USAGE:
%     [s,a]=slopes(chandata,ind)
% OUTPUT:
%    S, slopes
%    A, drainage area
% INPUT:
%    chandata, channel data matrix (profile.m)
%    ind, index flag (usgscontour.m, stepremover.m)

% slope calculated as central difference except at endpts, where it is 
% forward and backwards difference.  This is a way of smoothing, and keeps
% the slope column the same size.  ind is a vector of indices to calculate 
% slopes over; get it from usgscontour.m Uses unsmoothed chandata (or 
% whatever chandata you put in).

pelev=chandata(ind,2);
dfm=chandata(ind,7);
a=chandata(ind,3);

lg=length(dfm);

s=[];
s(1)=(pelev(2)-pelev(1))/(dfm(2)-dfm(1));

%keyboard
for i=2:lg-1,
    %central differencing:
    s(i)=(pelev(i+1)-pelev(i-1))/(dfm(i+1)-dfm(i-1));

    %forward differencing:
    %s(i)=(pelev(i+1)-pelev(i))/(dfm(i+1)-dfm(i));
end
s(lg)=(pelev(lg)-pelev(lg-1))/(dfm(lg)-dfm(lg-1));
s=s';