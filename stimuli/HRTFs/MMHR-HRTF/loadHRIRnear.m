function [ ir, locerr, location ] = loadHRIRnear(dbname, az, el, opt )
% loadHRIRnear    get HRIR nearest to given direction
%
%    ir  = loadHRIRnear( model, az, el )
%    ir  = loadHRIRnear( model, az, el, opt )
%    [ir, locerr] = loadHRIRnear(...)
%    [ir, locerr, location] = loadHRIRnear(...)
%
%    returns a MxN head related impulse response
%
%    Inputs
%       model   The name of the database.  Should be the full path
%               of a '.h5' file.
%
%        az     desired azimuth direction
%
%        el     desired elevation direction
%
%        opt    optional string containing keywords
%               recognized keywords: 
%                 'interaural'
%                   by default, coordinates are interpreted as
%                   vertical polar, that is, a constant elevation
%                   forms a cone on the top or the bottom of the head
%                   Alternatively, the 'interaural' option means that
%                   a constant azimuth forms a cone on one of the
%                   ears.
%                 '16k'
%                   return response downsampled to 16k and padded to
%                   1024 samples
%
%    Outputs 
%         ir    a HRIR of size 2x2205 or 8x2205
%               sampling rate of responses is 44100 Hz
%
%         locerr   inside angle difference of the desired direction
%               to the direction of the actual HRIR (NOT IMPLEMENTED)
%
%         location   2-element vector of azimuth and elevation
%               identifier in the database

if nargin<4
    opt = [];
end

plist = double(h5readatt(dbname, '/', 'azel_index'));

regular_coords = true;
resample16k = false;
if strfind(opt, 'interaural')
    regular_coords = false;
end
if strfind(opt, '16k')
    resample16k = true;
end

azr = deg2rad(az);
elr = deg2rad(el);

% place desired coordinates onto unit sphere
if regular_coords
    x = sin(azr)*cos(elr); % left-right
    y = cos(azr)*cos(elr); % front-back
    z = sin(elr);         % height
else
    % coordinates where the poles are aligned with the ears.  Not much
    % different in front but majorly on top of the head and left-right
    % extremes.
    x = sin(azr);
    y = cos(azr)*cos(elr);
    z = cos(azr)*sin(elr);
end
% in either coordinate system, az=0 el=0 means "nose direction" and is
% given the coordinates x=0 y=1 z=0.  We shift the desired point a smidgen
% into that direction so that in between points (eg. 5 degree points) will
% always get rounded towards there (that is, +5 -> +4, -5 -> -4).
x = .99*x;
y = .99*y + .01;
z = .99*z;

% put all points onto sphere (as recorded regular coords)
p2 = deg2rad(plist);
xp = sin(p2(1,:)).*cos(p2(2,:));
yp = cos(p2(1,:)).*cos(p2(2,:));
zp = sin(p2(2,:));

% compute cartesian dist
dist = sqrt((x-xp).^2 + (y-yp).^2 + (z-zp).^2);
[~, idx] = min(dist);
azdb = plist(1, idx);
eldb = plist(2, idx);

dsname = sprintf('/el%d/az%d', eldb, azdb);
%fprintf('Trying to load %s.\n', dsname);
ir = h5read(dbname, dsname);

if regular_coords
    location = [ azdb eldb ];
else
    azout = rad2deg(asin(xp(idx)));
    elout = rad2deg(atan(zp(idx)/yp(idx)));
    if yp(idx)<0
        elout = elout + 180;
    end
    location = [ azout elout ];
end
    
locerr = 0; % XXX TBD

if resample16k
    % pad to 2822 samples, then resample 44.1k -> 16k
    ir = resample([ir zeros(size(ir,1), 2822-size(ir,2))]', 160, 441);
end

% debug
% azback = rad2deg(atan(x/y));
% if y<0
%     azback = azback + 180;
% end
% if azback<0
%     azback = azback + 360;
% end
% fprintf('el: %f, az: %f (%f, %f)\n', rad2deg(asin(z)), azback);

