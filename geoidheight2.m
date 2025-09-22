function N = geoidheight( lat, lon, varargin )
%  GEOIDHEIGHT Implement a geopotential model to calculate geoid height
%   N = GEOIDHEIGHT( LAT, LON, MODEL ) calculates the geoid height as
%   determined from a selected geopotential model, MODEL. An array of M
%   geoid heights, N, are interpolated at M geodetic latitude, LAT, and M
%   longitude, LON, from a grid of point values in the tide-free system,
%   using the selected geopotential model, MODEL.
%
%   The GEOIDHEIGHT function calculates geoid heights to 0.01 meters for
%   EGM96 and custom.  The GEOIDHEIGHT function calculates geoid heights to
%   0.001 meters for EGM2008.
% 
%   The interpolation scheme wraps over the poles to allow for geoid height
%   calculations at and near these locations.
%
%   Alternate formats for calling geoid height are:
%   N = GEOIDHEIGHT( LAT, LON )
%   N = GEOIDHEIGHT( LAT, LON, ACTION )
%   N = GEOIDHEIGHT( LAT, LON, MODEL, ACTION )
%   N = GEOIDHEIGHT( LAT, LON, 'Custom', DATAFILE )
%   N = GEOIDHEIGHT( LAT, LON, 'Custom', DATAFILE, ACTION )
%
%   Inputs for GEOIDHEIGHT are:
%   LAT      :an array of M geodetic latitudes in degrees where north
%             latitude is positive, and south latitude is negative. LAT
%             must be of type single or double. If LAT is not in the range
%             of [-90,90] it is wrapped to be within the range.
%   LON      :an array of M longitude in degrees where east longitude is
%             positive, west is negative. LON must be of type single or
%             double. If LON is not in the range of [0,360] it is wrapped
%             to be within the range.
%   MODEL    :a string specifying the geopotential model: 'EGM2008' (Earth),
%             'EGM96' (Earth), or 'Custom'.  The default is 'EGM96'.  
%             'EGM96' uses a 15-minute grid of point values in the
%             tide-free system, using EGM96 Geopotential Model to degree
%             and order 360.  The EGM2008 MODEL uses a 2.5-minute grid of
%             point values in the tide-free system, using the EGM2008
%             Geopotential Model to degree and order 2159. The geoid
%             undulations for EGM96 and EGM2008 are with respect to the
%             WGS84 ellipsoid. 
%   DATAFILE :a mat-file containing an array of geodetic
%             latitude breakpoints, 'latbp', an array of longitude
%             breakpoints, 'lonbp', a table of geoid height values,'grid'
%             and a even integer scalar greater than 2 for number of
%             interpolation points, 'windowSize'.  This is only needed for
%             a 'Custom' geopotential model.
%   ACTION   :a string to determine action for out-of-range latitude or
%             longitude. Specify if out-of-range input invokes a 'Warning',
%             'Error', or no action ('None'). The default is 'Warning'.
%
%   Output for GEOIDHEIGHT is:
%   N        :an array of M geoid heights in meters with the same data type
%             as the input LAT.
%
%   Limitations:
%
%   This function using the 'EGM96' MODEL has the limitations of the 1996
%   Earth Geopotential Model.  For more information see the documentation.
%
%   The WGS84 EGM96 geoid undulations have an error range of +/- 0.5 to
%   +/- 1.0 meters worldwide.
%
%   This function using the 'EGM2008' MODEL has the limitations of the 2008
%   Earth Geopotential Model.  For more information see the documentation.
%
%   Examples:
%
%   Calculate the EGM96 geoid height at 42.4 degrees N latitude and 71.0 degrees 
%   W longitude with warning actions:
%       N = geoidheight( 42.4, -71.0 )
%
%   Calculate the EGM2008 geoid height at two different locations with
%   error actions.
%       N = geoidheight( [39.3, 33.4], [77.2, 36.5], 'egm2008','error')
%
%   Calculate a custom geoid height at two different locations with
%   no actions.
%       N = geoidheight( [39.3, 33.4], [-77.2, 36.5], 'custom', ...
%           'geoidegm96grid','none')
%
%   Note: This function uses geoid data that can be obtained using the
%   aeroDataPackage command.
%
%   See also GRAVITYWGS84, GRAVITYSPHERICALHARMONIC

%   Copyright 2010-2023 The MathWorks, Inc.

%   References:
%   [1] NIMA TR8350.2: "Department of Defense World Geodetic System
%       1984, Its Definition and Relationship with Local Geodetic Systems."
%   [2] NASA/TP-1998-206861: "The Development of the Joint NASA GSFC and NIMA
%       Geopotential Model EGM96"
%   [3] Pavlis, N.K., S.A. Holmes, S.C. Kenyon, and J.K. Factor, "An Earth
%       Gravitational Model to Degree 2160: EGM2008", presented at the 2008
%       General Assembly of the European Geosciences Union, Vienna,
%       Austria, April 13-18, 2008. 
%   [4] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.  
%
%   National Geospatial-Intelligence Agency (NGA) Office of Geomatics
%   website:
%   https://earth-info.nga.mil/index.php

persistent astgeoiddata astgeoidCustomError

narginchk(2, 5);

checkinputs();

% set default values
model  = 'EGM96';
action = 'warning';
datafile = 'geoidegm96grid.mat';

switch nargin                                                              
    case 3
        if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
            error(message('aero:geoidheight:inputTypeVar3'));
        end
        if strcmpi( varargin{1}, 'custom')
            narginchk(4, 5);
        else
            % check that the 2008 model is in the path	 
            checkegm2008( varargin{1} );
            % assign model or action
            modeloraction( varargin{1} );
        end       
    case 4
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} ) || ...
                (~ischar( varargin{2} ) && ~isstring( varargin{2} )))
            error(message('aero:geoidheight:inputTypeVar4'));
        end
        if strcmpi( varargin{1}, 'custom')
            % custom model with datafile
            model = lower( varargin{1} );
            datafile = varargin{2};
        else
            % check that the 2008 model is in the path	 
            checkegm2008( varargin{1} );
            % assign model and action
            checkmodel( varargin{1} );
            checkaction( varargin{2} );
        end
    case 5
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
                (~ischar( varargin{2} ) && ~isstring( varargin{2} )) || ...
                (~ischar( varargin{3} ) && ~isstring( varargin{3} ))
            error(message('aero:geoidheight:inputTypeVar5'));
        end
        % set action for custom model with datafile
        if strcmpi( varargin{1}, 'custom')
            % custom model with datafile
            model = lower( varargin{1} );
            datafile = varargin{2};
            checkaction( varargin{3} );
        else
            % This option is only for 'custom'
            error(message('aero:geoidheight:wrongModel'));
        end
end

checklatitude();
checklongitude();

if ( isempty(astgeoiddata) || ~strcmp( astgeoiddata.type, model ) || ...
                   ~strcmpi( astgeoiddata.file, datafile ) || ~isempty(astgeoidCustomError))
    % data needs to be initialized
    datafile = char(datafile);
    try
        data = load(datafile);
    catch MECustomFileLoad
        throwAsCaller(MECustomFileLoad)
    end
    data.type = model;
    data.file = datafile;
    if strcmp('custom', data.type)
        % check for the existence of correct variables in mat-file
        fieldsExist = ~isfield(data,{'latbp' 'lonbp' 'grid' 'windowSize'});
        
        if fieldsExist(1)
            error(message('aero:geoidheight:noLatitudeBP', datafile));
        end

        if fieldsExist(2)
            error(message('aero:geoidheight:noLongitudeBP', datafile));
        end

        if fieldsExist(3)
            error(message('aero:geoidheight:noGrid', datafile));
        end

        if fieldsExist(4)
            error(message('aero:geoidheight:noWindowSize', datafile));
        end

        % allow same custom filename to run through check if had previous error
        astgeoidCustomError = 1;
        
        % check data type and sizes of custom data
        if  ~isnumeric( data.latbp ) || ~isnumeric( data.lonbp ) || ...
                ~isnumeric( data.grid ) || ~isnumeric( data.windowSize )
            error(message('aero:geoidheight:notNumeric'))
        end
        if  ~isscalar( data.windowSize )
            error(message('aero:geoidheight:notScalar'))
        end
        if  ~isvector( data.latbp ) || ~isvector( data.lonbp ) 
            error(message('aero:geoidheight:not1DArray'))
        end
         if ~all( size( data.grid ) == [ numel(data.latbp) numel(data.lonbp) ] ) 
            error(message('aero:geoidheight:wrongMatrixSize'))
        end
        if  any(~isfinite( data.latbp )) || any(~isfinite( data.lonbp )) || ...
             any(any(~isfinite( data.grid ))) || any(~isfinite( data.windowSize ))
            error(message('aero:geoidheight:notFinite'))
        end
        if (mod(data.windowSize,2) || (data.windowSize <= 2))
            error(message('aero:geoidheight:notEvenInteger'))
        end
        % successful custom file read
        astgeoidCustomError = [];
    end
    
    if strcmpi(data.type,'egm96')
        astgeoiddata = [];
        astgeoiddata.type = data.type;
        astgeoiddata.file = data.file;
        astgeoiddata.grid = griddedInterpolant({data.latbp,data.lonbp},...
            data.grid,'spline','none');
    else
        astgeoiddata = data;
    end
end

% Get lat data type to set output later
latType = class(lat);

% Convert to single for interpolation
if ~isa(latType,'single')
    lat = single(lat);
    lon = single(lon);
end
    
% Interpolate to find geoid height
if strcmpi(astgeoiddata.type,'egm96')
    N = astgeoiddata.grid(lat,lon);
else
    n = size(lat);
    N  = zeros(n,'single');
    
    halfWindowSize = astgeoiddata.windowSize/2;

    findQuery = @(z,zq,hw) [find(z < zq,hw,'last');...
        find(z == zq); find(z > zq,hw,'first')];

    for k = 1:prod(n)
        % Find window centered at the current query point
        latInd = findQuery(astgeoiddata.latbp(:),lat(k),halfWindowSize);
        lonInd = findQuery(astgeoiddata.lonbp(:),lon(k),halfWindowSize);
        % Interpolate in the 2-D window
        G = griddedInterpolant({astgeoiddata.latbp(latInd),astgeoiddata.lonbp(lonInd)},...
            astgeoiddata.grid(latInd,lonInd),'spline','none');
        N(k) = G(lat(k),lon(k));
    end
end

% Expand to double for unit changes and truncating.
N = double(N);

if strcmp(astgeoiddata.type,'egm2008')
    N = fix(N.*1000)*0.001;
else
    N = fix(N.*100)*0.01;
end

% Cast to desired output data type
N = cast(N,latType);

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs()
        if ~isnumeric(lat)
            % Latitude should be a numeric array.  Otherwise error.
            error(message('aero:geoidheight:latitudeNotNumeric'));
        end
        if ~isnumeric(lon)
            % Altitude should be a numeric array.  Otherwise error.
            error(message('aero:geoidheight:longitudeNotNumeric'));
        end
        if ~all(size(lat) == size(lon))
            error(message('aero:geoidheight:arraySize'));
        end
        if ~isfloat(lat) || ~isfloat(lon)
            error(message('aero:geoidheight:notFloat'));
        end
        if ~strcmpi(class(lat),class(lon))
            error(message('aero:geoidheight:differentInputType'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checklatitude()
        % Wrap latitude and longitude if necessary
        [latwrapped,lat,lon] = ...
            Aero.internal.geodesy.wraplatitude(lat,lon,180);
        if latwrapped
            % Handle messages based on action

            switch action
                case 'none'
                case 'warning'
                    warning(message('aero:geoidheight:warnLatitudeWrap'));
                case 'error'
                    error(message('aero:geoidheight:latitudeWrap'));
            end
        end
    end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checklongitude()
        [lonwrapped, lon] = Aero.internal.geodesy.wraplongitude(lon, 360);

        if lonwrapped
            % Handle messages based on action

            switch action
                case 'none'
                case 'warning'
                    warning(message('aero:geoidheight:warnLongitudeWrap'));
                case 'error'
                    error(message('aero:geoidheight:longitudeWrap'));
            end
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmodel( str )
        switch lower( str )
            case 'egm2008'
                model = lower( str );
                datafile = 'geoidegm2008grid.mat';
            case 'egm96'
                model = lower( str );
                datafile = 'geoidegm96grid.mat';
            otherwise
                error(message('aero:geoidheight:unknownModel'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkaction( str )
        switch lower( str )
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:geoidheight:unknownAction'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function modeloraction( str )
        switch lower( str )
            case 'egm2008'
                model = lower( str );
                datafile = 'geoidegm2008grid.mat';
            case 'egm96'
                model = lower( str );
                datafile = 'geoidegm96grid.mat';
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:geoidheight:unknownString'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 function checkegm2008( str )	 
     if strcmp(str,'egm2008')
         if isempty(which('geoidegm2008grid.mat','-ALL'))
             error(message('aero:geoidheight:noDownloadData',...
                 '<a href="matlab:aeroDataPackage">aeroDataPackage</a>'))
         end
     end
 end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
