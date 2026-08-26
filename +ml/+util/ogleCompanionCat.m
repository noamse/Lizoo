function Cat = ogleCompanionCat(FileName, Args)
% Load the OGLE reference catalogue and map it into cut-out pixel coordinates.
% Returned in the [X, Y, Mag] form that ml.util.selectRefSources takes as its
% companion catalogue. The OGLE list reaches fainter and is denser than the
% matched source list, so it finds close companions that the matched list has
% merged into one entry.
% Input  : - Path of the .mat holding the OGLE AstroImage, e.g.
%            OGLEdata/OB260058/OB160058.mat. The first variable in the file is
%            used.
%          * ...,key,val,...
%            'XOffset','YOffset' - Offset subtracted from the OGLE corrX and
%                   corrY before scaling. Default is 230 for both.
%            'ScaleRatio' - OGLE pixel scale divided by the KMT one.
%                   Default is 0.26/0.4.
%            'CentrePix' - Pixel the offset position maps onto.
%                   Default is 150.
%            'ColNameX','ColNameY' - OGLE position columns.
%                   Default 'corrX','corrY'.
%            'ColNameMag' - OGLE magnitude column. Default 'I'.
%            'MaxMag' - Drop entries fainter than this, and with them the
%                   99 and 100 sentinels OGLE uses for a missing magnitude.
%                   Default is 21.
% Output : - An [X, Y, Mag] matrix in cut-out pixels.
% Author : ULTRASAT team (Aug 2026)
% <Unknown Bugs>: the mapping repeats the one written inline in
%   KMT_pipelineI.m, which is where the matched I_ogle and V_ogle come from.
%   The two have to be kept in step; if the pipeline's mapping changes, this
%   one must change with it.
% Example: Cat = ml.util.ogleCompanionCat('~/matlab/Lizoo/OGLEdata/OB260058/OB160058.mat');

    arguments
        FileName
        Args.XOffset     = 230;
        Args.YOffset     = 230;
        Args.ScaleRatio  = 0.26./0.4;
        Args.CentrePix   = 150;
        Args.ColNameX    = 'corrX';
        Args.ColNameY    = 'corrY';
        Args.ColNameMag  = 'I';
        Args.MaxMag      = 21;
    end

    Loaded = load(FileName);
    Vars   = fieldnames(Loaded);
    Obj    = Loaded.(Vars{1});
    if isa(Obj,'AstroImage')
        Tab = Obj.CatData.Table;
    elseif isa(Obj,'AstroCatalog')
        Tab = Obj.Table;
    else
        Tab = Obj;
    end

    X   = (double(Tab.(Args.ColNameX)) - Args.XOffset).*Args.ScaleRatio + Args.CentrePix;
    Y   = (double(Tab.(Args.ColNameY)) - Args.YOffset).*Args.ScaleRatio + Args.CentrePix;
    Mag = double(Tab.(Args.ColNameMag));

    Keep = isfinite(X) & isfinite(Y) & isfinite(Mag) & Mag <= Args.MaxMag;
    Cat  = [X(Keep), Y(Keep), Mag(Keep)];
end
