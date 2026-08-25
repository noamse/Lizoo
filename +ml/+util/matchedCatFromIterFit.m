function [Matched, ResAstrometry] = matchedCatFromIterFit(IF, Args)
% Attach sky coordinates to the sources of an IterFit solution.
% A MatchedSources produced by KMT_pipelineI carries no sky coordinates, the
% pipeline running with AddSkyCoo=false, so ml.scripts.gaiaAstrometryKMT has
% nothing to cross match against. This solves a WCS from the fitted source
% positions themselves with imProc.astrometry.astrometryCore, which pattern
% matches against an astrometric catalogue and fits the flip, rotation and
% scale rather than assuming them, and returns a catalogue with one row per
% IterFit source, in the same order.
% Input  : - An IterFit object holding a solution, i.e. with ParS populated.
%          * ...,key,val,...
%            'RefMag' - Per-source reference magnitude inserted as column 'I',
%                   which ml.scripts.gaiaAstrometryKMT reads. NaN when empty.
%                   Default is [].
%            'Scale' - Pixel scale [arcsec/pix] handed to astrometryCore.
%                   Default is 0.4.
%            'ColNameMag' - Magnitude field of the IterFit object, used as the
%                   catalogue magnitude for the pattern match.
%                   Default is 'MAG_PSF'.
%            'astrometryCoreArgs' - Cell array of arguments passed on to
%                   imProc.astrometry.astrometryCore. Default is {}.
% Output : - An AstroCatalog with X, Y, the magnitude, RA, Dec and I, one row
%            per IterFit source and in the same order, so that it can be
%            handed to ml.scripts.gaiaAstrometryKMT.
%          - The result structure returned by astrometryCore, whose Success
%            field and WCS.ResFit describe how well the match went.
% Author : ULTRASAT team (Aug 2026)
% Example: [Matched,Res] = ml.util.matchedCatFromIterFit(IFsys, 'RefMag',Info.SrcData.I_ogle);

    arguments
        IF;
        Args.RefMag             = [];
        Args.Scale              = 0.4;
        Args.ColNameMag         = 'MAG_PSF';
        Args.astrometryCoreArgs = {};
    end

    % The fitted mean positions, which correspond one to one with the sources
    X   = IF.ParS(1,:).';
    Y   = IF.ParS(2,:).';
    Mag = IF.medianFieldSource({Args.ColNameMag});

    Cat = AstroCatalog({[X, Y, Mag]}, 'ColNames', {'X','Y',Args.ColNameMag}, ...
                                      'ColUnits', {'pix','pix','mag'});

    RAdeg  = IF.CelestialCoo(1) .* 180 ./ pi;
    Decdeg = IF.CelestialCoo(2) .* 180 ./ pi;

    [ResAstrometry, Matched] = imProc.astrometry.astrometryCore(Cat, ...
        'RA', RAdeg, 'Dec', Decdeg, 'CooUnits', 'deg', ...
        'Scale', Args.Scale, 'CatColNamesMag', Args.ColNameMag, ...
        Args.astrometryCoreArgs{:});

    % astrometryCore leaves Success empty when the pattern match finds nothing,
    % and ~[] is empty, so testing Success alone lets a failure through
    IsSolved = ResAstrometry.Nsolutions > 0 && ~isempty(ResAstrometry.Success) && ResAstrometry.Success;
    if ~IsSolved
        error('matchedCatFromIterFit:NoWCS', ...
              'astrometryCore did not converge on a WCS (%d candidate solutions)', ResAstrometry.Nsolutions);
    end

    % gaiaAstrometryKMT reads the reference magnitude as column 'I'
    if isempty(Args.RefMag)
        RefMag = nan(size(X));
    else
        RefMag = Args.RefMag(:);
    end
    Matched.insertCol(RefMag, Inf, {'I'}, {'mag'});
end
