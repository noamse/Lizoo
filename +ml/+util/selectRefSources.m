function [Flag, Info] = selectRefSources(Obj, Args)
% Select a clean subset of stars to define the astrometric reference frame.
% The per-epoch transformation in @IterFit is fitted from every source, each
% weighted by the median of its weights over epochs. Faint and blended stars
% therefore help define the frame that everything is then measured against.
% This picks a small, bright, isolated subset instead: the frame is fitted from
% them alone while the source parameters are still solved for every star.
% Input  : - An MMS or IterFit object holding the sources.
%          * ...,key,val,...
%            'RefMag' - Reference magnitude per source, typically the OGLE I
%                   from SrcData.I_ogle. Sources without one cannot be
%                   selected. Mandatory.
%            'MagRange' - Keep sources whose reference magnitude falls in this
%                   range. Bright enough to be measured well, faint enough not
%                   to saturate. Default is [14 16].
%            'CompanionRadius' - A source is rejected if another star brighter
%                   than CompanionMaxMag lies within this distance [pix].
%                   Default is 5.
%            'CompanionMaxMag' - Companions fainter than this are ignored.
%                   Default is 18.
%            'CompanionCat' - An [X, Y, Mag] matrix of an external catalogue to
%                   search for companions, for instance the OGLE reference,
%                   which reaches fainter than the matched source list. The
%                   object's own sources are always searched as well.
%                   Default is [].
%            'SelfRadius' - A star appears in an external catalogue at its own
%                   position, offset by the matching accuracy, and must not be
%                   counted as its own companion. The nearest catalogue entry
%                   within this distance [pix] is taken to be the star itself
%                   and ignored. Should be about the radius the catalogue was
%                   matched at, 1.5 pix for the OGLE match in KMT_pipelineI.
%                   Default is 1.5.
%            'ColNameX','ColNameY' - Position fields. Default 'X','Y'.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - A logical row vector over the object's sources.
%          - A structure recording how many sources each criterion removed, the
%            positions and magnitudes used, and the nearest-companion distance
%            of every candidate.
% Author : ULTRASAT team (Aug 2026)
% Example: [Flag,Info] = ml.util.selectRefSources(Obj, 'RefMag',Info.SrcData.I_ogle, ...
%                                                 'CompanionCat',[Xo Yo Io]);

    arguments
        Obj;
        Args.RefMag;
        Args.MagRange         = [14 16];
        Args.CompanionRadius  = 5;
        Args.CompanionMaxMag  = 18;
        Args.CompanionCat     = [];
        Args.SelfRadius       = 1.5;
        Args.ColNameX         = 'X';
        Args.ColNameY         = 'Y';
        Args.Verbosity        = 0;
    end

    Nsrc = Obj.Nsrc;
    X    = median(Obj.Data.(Args.ColNameX), 1, 'omitnan').';
    Y    = median(Obj.Data.(Args.ColNameY), 1, 'omitnan').';
    Ref  = Args.RefMag(:);
    if numel(Ref) ~= Nsrc
        error('selectRefSources:BadRefMag','RefMag must hold one value per source (%d)', Nsrc);
    end

    % --- magnitude window -------------------------------------------------
    InMag = isfinite(Ref) & Ref >= Args.MagRange(1) & Ref <= Args.MagRange(2) & isfinite(X) & isfinite(Y);

    % --- companions -------------------------------------------------------
    % Searched both in the object's own sources and, when given, in an external
    % catalogue: the matched source list does not resolve every close blend.
    CompX = X; CompY = Y; CompM = Ref;
    if ~isempty(Args.CompanionCat)
        CompX = [CompX; Args.CompanionCat(:,1)];
        CompY = [CompY; Args.CompanionCat(:,2)];
        CompM = [CompM; Args.CompanionCat(:,3)];
    end
    IsComp = isfinite(CompM) & CompM <= Args.CompanionMaxMag & isfinite(CompX) & isfinite(CompY);
    CompX  = CompX(IsComp); CompY = CompY(IsComp);

    NearestD = inf(Nsrc,1);
    NSelf    = 0;
    Cand     = find(InMag);
    for Ic = 1:numel(Cand)
        Isrc = Cand(Ic);
        D = hypot(CompX - X(Isrc), CompY - Y(Isrc));
        % The star's own entry in the object sits at zero distance, and its
        % counterpart in an external catalogue sits within the matching
        % accuracy. Drop the nearest such entry: it is the star itself.
        D(D < 1e-6) = Inf;
        [Dmin, Imin] = min(D);
        if Dmin < Args.SelfRadius
            D(Imin) = Inf;
            NSelf   = NSelf + 1;
        end
        NearestD(Isrc) = min(D);
    end
    Isolated = NearestD > Args.CompanionRadius;

    Flag = (InMag & Isolated).';

    Info = struct('Nsrc',Nsrc, 'NinMag',sum(InMag), 'NafterIsolation',sum(Flag), ...
                  'NrejectedByCompanion',sum(InMag & ~Isolated), ...
                  'NearestCompanion',NearestD, 'X',X, 'Y',Y, 'RefMag',Ref, ...
                  'MagRange',Args.MagRange, 'CompanionRadius',Args.CompanionRadius, ...
                  'CompanionMaxMag',Args.CompanionMaxMag, 'NCompanionCat',sum(IsComp), ...
                  'SelfRadius',Args.SelfRadius, 'NSelfMatchDropped',NSelf);

    if Args.Verbosity > 0
        fprintf(['selectRefSources: %d of %d sources with %g < mag < %g; %d rejected for a companion ', ...
                 'brighter than %g within %g pix; %d kept as reference stars\n'], ...
                Info.NinMag, Nsrc, Args.MagRange(1), Args.MagRange(2), Info.NrejectedByCompanion, ...
                Args.CompanionMaxMag, Args.CompanionRadius, Info.NafterIsolation);
    end
end
