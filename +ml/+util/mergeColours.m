function [Colour, Info] = mergeColours(MSI, MSV, Args)
% Merge a KMT V-band colour with the OGLE colour into one consistent scale.
% The colour feeds the differential refraction bins, and refraction acts in the
% band the astrometry was measured in, so a KMT colour is the right variable
% where one exists. It does not exist everywhere: the V frames are far rarer
% and shallower than the I frames. OGLE reaches sources KMT V does not, and the
% reverse is also true, so the two together cover far more than either alone.
% They cannot simply be pooled. The KMT colour is instrumental and the OGLE one
% calibrated, so pooling them unscaled would place sources from one catalogue
% systematically in different bins from sources of the other, which is itself a
% colour dependent systematic. The OGLE colour is therefore carried onto the
% KMT scale by a straight line fitted on the sources that have both.
% Input  : - The I-band MatchedSources, whose SrcData holds the OGLE
%            magnitudes and whose sources define the output ordering.
%          - The V-band MatchedSources, or the path of a .mat holding one in a
%            variable named MSc.
%          * ...,key,val,...
%            'ColourFields' - SrcData fields {Blue, Red} of the OGLE colour.
%                   Default is {'V_ogle','I_ogle'}.
%            'MaxValidMag' - OGLE magnitudes at or above this are missing-value
%                   sentinels. A single one of these, left in, has enough
%                   leverage to flatten the fitted relation to nothing.
%                   Default is 99.
%            'ZP0' - Zero point of the instrumental magnitudes. Default is 25.
%            'MatchRadius' - Largest separation of a V source from an I source
%                   for the two to be the same star [pix]. Default is 1.5.
%            'Offset' - Offset [dx dy] of the V frame from the I frame [pix].
%                   Fitted from the data when empty. Default is [].
%            'MinOverlap' - Fewest sources with both colours for the transform
%                   to be fitted. Default is 20.
%            'ColourRange' - Colours outside this are not used. Default is
%                   [-1 5].
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - One colour per source of the I-band object, on the KMT scale, NaN
%            where neither catalogue reaches it.
%          - A structure holding the fitted transform, its scatter, how many
%            sources each catalogue supplied, and a per-source provenance of
%            'kmt', 'ogle' or 'none'.
% Author : ULTRASAT team (Sep 2026)
% <Unknown Bugs>: the transform is fitted on the sources bright enough for both
%   surveys and applied to sources fainter than any of them, so it is an
%   extrapolation. Measured on BLG01 its slope is 0.956 +- 0.092 with 0.059 mag
%   of scatter, near enough to unity that the extrapolation is mild, but a
%   colour dependent error in either survey at the faint end would pass through.
% Example: C = ml.util.mergeColours(MSc, '~/KMTdata/ResultsV/KMT_260058_BLG01_Vband_MSc.mat');
%          Obj = ml.util.mmsFromMatchedSources(MSc, 'ColourFields',{'Colour'});

    arguments
        MSI
        MSV
        Args.ColourFields = {'V_ogle','I_ogle'};
        Args.MaxValidMag  = 99;
        Args.ZP0          = 25;
        Args.MatchRadius  = 1.5;
        Args.Offset       = [];
        Args.MinOverlap   = 20;
        Args.ColourRange  = [-1 5];
        Args.Verbosity    = 0;
    end

    if ischar(MSV) || isstring(MSV)
        Loaded = load(char(MSV));
        Vars   = fieldnames(Loaded);
        IsMS   = cellfun(@(V) isa(Loaded.(V),'MatchedSources'), Vars);
        MSV    = Loaded.(Vars{find(IsMS,1)});
    end

    % --- the OGLE colour, sentinels removed --------------------------------
    Blue = MSI.SrcData.(Args.ColourFields{1})(:);
    Red  = MSI.SrcData.(Args.ColourFields{2})(:);
    Blue(Blue >= Args.MaxValidMag) = NaN;
    Red(Red   >= Args.MaxValidMag) = NaN;
    OgleCol = Blue - Red;
    OgleCol(OgleCol < Args.ColourRange(1) | OgleCol > Args.ColourRange(2)) = NaN;

    % --- the KMT colour, instrumental in both bands ------------------------
    % The pipeline leaves MAG_PSF empty, so the magnitudes come from the flux.
    Mv = Args.ZP0 - 2.5.*log10(median(MSV.Data.FLUX_PSF, 1, 'omitnan')).';
    Mi = Args.ZP0 - 2.5.*log10(median(MSI.Data.FLUX_PSF, 1, 'omitnan')).';
    Xv = median(MSV.Data.X, 1, 'omitnan').';  Yv = median(MSV.Data.Y, 1, 'omitnan').';
    Xi = median(MSI.Data.X, 1, 'omitnan').';  Yi = median(MSI.Data.Y, 1, 'omitnan').';

    % The two bands are registered separately, so their frames differ by a
    % shift. Take it from the median of the nearest-neighbour offsets.
    if isempty(Args.Offset)
        [~, J0] = pdistNearest(Xv, Yv, Xi, Yi, 0, 0);
        Args.Offset = [median(Xi(J0)-Xv, 'omitnan'), median(Yi(J0)-Yv, 'omitnan')];
    end
    [D, J] = pdistNearest(Xv, Yv, Xi, Yi, Args.Offset(1), Args.Offset(2));
    Matched = D < Args.MatchRadius;

    KmtCol = nan(MSI.Nsrc, 1);
    KmtCol(J(Matched)) = Mv(Matched) - Mi(J(Matched));
    KmtCol(KmtCol < Args.ColourRange(1) | KmtCol > Args.ColourRange(2)) = NaN;

    % --- carry OGLE onto the KMT scale -------------------------------------
    Both = isfinite(KmtCol) & isfinite(OgleCol);
    Info = struct('Offset',Args.Offset, 'NmatchedV',sum(Matched), 'NsrcV',MSV.Nsrc, ...
                  'Noverlap',sum(Both), 'Transform',[NaN NaN], 'Scatter',NaN, ...
                  'Slope',NaN, 'SlopeErr',NaN);
    if sum(Both) < Args.MinOverlap
        error('ml:util:mergeColours:TooFewOverlap', ...
              'Only %d sources carry both colours, %d are needed to fit the transform', ...
              sum(Both), Args.MinOverlap);
    end
    B = robustfit(OgleCol(Both), KmtCol(Both));
    Resid = KmtCol(Both) - (B(1) + B(2).*OgleCol(Both));
    Info.Transform = B(:).';
    Info.Slope     = B(2);
    Info.SlopeErr  = std(Resid)./(sqrt(sum(Both)).*std(OgleCol(Both)));
    Info.Scatter   = tools.math.stat.rstd(Resid);
    OgleOnKmt = B(1) + B(2).*OgleCol;

    % --- KMT where it exists, OGLE carried over elsewhere -------------------
    Colour = KmtCol;
    FromOgle = ~isfinite(Colour) & isfinite(OgleOnKmt);
    Colour(FromOgle) = OgleOnKmt(FromOgle);

    Prov = repmat({'none'}, MSI.Nsrc, 1);
    Prov(isfinite(KmtCol)) = {'kmt'};
    Prov(FromOgle)         = {'ogle'};
    Info.Provenance = Prov;
    Info.Nkmt   = sum(isfinite(KmtCol));
    Info.Nogle  = sum(FromOgle);
    Info.Nvalid = sum(isfinite(Colour));
    Info.FracValid = Info.Nvalid./MSI.Nsrc;
    report(Args.Verbosity, ['frame offset %+.2f,%+.2f pix; %d of %d V sources matched\n', ...
           'transform onto the KMT scale: slope %.3f +- %.3f, scatter %.3f mag on %d sources\n', ...
           'colour for %d of %d sources (%.1f%%): %d from KMT, %d from OGLE\n'], ...
           Args.Offset(1), Args.Offset(2), sum(Matched), MSV.Nsrc, ...
           Info.Slope, Info.SlopeErr, Info.Scatter, Info.Noverlap, ...
           Info.Nvalid, MSI.Nsrc, 100.*Info.FracValid, Info.Nkmt, Info.Nogle);
end


function [D, J] = pdistNearest(Xv, Yv, Xi, Yi, Dx, Dy)
    % Nearest I-band source to every V-band source, after a shift
    Nv = numel(Xv);
    D  = nan(Nv,1);
    J  = nan(Nv,1);
    for K = 1:Nv
        Dist = hypot(Xi - (Xv(K)+Dx), Yi - (Yv(K)+Dy));
        [D(K), J(K)] = min(Dist);
    end
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
