function [PM, Info] = fixedPMFromFit(IF, SrcInd, Args)
% Take the proper motions out of a long joint solution, to be held fixed in a
% shorter one.
% A fit spanning a single observing season constrains a proper motion barely at
% all: the time column of the source design matrix spans two thirds of a year
% while the parameter is a rate per year, so the fitted slope is dominated by
% whatever else drifts through the season, and it absorbs any real motion of
% the source along with it. The motions are far better determined by a fit over
% the whole decade, and this returns them in a form that can be handed back to
% the shorter fit and held.
% The values are relative motions in the registered pixel frame, and their
% overall level is a matter of gauge: the per-epoch transformation of any fit
% is free to absorb a motion shared by the field, together with a slow rotation
% or change of scale, and different fits divide that freedom differently. Only
% the part that varies from source to source is meaningful, so by default the
% shared part is removed here rather than left for the receiving fit to
% rediscover.
% Input  : - An IterFit object holding the joint solution, or the path of a
%            .mat file holding one in a variable named IFsys.
%          - Index of each of the fit's sources in the source list the result
%            should be reported against, normally Info.SrcInd from
%            ml.scripts.runIterDetrendMSc. When the first argument is a file
%            written by that chain the SrcInd stored beside it is used, and
%            this may be left empty.
%          * ...,key,val,...
%            'NsrcTotal' - Length of that source list, which is the source
%                   count of the MatchedSources the fit was made from. Default
%                   is max(SrcInd), which is right only when the last source
%                   survived the cuts.
%            'MinNepoch' - Least number of epochs a source must be measured in
%                   for its motion to be kept. Default is 100.
%            'MaxPM' - Motions faster than this [mas/yr] are dropped as
%                   implausible; the field disperses by about 4.5 mas/yr.
%                   Default is 50.
%            'MaxRstd' - Sources whose residual scatter [mas] exceeds this are
%                   dropped. Default is Inf.
%            'RemoveGauge' - Subtract the constant and the linear dependence on
%                   position that the receiving fit can absorb into its own
%                   per-epoch transformation. Default is true.
%            'GaugeMaxMag' - Only sources brighter than this take part in
%                   determining that shared part. Default is 17.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - A 2 by NsrcTotal matrix of proper motions in pixels per year, in
%            the row order [muX; muY] that ParS(3:4,:) uses, and NaN for every
%            source whose motion was not determined or did not survive the
%            cuts. A fit receiving it leaves those sources' motions free.
%          - A structure recording how many sources each cut removed, the
%            shared part that was subtracted, and the dispersion of what is
%            left.
% Author : ULTRASAT team (Aug 2026)
% <Unknown Bugs>: the gauge is fitted on the bright sources and applied to all,
%   so a real motion shared by the bright sources alone would be removed from
%   the faint ones as well. The two populations sit at different depths in the
%   bulge, so this is not merely hypothetical, but the size of it is well below
%   the scatter of the individual motions.
% Example: PM = ml.util.fixedPMFromFit( ...
%              '~/KMTdata/Results/JointPM1625/IFsys1625_260058_CTIO_BLG41.mat', ...
%              [], 'NsrcTotal',720, 'Verbosity',1);

    arguments
        IF
        SrcInd            = [];
        Args.NsrcTotal    = [];
        Args.MinNepoch    = 100;
        Args.MaxPM        = 50;
        Args.MaxRstd      = Inf;
        Args.RemoveGauge  = true;
        Args.GaugeMaxMag  = 17;
        Args.Verbosity    = 0;
    end

    if ischar(IF) || isstring(IF)
        Loaded = load(char(IF));
        if ~isfield(Loaded,'IFsys')
            error('ml:util:fixedPMFromFit:NoIFsys', ...
                  'The file holds no variable named IFsys');
        end
        if isempty(SrcInd) && isfield(Loaded,'Info') && isstruct(Loaded.Info) ...
                && isfield(Loaded.Info,'SrcInd')
            SrcInd = Loaded.Info.SrcInd;
        end
        IF = Loaded.IFsys;
    end
    if isempty(SrcInd)
        error('ml:util:fixedPMFromFit:NoSrcInd', ...
              'SrcInd is needed to say which source of the common list each fitted source is');
    end
    SrcInd = SrcInd(:).';
    if numel(SrcInd) ~= IF.Nsrc
        error('ml:util:fixedPMFromFit:BadSrcInd', ...
              'SrcInd holds %d entries for a fit of %d sources', numel(SrcInd), IF.Nsrc);
    end
    if isempty(Args.NsrcTotal)
        Args.NsrcTotal = max(SrcInd);
    end

    PMx = IF.ParS(3,:).';
    PMy = IF.ParS(4,:).';
    Mag = IF.medianFieldSource({'MAG_PSF'});
    Mag = Mag(:);
    PosX = IF.ParS(1,:).';
    PosY = IF.ParS(2,:).';
    Nep  = sum(isfinite(IF.Data.X),1).';
    [RstdX, RstdY] = IF.calculateRstd;
    Rstd = sqrt(RstdX(:).^2 + RstdY(:).^2);

    % --- which motions are worth passing on -------------------------------
    Keep = isfinite(PMx) & isfinite(PMy) & isfinite(PosX) & isfinite(PosY);
    Info.NsrcFit       = IF.Nsrc;
    Info.NrejNonFinite = sum(~Keep);
    Drop = Keep & Nep < Args.MinNepoch;
    Info.NrejNepoch    = sum(Drop);
    Keep = Keep & ~Drop;
    Drop = Keep & 400.*hypot(PMx, PMy) > Args.MaxPM;
    Info.NrejMaxPM     = sum(Drop);
    Keep = Keep & ~Drop;
    Drop = Keep & Rstd > Args.MaxRstd;
    Info.NrejMaxRstd   = sum(Drop);
    Keep = Keep & ~Drop;
    Info.Nkept         = sum(Keep);
    report(Args.Verbosity, ['%d of %d fitted sources keep their motion ', ...
           '(%d non-finite, %d too few epochs, %d too fast, %d too noisy)\n'], ...
           Info.Nkept, IF.Nsrc, Info.NrejNonFinite, Info.NrejNepoch, ...
           Info.NrejMaxPM, Info.NrejMaxRstd);

    % --- take out the part the receiving fit can absorb anyway -------------
    % A constant motion of the whole field, and a motion growing linearly with
    % position across it, are together a transformation drifting slowly in
    % time, which is exactly what a free per-epoch affine represents. Leaving
    % them in would ask the receiving fit to find that drift for itself.
    Info.Gauge = struct('Removed',false, 'ParX',[], 'ParY',[], 'Nused',0, ...
                        'AmpX',0, 'AmpY',0);
    if Args.RemoveGauge
        Bright = Keep & isfinite(Mag) & Mag < Args.GaugeMaxMag;
        if sum(Bright) < 20
            report(Args.Verbosity, ...
                'only %d sources brighter than %.1f, leaving the shared part in\n', ...
                sum(Bright), Args.GaugeMaxMag);
        else
            H  = [ones(sum(Bright),1), PosX(Bright)-median(PosX(Bright)), ...
                                       PosY(Bright)-median(PosY(Bright))];
            Bx = H\PMx(Bright);
            By = H\PMy(Bright);
            Hall = [ones(numel(PMx),1), PosX-median(PosX(Bright)), ...
                                        PosY-median(PosY(Bright))];
            GaugeX = Hall*Bx;
            GaugeY = Hall*By;
            PMx = PMx - GaugeX;
            PMy = PMy - GaugeY;
            Info.Gauge = struct('Removed',true, 'ParX',Bx, 'ParY',By, ...
                                'Nused',sum(Bright), ...
                                'AmpX',400.*std(GaugeX(Keep)), ...
                                'AmpY',400.*std(GaugeY(Keep)));
            report(Args.Verbosity, ['shared part removed on %d bright sources: ', ...
                   '%.2f/%.2f mas/yr constant, %.2f/%.2f mas/yr rms over the field\n'], ...
                   sum(Bright), 400.*Bx(1), 400.*By(1), Info.Gauge.AmpX, Info.Gauge.AmpY);
        end
    end

    % --- scatter onto the common source list ------------------------------
    PM = nan(2, Args.NsrcTotal);
    Ok = Keep(:).' & SrcInd >= 1 & SrcInd <= Args.NsrcTotal;
    PM(1, SrcInd(Ok)) = PMx(Ok);
    PM(2, SrcInd(Ok)) = PMy(Ok);

    Info.NsrcTotal  = Args.NsrcTotal;
    Info.Nassigned  = sum(isfinite(PM(1,:)));
    Info.DispMasYr  = [400.*tools.math.stat.rstd(PMx(Keep)), ...
                       400.*tools.math.stat.rstd(PMy(Keep))];
    Info.MedianMag  = median(Mag(Keep), 'omitnan');
    report(Args.Verbosity, ['%d of %d sources of the common list carry a motion, ', ...
           'dispersing by %.2f/%.2f mas/yr\n'], Info.Nassigned, Args.NsrcTotal, ...
           Info.DispMasYr(1), Info.DispMasYr(2));
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
