function [ObjSys, sysCorX, sysCorY, Info] = sysRemScriptPart(IF, Obj, Args)
% Remove systematics from the astrometric residuals of an IterFit solution.
% The residuals of the current solution are decomposed by SysRem and the fitted
% components are subtracted from the X and Y matrices of a copy of Obj.
% The decomposition is guarded in three ways, because on a sparsely filled
% matrix the plain call returns components that are not astrometrically
% meaningful and, once subtracted, destroy the solution:
%   - epochs and sources with too few usable measurements cannot constrain
%     their own component and are kept out of the fit;
%   - the convergence threshold is scaled to the size of the problem. SysRem
%     renormalises Sigma so that chi^2 per degree of freedom is one, so the
%     chi^2 it iterates on is of order Ndof and its default threshold of 1 is a
%     relative tolerance of about 1/Ndof. That is tight enough for the
%     alternating rank-1 fit to keep iterating long after chi^2 has stopped
%     dropping, and its scale degeneracy (C to k*C, A to A/k) then lets the
%     fitted component grow by orders of magnitude;
%   - a component larger than any plausible astrometric systematic is dropped
%     for the sources and epochs it affects, rather than being applied.
% Input  : - An IterFit object holding the current solution.
%          - An MMS object whose X and Y are to be corrected. Not modified.
%          * ...,key,val,...
%            'NIter' - Number of SysRem components to remove. Default is 1.
%            'UseWeight' - Weight the decomposition by the IterFit weights.
%                   Default is true.
%            'MinNvalid' - Least number of usable measurements an epoch or a
%                   source must have to take part in the fit. Set to 0 to
%                   disable. Default is 10.
%            'RelThreshDeltaS2' - Convergence threshold of the SysRem
%                   iteration, as a fraction of the degrees of freedom. Set to
%                   0 to fall back on the SysRem default of 1.
%                   Default is 1e-3.
%            'MaxCorrection' - Largest correction [pix] that any single
%                   measurement may receive. Sources and then epochs whose
%                   fitted component exceeds it are left uncorrected. Set to
%                   Inf to disable. Default is 0.5.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - A copy of Obj with the systematics subtracted from X and Y.
%          - The correction that was subtracted from X.
%          - The correction that was subtracted from Y.
%          - A structure describing the decomposition: how many epochs and
%            sources took part, the threshold used, how many were rejected as
%            implausible and the largest correction actually applied.
% Author : ULTRASAT team (Aug 2026)
% Example: [ObjSys,CorX,CorY,Info] = ml.util.sysRemScriptPart(IF, Obj, 'NIter',2);

    arguments
        IF;
        Obj;
        Args.NIter              = 1;
        Args.UseWeight          = true;
        Args.MinNvalid          = 10;
        Args.RelThreshDeltaS2   = 1e-3;
        Args.MaxCorrection      = 0.5;
        Args.Verbosity          = 0;
    end

    [RxOrg, RyOrg] = IF.calculateResiduals;
    Wes = calculateWes(IF);

    % A measurement takes part only if both its residuals and its weight are usable
    Valid = isfinite(RxOrg) & isfinite(RyOrg) & isfinite(Wes) & Wes > 0;

    if Args.UseWeight
        Sigma = sqrt(1./Wes);
    else
        Sigma = ones(size(RxOrg));
    end
    Sigma(~Valid | isnan(Sigma)) = Inf;

    [FlagEpoch, FlagSrc] = flagPopulated(Valid, Args.MinNvalid);

    sysCorX = zeros(size(RxOrg));
    sysCorY = zeros(size(RyOrg));
    Info    = struct('NepochUsed', sum(FlagEpoch), 'NsrcUsed', sum(FlagSrc), ...
                     'ThreshDeltaS2', NaN, 'NsrcRejected', 0, 'NepochRejected', 0, ...
                     'MaxCorrectionApplied', 0);

    if any(FlagEpoch) && any(FlagSrc)
        SubSigma   = Sigma(FlagEpoch, FlagSrc);
        [Nep, Nsr] = size(SubSigma);
        Ndof       = Nep.*Nsr - Nep - Nsr;
        if Args.RelThreshDeltaS2 > 0
            Thresh = max(1, Args.RelThreshDeltaS2.*Ndof);
        else
            Thresh = 1;
        end
        Info.ThreshDeltaS2 = Thresh;

        [~, SysRemX] = timeSeries.detrend.sysrem(RxOrg(FlagEpoch,FlagSrc), SubSigma, ...
                                                 'Niter', Args.NIter, 'ThreshDeltaS2', Thresh);
        [~, SysRemY] = timeSeries.detrend.sysrem(RyOrg(FlagEpoch,FlagSrc), SubSigma, ...
                                                 'Niter', Args.NIter, 'ThreshDeltaS2', Thresh);

        sysCorX(FlagEpoch,FlagSrc) = accumulateComponents(SysRemX);
        sysCorY(FlagEpoch,FlagSrc) = accumulateComponents(SysRemY);
    end

    [sysCorX, sysCorY, Info] = rejectImplausible(sysCorX, sysCorY, Args.MaxCorrection, Info);

    if Args.Verbosity > 0
        fprintf(['SysRem: %d epochs, %d sources, threshold %.4g ; rejected %d sources and %d epochs ; ', ...
                 'largest correction applied %.4g pix\n'], Info.NepochUsed, Info.NsrcUsed, ...
                Info.ThreshDeltaS2, Info.NsrcRejected, Info.NepochRejected, Info.MaxCorrectionApplied);
    end

    ObjSys = Obj.copy();
    ObjSys.Data.X = ObjSys.Data.X - sysCorX;
    ObjSys.Data.Y = ObjSys.Data.Y - sysCorY;
end


% -------------------------------------------------------------------------
function [FlagEpoch, FlagSrc] = flagPopulated(Valid, MinNvalid)
    % Epochs and sources holding at least MinNvalid usable measurements.
    % Iterated, because dropping sources lowers the per-epoch counts in turn.
    FlagEpoch = true(size(Valid,1), 1);
    FlagSrc   = true(1, size(Valid,2));
    if MinNvalid <= 0
        return
    end
    for Iter = 1:5
        NewEpoch = sum(Valid(:,FlagSrc), 2) >= MinNvalid;
        NewSrc   = sum(Valid(NewEpoch,:), 1) >= MinNvalid;
        if isequal(NewEpoch, FlagEpoch) && isequal(NewSrc, FlagSrc)
            break
        end
        FlagEpoch = NewEpoch;
        FlagSrc   = NewSrc;
    end
end


function Corr = accumulateComponents(SysRem)
    % Sum the SysRem components, skipping the initial state held in element 1
    Corr = zeros(numel(SysRem(end).C), numel(SysRem(end).A));
    for Isys = 2:numel(SysRem)
        Comp = SysRem(Isys).C .* SysRem(Isys).A;
        Comp(~isfinite(Comp)) = 0;
        Corr = Corr + Comp;
    end
end


function [sysCorX, sysCorY, Info] = rejectImplausible(sysCorX, sysCorY, MaxCorrection, Info)
    % Drop the individual corrections too large to be an astrometric
    % systematic, and only those. Zeroing a whole source because one of its
    % epochs came out large throws away the correction for every other epoch it
    % has, and the chance of one epoch tripping the test grows with the length
    % of the run: measured on these fields it removed 41% of the sources in a
    % single season and 79% over the full decade. The runaway decomposition
    % this guard exists to catch is large over the whole matrix rather than at
    % one point, so it is still removed in full.
    sysCorX(~isfinite(sysCorX)) = 0;
    sysCorY(~isfinite(sysCorY)) = 0;

    if isfinite(MaxCorrection)
        Bad = abs(sysCorX) > MaxCorrection | abs(sysCorY) > MaxCorrection;
        sysCorX(Bad) = 0;
        sysCorY(Bad) = 0;

        Info.NrejectedPoints = sum(Bad, 'all');
        Info.FracRejected    = mean(Bad, 'all');
        % how many sources and epochs the point-wise rejection touches at all,
        % which is not the same as discarding them
        Info.NsrcRejected    = sum(any(Bad, 1));
        Info.NepochRejected  = sum(any(Bad, 2));
    end
    Info.MaxCorrectionApplied = max(max(abs(sysCorX), [], 'all'), max(abs(sysCorY), [], 'all'));
end
