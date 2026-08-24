function [IFsys, Obj, IFsysB, Info] = runIterDetrendMSc(MS, Args)
% Run the Lizoo iterative astrometric detrending on a KMT_pipelineI MatchedSources.
% This is the ml.scripts.runOverFields analogue for input that has already
% been matched by KMT_pipelineI: the per-epoch AstCat files and the msMatch
% matching step are replaced by ml.util.mmsFromMatchedSources, and the rest of
% the chain (detrend, SysRem, final detrend, per-source tables) is identical.
% Input  : - A MatchedSources object as produced by KMT_pipelineI, or a
%            char/string with the path of a .mat file containing one.
%          * ...,key,val,...
%            'mmsFromMatchedSourcesArgs' - Cell array of arguments passed on
%                   to ml.util.mmsFromMatchedSources. Default is {}.
%            'CelestialCoo' - Field coordinates [RA, Dec] in radians. If
%                   empty, the value returned by the converter is used.
%                   Default is [].
%            'HALat' - Fit the airmass/parallactic-angle terms. Default is true.
%            'AnnualEffect' - Fit the annual terms. Default is true.
%            'UseWeights' - Use weights in the iterations. Default is true.
%            'PixPhase' - Apply the pixel-phase correction. Requires the
%                   per-epoch shifts to have been stored by the pipeline.
%                   Default is false.
%            'NiterNoWeightsBeforeSys' - Unweighted iterations of the first
%                   pass. Default is 2.
%            'NiterWeightsBeforeSys' - Weighted iterations of the first pass.
%                   Default is 10.
%            'NiterWeightsAfterSys' - Weighted iterations of the final pass.
%                   Default is 4.
%            'ChromaicHighOrder' - Use the high-order airmass terms of
%                   @IterFit/generateHALatDesignMat. Their columns are powers
%                   of secz*sin(pa) and secz*cos(pa) and so are strongly
%                   collinear; on a sparse matrix the normal matrix becomes
%                   ill-conditioned enough that the bicg solve inside runIter
%                   stops converging. Set false for the well-conditioned
%                   linear form. Default is true, as in the current pipeline.
%            'RunSysRem' - Run the SysRem step between the two passes.
%                   Default is true.
%            'MaxSysRemShift' - Largest astrometric correction [pix] that the
%                   SysRem step is allowed to apply. A larger correction, or
%                   any non-finite one, means the sysrem decomposition did not
%                   behave, and the step is rejected rather than applied.
%                   Default is 1.
%            'NIterSysRem' - Number of SysRem iterations. Default is 2.
%            'PerSourcesTargetPath' - Directory for the per-source csv tables
%                   written by ml.scripts.IterFitToPerSourceFormat. No tables
%                   are written when empty. Default is ''.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - The IterFit object of the final pass.
%          - The MMS object handed to the final pass.
%          - The IterFit object of the first pass, before SysRem.
%          - A structure describing the conversion and the run, including
%            SrcInd/EpochInd back into the input object and, when tables were
%            written, their path in Info.CatsPath.
% Author : ULTRASAT team (Aug 2026)
% Example: [IF,Obj,IFB,Info] = ml.scripts.runIterDetrendMSc(...
%              '/bigdata3/projects/KMTdata/Results16_26/KMT_260058_BLG41_MSc.mat');

    arguments
        MS
        Args.mmsFromMatchedSourcesArgs = {};
        Args.CelestialCoo              = [];
        Args.HALat                     = true;
        Args.AnnualEffect              = true;
        Args.UseWeights                = true;
        Args.PixPhase                  = false;
        Args.NiterNoWeightsBeforeSys   = 2;
        Args.NiterWeightsBeforeSys     = 10;
        Args.NiterWeightsAfterSys      = 4;
        Args.ChromaicHighOrder         = true;
        Args.RunSysRem                 = true;
        Args.MaxSysRemShift            = 1;
        Args.NIterSysRem               = 2;
        Args.PerSourcesTargetPath      = '';
        Args.Verbosity                 = 0;
    end

    % --- input conversion ---------------------------------------------------
    if isa(MS,'MMS')
        Obj          = MS;
        CelestialCoo = Args.CelestialCoo;
        Info         = struct('SrcInd',(1:MS.Nsrc)', 'EpochInd',(1:MS.Nepoch)', ...
                              'PixPhaseAvailable', isfield(MS.Data,'Xphase'));
    else
        [Obj, CelestialCoo, Info] = ml.util.mmsFromMatchedSources(MS, ...
            'Verbosity', Args.Verbosity, Args.mmsFromMatchedSourcesArgs{:});
    end
    if ~isempty(Args.CelestialCoo)
        CelestialCoo = Args.CelestialCoo;
    end

    if Args.PixPhase && ~Info.PixPhaseAvailable
        error('runIterDetrendMSc:NoPixPhase', ...
              ['PixPhase was requested but the pixel phase is not recoverable: X and Y are the ', ...
               'registered coordinates and the per-epoch shifts were not stored by the pipeline.']);
    end

    % --- first pass, before SysRem -----------------------------------------
    report(Args.Verbosity, 'Detrending pass 1 (%d epochs, %d sources)\n', Obj.Nepoch, Obj.Nsrc);
    [IFsysB, MMSsysB] = ml.scripts.runIterDetrend(Obj.copy(), ...
        'CelestialCoo',   CelestialCoo, ...
        'HALat',          Args.HALat, ...
        'UseWeights',     Args.UseWeights, ...
        'Plx',            false, ...
        'PixPhase',       Args.PixPhase, ...
        'AnnualEffect',   Args.AnnualEffect, ...
        'ChromaicHighOrder', Args.ChromaicHighOrder, ...
        'NiterWeights',   Args.NiterWeightsBeforeSys, ...
        'NiterNoWeights', Args.NiterNoWeightsBeforeSys);

    % --- SysRem on the astrometric residuals --------------------------------
    % ml.util.sysRemScriptPart builds Sigma as sqrt(1./Wes), which is Inf at
    % every masked point. On a sparse matrix an epoch or a source can be left
    % with no unmasked point at all, and the sysrem normal sums then divide by
    % zero; the resulting correction is applied to X and Y regardless. The
    % correction is therefore checked before it is accepted.
    Info.SysRemApplied      = false;
    Info.SysRemRejectReason = '';
    ObjSysAfter             = MMSsysB;
    if Args.RunSysRem
        report(Args.Verbosity, 'SysRem (%d iterations)\n', Args.NIterSysRem);
        try
            [ObjCandidate, SysCorX, SysCorY] = ml.util.sysRemScriptPart(IFsysB, MMSsysB, ...
                'UseWeight', Args.UseWeights, 'NIter', Args.NIterSysRem);
            Info.SysRemRejectReason = checkSysRem(MMSsysB, ObjCandidate, SysCorX, SysCorY, Args.MaxSysRemShift);
            if isempty(Info.SysRemRejectReason)
                ObjSysAfter        = ObjCandidate;
                Info.SysCorX       = SysCorX;
                Info.SysCorY       = SysCorY;
                Info.SysRemApplied = true;
            else
                fprintf('SysRem correction rejected and skipped: %s\n', Info.SysRemRejectReason);
            end
        catch ME
            Info.SysRemRejectReason = ME.message;
            fprintf('SysRem failed and was skipped: %s\n', ME.message);
        end
    end

    % Sources left without a position after SysRem are rare but do occur.
    Xguess  = median(ObjSysAfter.Data.X, 'omitnan').';
    Yguess  = median(ObjSysAfter.Data.Y, 'omitnan').';
    FlagNan = ~(isnan(Xguess) | isnan(Yguess));

    IFinit = IFsysB.copy();
    if ~all(FlagNan)
        % Drop those sources from the object and from the parameters of the
        % first-pass solution that is being reused, so the two stay aligned.
        ObjSysAfter.Data = ml.util.flag_struct_field(ObjSysAfter.Data, FlagNan, 'FlagByCol', true);
        IFinit.ParS      = IFinit.ParS(:, FlagNan);
        Info.SrcInd      = Info.SrcInd(FlagNan);
        if isfield(Info,'SrcData') && isstruct(Info.SrcData)
            SrcFields = fieldnames(Info.SrcData);
            for Ifield = 1:numel(SrcFields)
                Info.SrcData.(SrcFields{Ifield}) = Info.SrcData.(SrcFields{Ifield})(FlagNan);
            end
        end
        report(Args.Verbosity, 'Dropped %d sources left without a position by SysRem\n', sum(~FlagNan));
    end
    Info.NsrcFinal = numel(Info.SrcInd);

    % --- final pass ---------------------------------------------------------
    report(Args.Verbosity, 'Detrending pass 2 (final)\n');
    [IFsys, Obj] = ml.scripts.runIterDetrend(ObjSysAfter, ...
        'IF',             IFinit, ...
        'CelestialCoo',   CelestialCoo, ...
        'HALat',          Args.HALat, ...
        'UseWeights',     Args.UseWeights, ...
        'Plx',            false, ...
        'PixPhase',       Args.PixPhase, ...
        'AnnualEffect',   Args.AnnualEffect, ...
        'ChromaicHighOrder', Args.ChromaicHighOrder, ...
        'NiterWeights',   Args.NiterWeightsAfterSys, ...
        'NiterNoWeights', Args.NiterNoWeightsBeforeSys, ...
        'FinalStep',      true);

    % --- per-source tables --------------------------------------------------
    Info.CatsPath = '';
    if ~isempty(Args.PerSourcesTargetPath)
        try
            [Info.CatsPath, Info.TableSrcInd] = ml.scripts.IterFitToPerSourceFormat(IFsys, Args.PerSourcesTargetPath);
            report(Args.Verbosity, 'Wrote per-source tables to %s\n', Info.CatsPath);
        catch ME
            fprintf('Failed to write the per-source tables: %s\n', ME.message);
        end
    end
end


function Reason = checkSysRem(ObjBefore, ObjAfter, SysCorX, SysCorY, MaxShift)
    % Empty when the SysRem correction is safe to apply, otherwise the reason
    Reason = '';
    if ~all(isfinite(SysCorX(:))) || ~all(isfinite(SysCorY(:)))
        Reason = 'the correction contains non-finite values';
        return
    end
    MaxCorr = max(max(abs(SysCorX(:))), max(abs(SysCorY(:))));
    if MaxCorr > MaxShift
        Reason = sprintf('the correction reaches %.3g pix, above MaxSysRemShift = %.3g', MaxCorr, MaxShift);
        return
    end
    LostX = isfinite(ObjBefore.Data.X) & ~isfinite(ObjAfter.Data.X);
    LostY = isfinite(ObjBefore.Data.Y) & ~isfinite(ObjAfter.Data.Y);
    if any(LostX(:)) || any(LostY(:))
        Reason = sprintf('%d measurements became non-finite', sum(LostX(:)) + sum(LostY(:)));
    end
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
