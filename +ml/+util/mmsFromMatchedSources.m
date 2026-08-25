function [Obj, CelestialCoo, Info] = mmsFromMatchedSources(MS, Args)
% Build an MMS object ready for IterFit detrending from a MatchedSources object.
% Input  : - A MatchedSources object, as produced by KMT_pipelineI, or a
%            char/string holding the path of a .mat file containing such an
%            object. The input object is never modified.
%          * ...,key,val,...
%            'RA' - J2000 RA [deg] of the cut-out centre, used to derive the
%                   per-epoch airmass and parallactic angle. Overridden by
%                   MS.UserData.RA when present. Default is the KB260058
%                   BLG41 cut-out centre.
%            'Dec' - J2000 Dec [deg] of the cut-out centre. Overridden by
%                   MS.UserData.Dec when present. Default is the KB260058
%                   BLG41 cut-out centre.
%            'GeoPos' - Observatory geodetic position [Lon, Lat] in deg,
%                   East longitude positive. Default is CTIO (KMTNet-South).
%            'KeepFields' - Cell array of MS.Data fields to carry over, or
%                   'all'. Each field costs Nepoch*Nsrc doubles and is copied
%                   several times downstream, so the default is the lean set
%                   actually used by IterFit. Default is
%                   {'X','Y','FLUX_PSF','MAGERR_PSF','PSF_CHI2DOF','SN'}.
%            'ColNameFlux' - Flux field from which the magnitude is built.
%                   Default is 'FLUX_PSF'.
%            'ColNameMag' - Name of the magnitude field to create.
%                   Default is 'MAG_PSF'.
%            'ZP0' - Zero point of the instrumental magnitude, i.e.
%                   Mag = ZP0 - 2.5*log10(Flux). Default is 25.
%            'MaxSecz' - Epochs with secz above this value are removed.
%                   Set to Inf to disable. Default is 1.6.
%            'MinNdetSrc' - Sources detected in fewer epochs than this are
%                   removed. Default is 200.
%            'MinDetFrac' - Sources detected in a smaller fraction of the
%                   surviving epochs than this are removed. Acts together with
%                   MinNdetSrc, whichever is stricter, so that the cut follows
%                   the number of epochs instead of being tied to one dataset.
%                   The detrending needs a reasonably filled matrix: the
%                   sysrem step divides by per-epoch and per-source weight
%                   sums, and the airmass normal matrix becomes ill
%                   conditioned when the matrix is mostly empty.
%                   Default is 0.3.
%            'MaxSrcStdXY' - Sources whose X or Y scatter exceeds this value
%                   [pix] are removed. Together with the detection cut this is
%                   the FlagGoodSrc criterion of KMT_pipelineI, recomputed
%                   here from the data for objects whose UserData does not
%                   carry it. Set to Inf to disable. Default is 1.
%            'MinNsrcEpoch' - Epochs with fewer detected sources than this
%                   are removed. Default is 50.
%            'RemoveDuplicateSrc' - Remove sources that are exact duplicates
%                   of another source, i.e. share its median X and Y to the
%                   last bit. The unified source list of KMT_pipelineI does
%                   contain a few of these, and a single repeated median
%                   magnitude is enough to break the outlier rejection that
%                   @IterFit/calculateWes relies on. Default is true.
%            'RequireFiniteMag' - Remove sources whose median magnitude is not
%                   finite, for the same reason. Default is true.
%            'ApplyPipelineFlags' - Apply the FlagGoodEpoch/FlagGoodSrc
%                   selections that KMT_pipelineI recorded in UserData, when
%                   they are present and have not already been applied to the
%                   object. Default is true.
%            'MinNormPeakCorr' - Reject epochs whose registration
%                   cross-correlation peak, as recorded in UserData, is below
%                   this value. No cut when empty. Default is [].
%            'ApplyRefZP' - Fit and apply a per-epoch zero point against the
%                   reference magnitude. Default is true.
%            'RefMagField' - SrcData field holding the reference magnitude
%                   used to anchor the zero point. Default is 'I_ogle'.
%            'fitRefZPArgs' - Cell array of arguments passed to MMS/fitRefZP.
%                   Default is {'ZPFun',@median,'ZPFunArgs',{'omitnan'}}.
%            'ApplySysRemPhot' - Run MMS/SysRemPhotometry on the magnitudes.
%                   Default is true.
%            'SysRemPhotometryArgs' - Cell array of arguments passed to
%                   MMS/SysRemPhotometry.
%                   Default is {'ThreshDeltaS2',1,'Niter',10}.
%            'ColourFields' - SrcData fields {Blue, Red} whose difference is
%                   the colour. Default is {'V_ogle','I_ogle'}.
%            'ColourMode' - How to treat sources without a valid colour:
%                   'fill'     - assign the median colour of the field.
%                   'fillcmd'  - assign the colour predicted by a robust
%                                colour-magnitude relation fitted to the
%                                sources that do have a colour.
%                   'restrict' - keep only sources with a valid colour.
%                   'ownbin'   - keep all sources, but place the colourless
%                                ones below every valid colour so that they
%                                occupy the lowest bins instead of being mixed
%                                in with real colours. The bins IterFit builds
%                                are equal-population quantiles, so they will
%                                span more than one bin.
%                   Default is 'ownbin'. Measured on BLG41 at 4330 epochs,
%                   where 87% of the sources have a real colour and so the
%                   fill barely invents anything, the four modes give a bright
%                   source rstd of X/Y 8.70/7.85 ('fill'), 8.91/8.02
%                   ('fillcmd'), 7.41/7.80 ('restrict', but 515 sources rather
%                   than 592) and 7.43/7.26 ('ownbin'). 'ownbin' keeps every
%                   source and is the best or close to it throughout, most
%                   clearly between 16th and 17th magnitude, 11.8 mas against
%                   15.4 for 'fill'.
%            'ColourRange' - Colours outside this range are treated as
%                   invalid. Default is [0 5].
%            'MaxValidRefMag' - Reference magnitudes at or above this value
%                   are missing-value sentinels. Default is 99.
%            'NColourBins' - Number of colour bins IterFit will use; only
%                   needed so that tied colours can be separated before
%                   @IterFit/generateBins takes quantiles of them.
%                   Default is 6.
%            'CentreColour' - Subtract the median colour, as
%                   ml.util.loadAstCatMatch does. Default is true.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - An MMS object whose Data holds X, Y, the magnitude, C, secz, pa,
%            ha and alt (plus Xphase/Yphase when the per-epoch shifts were
%            stored by the pipeline), ready for ml.scripts.runIterDetrend.
%          - Celestial coordinates [RA, Dec] in radians, to be passed as the
%            'CelestialCoo' argument of ml.scripts.runIterDetrend.
%          - A structure describing the conversion: the surviving epoch and
%            source indices into the input object, the colour statistics and
%            which optional steps were actually applied.
% Author : ULTRASAT team (Aug 2026)
% Example: [Obj,Coo,Info] = ml.util.mmsFromMatchedSources(...
%              '/bigdata3/projects/KMTdata/Results16_26/KMT_260058_BLG41_MSc.mat');

    arguments
        MS
        Args.RA                       = celestial.coo.convertdms('17:54:16.84','gH','d');
        Args.Dec                      = celestial.coo.convertdms('-31:08:44.3','gD','d');
        Args.GeoPos                   = [-70.80399722, -30.16717778];  % CTIO [Lon Lat] deg
        Args.KeepFields               = {'X','Y','FLUX_PSF','MAGERR_PSF','PSF_CHI2DOF','SN'};
        Args.ColNameFlux              = 'FLUX_PSF';
        Args.ColNameMag               = 'MAG_PSF';
        Args.ZP0                      = 25;
        Args.MaxSecz                  = 1.6;
        Args.MinNdetSrc               = 200;
        Args.MinDetFrac               = 0.3;
        Args.MaxSrcStdXY              = 1;
        Args.MinNsrcEpoch             = 50;
        Args.RemoveDuplicateSrc       = true;
        Args.RequireFiniteMag         = true;
        Args.ApplyPipelineFlags       = true;
        Args.MinNormPeakCorr          = [];
        Args.ApplyRefZP               = true;
        Args.RefMagField              = 'I_ogle';
        Args.fitRefZPArgs             = {'ZPFun',@median,'ZPFunArgs',{'omitnan'}};
        Args.ApplySysRemPhot          = true;
        Args.SysRemPhotometryArgs     = {'ThreshDeltaS2',1,'Niter',10};
        Args.ColourFields             = {'V_ogle','I_ogle'};
        Args.ColourMode               = 'ownbin';
        Args.ColourRange              = [0 5];
        Args.MaxValidRefMag           = 99;
        Args.NColourBins              = 6;
        Args.CentreColour             = true;
        Args.Verbosity                = 0;
    end

    Info = struct();

    % --- read the MatchedSources object ------------------------------------
    if ischar(MS) || isstring(MS)
        Info.SourceFile = char(MS);
        MS = readMatchedSourcesFile(Info.SourceFile);
    else
        Info.SourceFile = '';
    end
    if ~isa(MS,'MatchedSources')
        error('mmsFromMatchedSources:BadInput','First input must be a MatchedSources object or a path to a .mat file holding one');
    end

    Nepoch0 = MS.Nepoch;
    Nsrc0   = MS.Nsrc;
    Info.Nepoch0 = Nepoch0;
    Info.Nsrc0   = Nsrc0;

    % --- pointing and site, preferring what the pipeline stored ------------
    UserData = MS.UserData;
    if isstruct(UserData)
        if isfield(UserData,'RA') && ~isempty(UserData.RA)
            Args.RA = UserData.RA;
        end
        if isfield(UserData,'Dec') && ~isempty(UserData.Dec)
            Args.Dec = UserData.Dec;
        end
        if isfield(UserData,'GeoPos') && ~isempty(UserData.GeoPos)
            Args.GeoPos = UserData.GeoPos;
        end
    end
    Info.RA     = Args.RA;
    Info.Dec    = Args.Dec;
    Info.GeoPos = Args.GeoPos;
    CelestialCoo = [Args.RA, Args.Dec]./180.*pi;

    % --- assemble a lean MMS ------------------------------------------------
    Obj    = MMS;
    Obj.JD = MS.JD(:);
    if ischar(Args.KeepFields) || isstring(Args.KeepFields)
        Fields = fieldnames(MS.Data);
    else
        Fields = intersect(Args.KeepFields, fieldnames(MS.Data), 'stable');
        Missing = setdiff(Args.KeepFields, fieldnames(MS.Data));
        if ~isempty(Missing)
            error('mmsFromMatchedSources:MissingField','Requested field(s) absent from the input Data: %s', strjoin(Missing,', '));
        end
    end
    for Ifield = 1:numel(Fields)
        Obj.Data.(Fields{Ifield}) = MS.Data.(Fields{Ifield});
    end
    Required = [{'X','Y'}, {Args.ColNameFlux}];
    Absent   = Required(~isfield(Obj.Data, Required));
    if ~isempty(Absent)
        error('mmsFromMatchedSources:MissingRequired','Field(s) %s must be among KeepFields', strjoin(Absent,', '));
    end

    % Source-level quantities are carried separately: they are indexed by
    % source only, so they must not go into Data (whose fields are indexed by
    % [epoch, source]) before they have been expanded.
    SrcData = MS.SrcData;

    % --- chronological order ------------------------------------------------
    % The epoch order follows the file listing used by the pipeline, which is
    % not strictly chronological. @IterFit/calculateWes and calculateRstd pass
    % the JD to isoutlier as 'SamplePoints', which requires sorted points.
    EpochMap = (1:Nepoch0).';
    Info.EpochSorted = ~issorted(Obj.JD);
    if Info.EpochSorted
        [~, IsortJD] = sort(Obj.JD);
        Obj.selectByEpoch(IsortJD, 'CreateNewObj', false);
        EpochMap = EpochMap(IsortJD);
    end

    % --- per-epoch observing geometry, derived from the JD -----------------
    % The saved MatchedSources holds no header metadata, but secz and the
    % parallactic angle follow from the mid-exposure JD, the cut-out centre
    % and the site. pa uses the atan convention of @ImRed/populateMetaData.
    [~, Alt, ~, PA, HA] = celestial.coo.radec2azalt(Obj.JD, Args.RA, Args.Dec, ...
                                                    'GeoCoo', Args.GeoPos, ...
                                                    'InUnits', 'deg', ...
                                                    'OutUnits', 'rad');
    Secz = 1./sin(Alt);
    HA   = mod(HA + pi, 2.*pi) - pi;

    % --- epoch selection ----------------------------------------------------
    FlagEpoch = Secz < Args.MaxSecz & isfinite(Secz);
    [~, ~, IndUnique] = unique(Obj.JD, 'stable');
    CountJD   = accumarray(IndUnique, 1);
    FlagEpoch = FlagEpoch & CountJD(IndUnique) == 1;
    Info.NepochRejectedSecz = sum(~FlagEpoch);

    % Epoch quality recorded by the pipeline, when it is available
    Info.NepochRejectedPipeline = 0;
    [NormPeakCorr, HasNPC] = alignEpochVector(getUserField(UserData,'NormPeakCorr'), UserData, Nepoch0);
    if HasNPC && ~isempty(Args.MinNormPeakCorr)
        FlagEpoch = FlagEpoch & NormPeakCorr >= Args.MinNormPeakCorr;
    end
    if Args.ApplyPipelineFlags
        FlagGoodEpochUD = getUserField(UserData,'FlagGoodEpoch');
        if numel(FlagGoodEpochUD) == Nepoch0
            Before    = sum(FlagEpoch);
            FlagEpoch = FlagEpoch & logical(FlagGoodEpochUD(:));
            Info.NepochRejectedPipeline = Before - sum(FlagEpoch);
        end
    end

    [Obj, Secz, PA, HA, Alt] = selectEpochs(Obj, FlagEpoch, Secz, PA, HA, Alt);
    EpochMap = EpochMap(FlagEpoch);

    % --- source selection ---------------------------------------------------
    Ndet     = sum(~isnan(Obj.Data.X), 1);
    MinNdet  = max(Args.MinNdetSrc, ceil(Args.MinDetFrac.*Obj.Nepoch));
    FlagSrc  = Ndet >= MinNdet;
    Info.MinNdetUsed      = MinNdet;
    Info.NsrcRejectedNdet = sum(~FlagSrc);

    % Scatter cut: a source whose position wanders by more than a pixel is
    % either blended, mismatched, or genuinely variable in position, and it
    % degrades the per-epoch affine solution that every other source shares.
    Info.NsrcRejectedStd = 0;
    if isfinite(Args.MaxSrcStdXY)
        StdX    = std(Obj.Data.X, [], 1, 'omitnan');
        StdY    = std(Obj.Data.Y, [], 1, 'omitnan');
        Before  = sum(FlagSrc);
        FlagSrc = FlagSrc & StdX < Args.MaxSrcStdXY & StdY < Args.MaxSrcStdXY;
        Info.NsrcRejectedStd = Before - sum(FlagSrc);
    end
    Info.NsrcRejectedPipeline = 0;
    if Args.ApplyPipelineFlags
        FlagGoodSrcUD = getUserField(UserData,'FlagGoodSrc');
        if numel(FlagGoodSrcUD) == Nsrc0
            Before  = sum(FlagSrc);
            FlagSrc = FlagSrc & logical(FlagGoodSrcUD(:)).';
            Info.NsrcRejectedPipeline = Before - sum(FlagSrc);
        end
    end
    % The unified source list can hold exact duplicates of a source: two
    % columns with bit-identical positions and fluxes. They would give two
    % identical median magnitudes, and ml.util.IterativeMovingMedianEpoch
    % passes those medians to isoutlier as 'SamplePoints', which rejects
    % duplicates - leaving @IterFit/calculateWes to fall back to unit weights.
    Info.NsrcRejectedDuplicate = 0;
    if Args.RemoveDuplicateSrc
        FlagUnique = flagUniqueSources(Obj);
        Info.NsrcRejectedDuplicate = sum(FlagSrc & ~FlagUnique);
        FlagSrc    = FlagSrc & FlagUnique;
    end

    Obj.selectBySrcIndex(FlagSrc, 'CreateNewObj', false);
    SrcData  = selectSrcData(SrcData, FlagSrc);
    SrcInd   = find(FlagSrc);

    % A second epoch pass: dropping sources can leave epochs too empty to
    % constrain their own affine transformation.
    Nsrcep   = sum(~isnan(Obj.Data.X), 2);
    FlagEp2  = Nsrcep >= Args.MinNsrcEpoch;
    Info.NepochRejectedNsrc = sum(~FlagEp2);
    [Obj, Secz, PA, HA, Alt] = selectEpochs(Obj, FlagEp2, Secz, PA, HA, Alt);
    EpochMap = EpochMap(FlagEp2);

    if Obj.Nepoch < 2 || Obj.Nsrc < 2
        error('mmsFromMatchedSources:EmptySelection','No data left after the epoch/source cuts (Nepoch=%d, Nsrc=%d)', Obj.Nepoch, Obj.Nsrc);
    end

    Info.EpochInd = EpochMap;
    Info.SrcInd   = SrcInd;
    Info.Nepoch   = Obj.Nepoch;
    Info.Nsrc     = Obj.Nsrc;
    Info.FracNaN = mean(isnan(Obj.Data.X(:)));
    report(Args.Verbosity, 'Kept %d/%d epochs and %d/%d sources (Ndet>=%d, fill %.1f%%)\n', ...
           Obj.Nepoch, Nepoch0, Obj.Nsrc, Nsrc0, MinNdet, 100.*(1-Info.FracNaN));

    % --- expand the per-epoch geometry to [Nepoch, Nsrc] --------------------
    % @IterFit/generateHALatDesignMat reads Data.pa(:,1) and Data.secz(:,1),
    % while MatchedSources derives Nsrc/Nepoch from the size of the first Data
    % field and ml.util.flag_struct_field indexes every field by column, so
    % these must be full matrices rather than column vectors.
    Ones            = ones(1, Obj.Nsrc);
    Obj.Data.secz   = Secz * Ones;
    Obj.Data.pa     = PA   * Ones;
    Obj.Data.ha     = HA   * Ones;
    Obj.Data.alt    = Alt  * Ones;

    % --- pixel phase, only when the pipeline stored the per-epoch shifts ----
    % X and Y are the registered (shifted) coordinates, so the detector pixel
    % phase is recoverable only if the shift that was applied is known.
    [Obj, Info.PixPhaseAvailable] = addPixelPhase(Obj, UserData, EpochMap, Nepoch0);

    % --- magnitudes ---------------------------------------------------------
    [Obj, Info.Phot] = buildMagnitudes(Obj, SrcData, Args);

    % Sources whose median magnitude is not finite, or which still share it
    % with another source, would break the same outlier rejection.
    MedMag = Obj.medianFieldSource({Args.ColNameMag}).';
    if Args.RequireFiniteMag
        FlagMag = isfinite(MedMag);
        Info.NsrcRejectedMag = sum(~FlagMag);
        if any(~FlagMag)
            Obj.selectBySrcIndex(FlagMag, 'CreateNewObj', false);
            SrcData     = selectSrcData(SrcData, FlagMag);
            Info.SrcInd = Info.SrcInd(FlagMag);
            Info.Nsrc   = Obj.Nsrc;
            MedMag      = MedMag(FlagMag);
        end
    else
        Info.NsrcRejectedMag = 0;
    end
    [Obj, Info.MagTiesBroken] = breakMagTies(Obj, MedMag, Args.ColNameMag);

    % --- colour -------------------------------------------------------------
    [Obj, SrcData, Info.Colour, FlagColourSrc] = buildColour(Obj, SrcData, Args);
    if ~isempty(FlagColourSrc)
        Info.SrcInd = Info.SrcInd(FlagColourSrc);
        Info.Nsrc   = Obj.Nsrc;
    end

    Info.SrcData = SrcData;
end


% -------------------------------------------------------------------------
function MS = readMatchedSourcesFile(FileName)
    % Load the first MatchedSources variable found in a .mat file
    Loaded = load(FileName);
    Vars   = fieldnames(Loaded);
    IsMS   = cellfun(@(V) isa(Loaded.(V),'MatchedSources'), Vars);
    if ~any(IsMS)
        error('mmsFromMatchedSources:NoObject','No MatchedSources object found in %s', FileName);
    end
    MS = Loaded.(Vars{find(IsMS,1)});
end


function [Obj, Secz, PA, HA, Alt] = selectEpochs(Obj, Flag, Secz, PA, HA, Alt)
    % Apply an epoch flag to the object and to the per-epoch geometry vectors
    if all(Flag)
        return
    end
    Obj.selectByEpoch(Flag, 'CreateNewObj', false);
    Secz = Secz(Flag);
    PA   = PA(Flag);
    HA   = HA(Flag);
    Alt  = Alt(Flag);
end


function SrcData = selectSrcData(SrcData, Flag)
    % Apply a source flag to every field of a SrcData structure
    if isempty(SrcData) || ~isstruct(SrcData)
        SrcData = struct();
        return
    end
    Fields = fieldnames(SrcData);
    for Ifield = 1:numel(Fields)
        SrcData.(Fields{Ifield}) = SrcData.(Fields{Ifield})(Flag);
    end
end


function [Obj, IsAvailable] = addPixelPhase(Obj, UserData, EpochInd, Nepoch0)
    % Restore the detector pixel phase from the registration shifts.
    % X and Y are the registered coordinates, so the phase is recoverable only
    % when the shift that was applied to each epoch is known.
    [ShiftX, HasX] = alignEpochVector(getUserField(UserData,'ShiftX'), UserData, Nepoch0);
    [ShiftY, HasY] = alignEpochVector(getUserField(UserData,'ShiftY'), UserData, Nepoch0);
    IsAvailable    = HasX && HasY;
    if ~IsAvailable
        return
    end
    Obj.Data.Xphase = mod(Obj.Data.X + ShiftX(EpochInd), 1) - 0.5;
    Obj.Data.Yphase = mod(Obj.Data.Y + ShiftY(EpochInd), 1) - 0.5;
end


function FlagUnique = flagUniqueSources(Obj)
    % Flag the first occurrence of every distinct median position
    MedXY = [median(Obj.Data.X, 1, 'omitnan').', median(Obj.Data.Y, 1, 'omitnan').'];
    [~, IndFirst] = unique(MedXY, 'rows', 'stable');
    FlagUnique = false(1, size(MedXY,1));
    FlagUnique(IndFirst) = true;
end


function [Obj, WasBroken] = breakMagTies(Obj, MedMag, ColNameMag)
    % Separate sources that share a median magnitude, by an amount far below
    % any photometric significance, so that it can be used as sample points.
    WasBroken = false;
    [~, ~, IndUnique] = unique(MedMag);
    Count = accumarray(IndUnique, 1);
    IsTied = Count(IndUnique) > 1;
    if ~any(IsTied)
        return
    end
    Offset         = zeros(1, numel(MedMag));
    Offset(IsTied) = (1:sum(IsTied)).*1e-9;
    Obj.Data.(ColNameMag) = Obj.Data.(ColNameMag) + Offset;
    WasBroken = true;
end


function Value = getUserField(UserData, FieldName)
    % Field of a UserData structure, or empty when it is not there
    if isstruct(UserData) && isfield(UserData, FieldName)
        Value = UserData.(FieldName);
    else
        Value = [];
    end
end


function [Vec, IsAligned] = alignEpochVector(Vec, UserData, Nepoch0)
    % Bring a per-epoch vector stored in UserData onto the epochs of the
    % object, which may or may not already have had FlagGoodEpoch applied.
    Vec       = Vec(:);
    IsAligned = numel(Vec) == Nepoch0;
    if IsAligned || isempty(Vec)
        IsAligned = IsAligned && ~isempty(Vec);
        return
    end
    FlagGoodEpoch = getUserField(UserData, 'FlagGoodEpoch');
    if numel(FlagGoodEpoch) == numel(Vec) && sum(logical(FlagGoodEpoch)) == Nepoch0
        Vec       = Vec(logical(FlagGoodEpoch));
        IsAligned = true;
    end
end


function [Obj, Phot] = buildMagnitudes(Obj, SrcData, Args)
    % Build a magnitude field from the flux, anchor its zero point on the
    % reference magnitude and remove the remaining photometric systematics.
    Phot = struct('InstrumentalOnly', true, 'RefZPApplied', false, 'SysRemApplied', false, 'Message', '');

    Flux = Obj.Data.(Args.ColNameFlux);
    Flux(Flux <= 0) = NaN;
    Obj.Data.(Args.ColNameMag) = Args.ZP0 - 2.5.*log10(Flux);

    HasRef = isstruct(SrcData) && isfield(SrcData, Args.RefMagField);
    if Args.ApplyRefZP && HasRef
        RefMag = SrcData.(Args.RefMagField)(:).';
        RefMag(RefMag >= Args.MaxValidRefMag) = NaN;
        Obj.Data.RefMag = ones(Obj.Nepoch,1) * RefMag;
        try
            ZP = Obj.fitRefZP('ColNameMag', Args.ColNameMag, 'ColNameRefMag', 'RefMag', Args.fitRefZPArgs{:});
            Obj.applyZP(ZP, 'ApplyToMagField', Args.ColNameMag);
            Obj.ZP = ZP;
            Obj.Data = rmfield(Obj.Data, 'RefMag');
            Phot.RefZPApplied     = true;
            Phot.InstrumentalOnly = false;
            report(Args.Verbosity, 'Applied per-epoch zero point anchored on %s\n', Args.RefMagField);
        catch ME
            Phot.Message = sprintf('fitRefZP failed (%s); magnitudes left instrumental', ME.message);
            fprintf('%s\n', Phot.Message);
            if isfield(Obj.Data, 'RefMag')
                Obj.Data = rmfield(Obj.Data, 'RefMag');
            end
        end
    elseif Args.ApplyRefZP
        Phot.Message = sprintf('SrcData has no field "%s"; magnitudes left instrumental', Args.RefMagField);
        fprintf('%s\n', Phot.Message);
    end

    if Args.ApplySysRemPhot
        try
            CorrectedMag = SysRemPhotometry(Obj, 'ColNameMag', Args.ColNameMag, Args.SysRemPhotometryArgs{:});
            Obj.Data.(Args.ColNameMag) = CorrectedMag;
            Phot.SysRemApplied = true;
            report(Args.Verbosity, 'Applied SysRem to the magnitudes\n');
        catch ME
            Phot.Message = sprintf('SysRemPhotometry failed (%s); magnitudes left uncorrected', ME.message);
            fprintf('%s\n', Phot.Message);
        end
    end
end


function [Obj, SrcData, Colour, FlagColourSrc] = buildColour(Obj, SrcData, Args)
    % Build the per-source colour used by the IterFit colour bins
    Colour        = struct();
    FlagColourSrc = [];

    Blue = Args.ColourFields{1};
    Red  = Args.ColourFields{2};
    if ~isstruct(SrcData) || ~isfield(SrcData, Blue) || ~isfield(SrcData, Red)
        error('mmsFromMatchedSources:NoColour','SrcData must hold the colour fields "%s" and "%s"', Blue, Red);
    end

    BlueMag = SrcData.(Blue)(:).';
    RedMag  = SrcData.(Red)(:).';
    BlueMag(BlueMag >= Args.MaxValidRefMag) = NaN;
    RedMag(RedMag   >= Args.MaxValidRefMag) = NaN;
    C       = BlueMag - RedMag;
    IsValid = isfinite(C) & C >= Args.ColourRange(1) & C <= Args.ColourRange(2);

    Colour.NValid   = sum(IsValid);
    Colour.FracValid = mean(IsValid);
    Colour.Median   = median(C(IsValid));
    Colour.Mode     = Args.ColourMode;
    if Colour.NValid < 2
        error('mmsFromMatchedSources:NoValidColour','Fewer than two sources have a valid colour');
    end
    report(Args.Verbosity, 'Valid colour for %d/%d sources (%.1f%%)\n', Colour.NValid, numel(C), 100.*Colour.FracValid);

    switch lower(Args.ColourMode)
        case 'restrict'
            FlagColourSrc = IsValid;
            Obj.selectBySrcIndex(FlagColourSrc, 'CreateNewObj', false);
            SrcData = selectSrcData(SrcData, FlagColourSrc);
            C       = C(FlagColourSrc);

        case 'fill'
            C(~IsValid) = Colour.Median;

        case 'fillcmd'
            C = fillColourFromMag(C, IsValid, Obj, Args);

        case 'ownbin'
            % One bin width below every valid colour, so that discretize
            % gathers the colourless sources into the lowest bin.
            ValidSpread = spread(C(IsValid));
            C(~IsValid) = min(C(IsValid)) - max(ValidSpread./Args.NColourBins, 0.1);

        otherwise
            error('mmsFromMatchedSources:BadColourMode','Unknown ColourMode "%s"', Args.ColourMode);
    end

    % @IterFit/generateBins takes quantiles of the colour and feeds them to
    % discretize, which rejects repeated edges. Any mode that assigns one
    % value to many sources therefore needs those ties separated first.
    [C, Colour.TiesBroken] = breakColourTies(C, Args.NColourBins);

    if Args.CentreColour
        C = C - median(C, 'omitnan');
    end
    Colour.C = C;

    Obj.Data.C = ones(Obj.Nepoch,1) * C;
end


function C = fillColourFromMag(C, IsValid, Obj, Args)
    % Predict the missing colours from a robust colour-magnitude relation
    Mag = Obj.medianFieldSource({Args.ColNameMag}).';
    Fit = isfinite(Mag) & IsValid;
    if sum(Fit) < 10
        C(~IsValid) = median(C(IsValid));
        return
    end
    Par = robustfit(Mag(Fit), C(Fit));
    Missing = ~IsValid & isfinite(Mag);
    C(Missing) = Par(1) + Par(2).*Mag(Missing);
    C(~IsValid & ~isfinite(Mag)) = median(C(IsValid));
    C = min(max(C, Args.ColourRange(1)), Args.ColourRange(2));
end


function [C, WasBroken] = breakColourTies(C, NBins)
    % Add a monotone, deterministic ramp so that the quantile bin edges used
    % by @IterFit/generateBins are strictly increasing.
    WasBroken = false;
    Edges     = quantile(C, linspace(0, 1, NBins+1));
    if all(diff(Edges) > 0)
        return
    end
    Span = spread(C);
    if ~isfinite(Span) || Span == 0
        Span = 1;
    end
    Step        = 1e-6.*Span;
    [~, Isort]  = sort(C, 'ascend');
    Ramp        = zeros(size(C));
    Ramp(Isort) = (0:numel(C)-1).*Step;
    C           = C + Ramp;
    WasBroken   = true;
end


function Result = spread(Vec)
    % Range of the finite elements of a vector
    Vec    = Vec(isfinite(Vec));
    if isempty(Vec)
        Result = 0;
    else
        Result = max(Vec) - min(Vec);
    end
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
