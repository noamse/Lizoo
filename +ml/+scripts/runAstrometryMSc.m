function [OutputFileName, File] = runAstrometryMSc(MS, Args)
% Astrometric detrending of a KMT_pipelineI MatchedSources, with optional Gaia calibration.
% This is the ml.scripts.runAstrometryField analogue for input that has
% already been matched by KMT_pipelineI. Step 1 is the detrending chain of
% ml.scripts.runIterDetrendMSc; step 2, the Gaia proper-motion calibration, is
% off by default because the saved MatchedSources carries no sky coordinates,
% and when switched on it first solves a WCS from the source list itself.
% Input  : - A MatchedSources object as produced by KMT_pipelineI, or a
%            char/string with the path of a .mat file containing one.
%          * ...,key,val,...
%            'EventNum' - Event number recorded in the output. Default is [].
%            'Site' - Site name recorded in the output. Default is 'CTIO'.
%            'Field' - Field name recorded in the output. Default is ''.
%            'runIterDetrendMScArgs' - Cell array of arguments passed on to
%                   ml.scripts.runIterDetrendMSc. Default is {}.
%            'PerSourcesTargetPath' - Directory for the per-source csv tables.
%                   No tables are written when empty. Default is ''.
%            'GaiaCalib' - Run the Gaia proper-motion calibration. Requires
%                   catsHTM and a successful WCS solution. Default is false.
%            'Scale' - Pixel scale [arcsec/pix] handed to astrometryCore.
%                   Default is 0.4.
%            'astrometryCoreArgs' - Cell array of arguments passed on to
%                   imProc.astrometry.astrometryCore. Default is {}.
%            'GaiaCatMatchedFile' - Passed on to ml.scripts.gaiaAstrometryKMT.
%                   Default is [].
%            'RstdTopPrc' - Passed on to ml.scripts.gaiaAstrometryKMT.
%                   Default is 33.
%            'MaxMag' - Passed on to ml.scripts.gaiaAstrometryKMT. Default is 18.
%            'RefMagField' - SrcData field holding the reference magnitude
%                   inserted as column 'I' of the matched catalog.
%                   Default is 'I_ogle'.
%            'Save' - Save the output structure. Default is true.
%            'OutputDir' - Output directory. When empty and the input was a
%                   file, an AstrometryMSc directory next to it is used.
%                   Default is ''.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - Full path of the saved file, empty when nothing was saved.
%          - The output structure, with the same field names as the one
%            written by ml.scripts.runAstrometryField.
% Author : ULTRASAT team (Aug 2026)
% Example: OutFile = ml.scripts.runAstrometryMSc(...
%              '/bigdata3/projects/KMTdata/Results16_26/KMT_260058_BLG41_MSc.mat', ...
%              'EventNum',260058, 'Field','BLG41', 'Verbosity',1);

    arguments
        MS
        Args.EventNum                  = [];
        Args.Site                      = 'CTIO';
        Args.Field                     = '';
        Args.runIterDetrendMScArgs     = {};
        Args.PerSourcesTargetPath      = '';
        Args.GaiaCalib                 = false;
        Args.Scale                     = 0.4;
        Args.astrometryCoreArgs        = {};
        Args.GaiaCatMatchedFile        = [];
        Args.RstdTopPrc                = 33;
        Args.MaxMag                    = 18;
        Args.RefMagField               = 'I_ogle';
        Args.Save                      = true;
        Args.OutputDir                 = '';
        Args.Verbosity                 = 0;
    end

    OutputFileName = '';
    if ischar(MS) || isstring(MS)
        SourceFile = char(MS);
    else
        SourceFile = '';
    end

    % --- Step 1: convert and detrend ---------------------------------------
    [IFsys, Obj, IFsysB, Info] = ml.scripts.runIterDetrendMSc(MS, ...
        'PerSourcesTargetPath', Args.PerSourcesTargetPath, ...
        'Verbosity', Args.Verbosity, ...
        Args.runIterDetrendMScArgs{:});

    % --- Step 2: Gaia proper-motion calibration (opt-in) --------------------
    Matched         = AstroCatalog;
    ResAstrometry   = [];
    ParScalibrated  = [];
    T               = [];
    DeltaPM_KMT_GAIA = [];
    OutLiersRMSvsMag = [];
    PMRA_kmt_to_gaia_fit  = [];
    PMDec_kmt_to_gaia_fit = [];

    if Args.GaiaCalib
        try
            [Matched, ResAstrometry] = buildMatchedCatalog(IFsys, Info, Args);
            [ParScalibrated, T, DeltaPM_KMT_GAIA, OutLiersRMSvsMag, ...
             PMRA_kmt_to_gaia_fit, PMDec_kmt_to_gaia_fit] = ...
                ml.scripts.gaiaAstrometryKMT(IFsys, Matched, ...
                    'GaiaCatMatchedFile', Args.GaiaCatMatchedFile, ...
                    'RstdTopPrc', Args.RstdTopPrc, ...
                    'MaxMag', Args.MaxMag);
            report(Args.Verbosity, 'Gaia calibration completed\n');
        catch ME
            fprintf('Gaia calibration failed: %s\n', ME.message);
        end
    end

    % --- Step 3: reference source and field centre --------------------------
    IndForPhotRefernce = IFsys.findClosestSource([150, 150]);
    Coo = IFsys.CelestialCoo .* 180 ./ pi;

    % --- Step 4: package ----------------------------------------------------
    File                       = struct();
    File.EventNum              = Args.EventNum;
    File.Site                  = Args.Site;
    File.Field                 = Args.Field;
    File.SourceFile            = SourceFile;
    File.CatsPath              = Info.CatsPath;
    File.IFsys                 = IFsys;
    File.IFsysBeforeSysRem     = IFsysB;
    File.Obj                   = Obj;
    File.Info                  = Info;
    File.Matched               = Matched;
    File.ResAstrometry         = ResAstrometry;
    File.ParScalibrated        = ParScalibrated;
    File.GaiaTable             = T;
    File.DeltaPM_KMT_GAIA      = DeltaPM_KMT_GAIA;
    File.OutLiersRMSvsMag      = OutLiersRMSvsMag;
    File.PMRA_kmt_to_gaia_fit  = PMRA_kmt_to_gaia_fit;
    File.PMDec_kmt_to_gaia_fit = PMDec_kmt_to_gaia_fit;
    File.IndForPhotRefernce    = IndForPhotRefernce;
    File.FieldCenterDeg        = Coo;

    % --- Step 5: save -------------------------------------------------------
    if Args.Save
        OutputDir = Args.OutputDir;
        if isempty(OutputDir)
            if isempty(SourceFile)
                error('runAstrometryMSc:NoOutputDir','OutputDir must be given when the input is an object rather than a file');
            end
            OutputDir = fullfile(fileparts(SourceFile), 'AstrometryMSc');
        end
        if ~isfolder(OutputDir)
            mkdir(OutputDir);
        end
        OutputFileName = fullfile(OutputDir, sprintf('AstrometryMSc_%s_%s_%s.mat', ...
                                                     num2str(Args.EventNum), Args.Site, Args.Field));
        if isfile(OutputFileName)
            delete(OutputFileName);
        end
        save(OutputFileName, 'File', '-v7.3');
        fprintf('Saved astrometry file to %s\n', OutputFileName);
    end
end


% -------------------------------------------------------------------------
function [Matched, ResAstrometry] = buildMatchedCatalog(IFsys, Info, Args)
    % Solve a WCS from the fitted source positions and return a catalog with
    % RA, Dec and the reference magnitude, one row per IterFit source.
    X   = IFsys.ParS(1,:).';
    Y   = IFsys.ParS(2,:).';
    Mag = IFsys.medianFieldSource({'MAG_PSF'});

    Cat = AstroCatalog({[X, Y, Mag]}, 'ColNames', {'X','Y','MAG_PSF'}, ...
                                      'ColUnits', {'pix','pix','mag'});

    RAdeg  = IFsys.CelestialCoo(1) .* 180 ./ pi;
    Decdeg = IFsys.CelestialCoo(2) .* 180 ./ pi;

    [ResAstrometry, Matched] = imProc.astrometry.astrometryCore(Cat, ...
        'RA', RAdeg, 'Dec', Decdeg, 'CooUnits', 'deg', ...
        'Scale', Args.Scale, 'CatColNamesMag', 'MAG_PSF', ...
        Args.astrometryCoreArgs{:});

    if ~ResAstrometry.Success
        error('runAstrometryMSc:NoWCS','astrometryCore did not converge on a WCS (%d candidate solutions)', ResAstrometry.Nsolutions);
    end

    % gaiaAstrometryKMT reads the reference magnitude as column 'I'
    if isstruct(Info.SrcData) && isfield(Info.SrcData, Args.RefMagField)
        RefMag = Info.SrcData.(Args.RefMagField)(:);
    else
        RefMag = nan(size(X));
    end
    Matched.insertCol(RefMag, Inf, {'I'}, {'mag'});
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
