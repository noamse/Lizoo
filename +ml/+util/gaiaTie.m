function [Tie, Info] = gaiaTie(IF, Args)
% Tie a KMT solution to Gaia, through OGLE, and calibrate its proper motions.
% The astrometry of this package is relative: a free per-epoch transformation
% absorbs any motion shared by the field, so proper motions come out measured
% against the mean motion of whichever sources happen to be in the cut-out, and
% the plate scale and the orientation on the sky are assumed rather than known.
% Gaia fixes all three. What it does not fix is the per-source precision, which
% is set by the per-epoch noise: measured on these fields the scatter about the
% Gaia relation is 3.4 to 3.8 mas/yr against a Gaia spread of 2.9 to 3.3, so
% most of the difference is ours.
% The match runs through OGLE rather than straight to Gaia. A KMT detection is
% usually a blend, a median of two OGLE stars inside half a seeing disc, so its
% centroid sits between real stars and no pattern matcher will place it: tried
% directly, 64 of 594 sources match and the two fields disagree about the
% rotation by 47 degrees. OGLE resolves the field, Gaia resolves the field, and
% KMT is already tied to OGLE by the pipeline's own matching, so chaining the
% three works: 532 of 594 sources match and the fields agree to 0.2 degrees.
% Input  : - A fitted IterFit object, or the path of a .mat holding one in
%            IFsys.
%          * ...,key,val,...
%            'Field' - 'BLG41' or 'BLG01'. Selects the OGLE offset. Required.
%            'RA','Dec' - Centre of the cone search [deg]. Default is the
%                   OGLE-2026-BLG-0058 position. These must be the position of
%                   the field, not the telescope pointing: they differ by 44
%                   arcmin here, and the wrong one returns a patch of empty sky
%                   that matches nothing.
%            'SearchRadius' - Cone radius [arcsec]. Default is 120.
%            'OgleFile' - The OGLE reference catalogue.
%            'OgleScale' - OGLE pixel scale [arcsec]. Default is 0.26.
%            'KmtScale' - KMT pixel scale [arcsec]. Default is 0.4.
%            'MatchRadiusOgle' - OGLE to Gaia tolerance [arcsec]. Default 0.5.
%            'MatchRadiusKmt' - KMT to Gaia tolerance [arcsec]. Default 0.6.
%            'MaxPMErr' - Gaia sources with a larger proper-motion error take
%                   no part in the calibration [mas/yr]. Default is 0.4.
%            'MaxMag' - KMT sources fainter than this take no part. Default 18.
%            'OgleSolution' - A previous solve's Tie.OgleSolution, reused
%                   instead of repeating the pattern search. That search is the
%                   expensive part and does not depend on the KMT field, so on
%                   the second field it is pure waste. Default is [], solve it.
%            'SearchMaxG' - Only Gaia sources brighter than this drive the
%                   pattern search, which is quadratic in the catalogue sizes.
%                   The refinement afterwards uses them all. Default is 19.5.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - A structure holding the calibration: the 2 by 2 matrix and the
%            offset carrying this solution's proper motions onto Gaia's, the
%            implied rotation and scale, and the absolute proper motion of every
%            matched source.
%          - A structure recording how many sources matched at each step and the
%            residual left at each, so that a bad tie is visible rather than
%            silent.
% Author : ULTRASAT team (Sep 2026)
% <Unknown Bugs>: the OGLE to Gaia step assumes the OGLE pixel grid is a pure
%   rotation and scaling of the sky, which over 2 arcmin it is to well within
%   the 105 to 123 mas that step leaves behind. That residual, not Gaia, sets
%   the floor of the tie.
% Example: Tie = ml.util.gaiaTie(IFsys,'Field','BLG41','Verbosity',1);

    arguments
        IF
        Args.Field
        Args.RA              = celestial.coo.convertdms('17:52:38.09','gH','d');
        Args.Dec             = celestial.coo.convertdms('-31:47:36.1','gD','d');
        Args.SearchRadius    = 120;
        Args.OgleFile        = '/home/ocs/matlab/Lizoo/OGLEdata/OB260058/OB160058.mat';
        Args.OgleScale       = 0.26;
        Args.KmtScale        = 0.4;
        Args.MatchRadiusOgle = 0.5;
        Args.MatchRadiusKmt  = 0.6;
        Args.MaxPMErr        = 0.4;
        Args.MaxMag          = 18;
        Args.OgleSolution    = [];   % reuse a previous OGLE-to-Gaia solve
        Args.SearchMaxG      = 19.5; % only these Gaia sources drive the pattern search
        Args.Verbosity       = 0;
    end

    if ischar(IF) || isstring(IF)
        L = load(char(IF), 'IFsys');
        IF = L.IFsys;
    end
    switch upper(Args.Field)
        case 'BLG41', OgleOff = [229.202, 229.283];
        case 'BLG01', OgleOff = [227.542, 229.689];
        otherwise
            error('ml:util:gaiaTie:UnknownField', ...
                  'No OGLE offset is recorded for field %s', Args.Field);
    end

    % --- Gaia ------------------------------------------------------------
    [GC, Col] = catsHTM.cone_search('GAIADR3', Args.RA./180.*pi, Args.Dec./180.*pi, Args.SearchRadius);
    col = @(n) GC(:, strcmp(Col, n));
    Xi  = (col('RA').*180./pi  - Args.RA ).*cosd(Args.Dec).*3600;
    Eta = (col('Dec').*180./pi - Args.Dec).*3600;
    PMRA = col('PMRA');  PMDec = col('PMDec');
    EPMRA= col('ErrPMRA'); EPMDec = col('ErrPMDec');
    Info.NGaia = size(GC,1);

    % --- OGLE, and its tie to Gaia ---------------------------------------
    Loaded = load(Args.OgleFile);
    Vars   = fieldnames(Loaded);
    Tab    = Loaded.(Vars{1}).CatData.Table;
    Oxp = double(Tab.corrX);  Oyp = double(Tab.corrY);  OI = double(Tab.I);
    Keep = isfinite(Oxp) & isfinite(Oyp) & isfinite(OI) & OI < 19;
    OxCen = mean(Oxp(Keep));  OyCen = mean(Oyp(Keep));
    Ox = (Oxp(Keep)-OxCen).*Args.OgleScale;
    Oy = (Oyp(Keep)-OyCen).*Args.OgleScale;
    if isempty(Args.OgleSolution)
        [ROgle, Info.Ogle] = solveOgleGaia(Ox, Oy, Xi, Eta, col('phot_g_mean_mag'), Args);
    else
        ROgle = Args.OgleSolution;
        Info.Ogle = struct('Reused',true);
        report(Args.Verbosity, 'reusing a previous OGLE-to-Gaia solution\n');
    end
    Tie.OgleSolution = ROgle;
    if isfield(Info.Ogle,'Nmatch')
        report(Args.Verbosity, 'OGLE to Gaia: %d of %d matched, residual %.0f/%.0f mas\n', ...
               Info.Ogle.Nmatch, numel(Ox), 1000.*Info.Ogle.RmsX, 1000.*Info.Ogle.RmsY);
    end

    % --- KMT, carried through OGLE onto the Gaia frame --------------------
    Xk = IF.ParS(1,:).';  Yk = IF.ParS(2,:).';
    Mag = IF.medianFieldSource({'MAG_PSF'});  Mag = Mag(:);
    OxK = (Xk-150).*Args.KmtScale./Args.OgleScale + OgleOff(1);
    OyK = (Yk-150).*Args.KmtScale./Args.OgleScale + OgleOff(2);
    AxK = (OxK-OxCen).*Args.OgleScale;
    AyK = (OyK-OyCen).*Args.OgleScale;
    W   = ROgle.M \ [AxK.'-ROgle.c(1); AyK.'-ROgle.c(2)];
    KXi = W(1,:).';  KEta = W(2,:).';

    [D, J] = nearest(KXi, KEta, Xi, Eta);
    Matched = D < Args.MatchRadiusKmt;
    Info.NKmt      = numel(Xk);
    Info.NMatched  = sum(Matched);
    Info.MedianSep = median(D(Matched));
    report(Args.Verbosity, 'KMT to Gaia: %d of %d matched within %.1f arcsec, median %.3f\n', ...
           Info.NMatched, numel(Xk), Args.MatchRadiusKmt, Info.MedianSep);

    % --- calibrate the proper motions -------------------------------------
    PMx = 400.*IF.ParS(3,:).';  PMy = 400.*IF.ParS(4,:).';
    Use = Matched & isfinite(PMRA(J)) & isfinite(PMDec(J)) & ...
          0.5.*(EPMRA(J)+EPMDec(J)) < Args.MaxPMErr & Mag < Args.MaxMag;
    Info.NCalib = sum(Use);
    if Info.NCalib < 30
        error('ml:util:gaiaTie:TooFewCalibrators', ...
              'Only %d sources carry both a KMT and a good Gaia proper motion', Info.NCalib);
    end
    G = [ones(sum(Use),1), PMRA(J(Use)), PMDec(J(Use))];
    ax = G\PMx(Use);  ay = G\PMy(Use);
    Rx = PMx(Use)-G*ax;  Ry = PMy(Use)-G*ay;
    Tie.A      = [ax(2) ax(3); ay(2) ay(3)];
    Tie.Offset = [ax(1); ay(1)];
    Tie.Rotation = atan2d(ay(2)-ax(3), ax(2)+ay(3));
    Tie.Scale    = sqrt(abs(det(Tie.A)));
    Tie.ResidX   = tools.math.stat.rstd(Rx);
    Tie.ResidY   = tools.math.stat.rstd(Ry);
    Tie.Matched  = Matched;
    Tie.GaiaInd  = J;
    % every matched source, carried onto the Gaia system
    Abs = nan(numel(Xk),2);
    Ai  = Tie.A \ ([PMx(Matched) PMy(Matched)].' - Tie.Offset);
    Abs(Matched,:) = Ai.';
    Tie.PMabs = Abs;
    Tie.GaiaPM = [PMRA(J), PMDec(J)];
    report(Args.Verbosity, ['proper-motion tie on %d sources: rotation %+.2f deg, scale %.3f, ', ...
           'gauge %+.2f/%+.2f mas/yr, residual %.2f/%.2f\n'], Info.NCalib, Tie.Rotation, Tie.Scale, ...
           Tie.Offset(1), Tie.Offset(2), Tie.ResidX, Tie.ResidY);
end


function [R, Rep] = solveOgleGaia(Ox, Oy, Xi, Eta, Gmag, Args)
    % Rotation, parity and shift carrying Gaia standard coordinates onto the
    % OGLE grid, found from the histogram of pairwise offsets, then refined.
    Bright = isfinite(Gmag) & Gmag < Args.SearchMaxG;   % the search is quadratic
    Xs = Xi(Bright);  Es = Eta(Bright);
    Best = struct('n',-1);
    for Flip = [1 -1]
        for Th = -180:0.5:179.5
            Rot = [cosd(Th) -sind(Th); sind(Th) cosd(Th)];
            P   = Rot*[Flip.*Xs.'; Es.'];
            Du  = Ox-P(1,:);  Dv = Oy-P(2,:);
            Du  = Du(:);  Dv = Dv(:);
            In  = abs(Du)<60 & abs(Dv)<60;
            if sum(In)<50, continue; end
            H = histcounts2(Du(In), Dv(In), -60:0.4:60, -60:0.4:60);
            [N, Ix] = max(H(:));
            if N > Best.n
                [a,b] = ind2sub(size(H), Ix);
                Best = struct('n',N,'Th',Th,'Flip',Flip, ...
                              'u0',-60+0.4.*(a-0.5),'v0',-60+0.4.*(b-0.5),'Npair',sum(In));
            end
        end
    end
    Rot = [cosd(Best.Th) -sind(Best.Th); sind(Best.Th) cosd(Best.Th)];
    P   = Rot*[Best.Flip.*Xi.'; Eta.'];
    [D, J] = nearest(Ox, Oy, P(1,:).'+Best.u0, P(2,:).'+Best.v0);
    Ok = D < Args.MatchRadiusOgle;
    A  = [ones(sum(Ok),1), Xi(J(Ok)), Eta(J(Ok))];
    px = A\Ox(Ok);  py = A\Oy(Ok);
    R.M = [px(2) px(3); py(2) py(3)];
    R.c = [px(1); py(1)];
    Rep.Nmatch = sum(Ok);
    Rep.RmsX   = rms(Ox(Ok)-A*px);
    Rep.RmsY   = rms(Oy(Ok)-A*py);
    Rep.Peak   = Best.n;
    Rep.Background = Best.Npair./(300.*300);
end


function [D, J] = nearest(X, Y, Xr, Yr)
    % Nearest (Xr,Yr) to every (X,Y)
    D = nan(numel(X),1);  J = nan(numel(X),1);
    for K = 1:numel(X)
        [D(K), J(K)] = min(hypot(Xr-X(K), Yr-Y(K)));
    end
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
