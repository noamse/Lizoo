function [ObjOut, CorrX, CorrY, Info] = correctSourceDCR(IF, Obj, Args)
% Remove a per-source differential chromatic refraction term from an IterFit solution.
% @IterFit fits the airmass terms in colour bins, six of them by default, so a
% source whose colour differs from its bin keeps a residual that follows
% sin(pa)*secz and cos(pa)*secz. Because the parallactic angle sweeps through a
% season, that residual appears as a seasonal wander: on KMT-2026-BLG-0058 the
% June 2019 monthly median reached -50 mas with no photometric change, and a
% per-source chromatic fit removed all but a few mas of it.
% This fits the two chromatic terms separately for every source and subtracts
% them. Only those two are subtracted; the constant is fitted so that it
% absorbs any offset rather than leaking into the chromatic coefficients, but
% it is not applied, which leaves each source's mean position untouched.
% Input  : - An IterFit object holding a solution.
%          - An MMS object whose X and Y are to be corrected. Not modified.
%          * ...,key,val,...
%            'FitJDRange' - [JDmin JDmax] of the epochs the fit may use. All
%                   epochs are corrected regardless, so the model is
%                   extrapolated outside the range. Leave empty to fit on
%                   everything, which is only safe when no signal is expected.
%                   Default is [].
%            'MinNfit' - Least number of usable epochs inside FitJDRange for a
%                   source to be corrected at all. Default is 200.
%            'MaxCorrection' - Largest correction [pix] any measurement may
%                   receive. A source whose model exceeds it is left alone
%                   rather than moved by an implausible amount. Default is 0.5.
%            'UseWeights' - Weight the fit by the IterFit weights.
%                   Default is true.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - A copy of Obj with the chromatic term subtracted from X and Y.
%          - The correction subtracted from X, [Nepoch, Nsrc].
%          - The correction subtracted from Y.
%          - A structure with the per-source coefficients, the number of epochs
%            each fit used, and which sources were skipped or capped.
% Author : ULTRASAT team (Aug 2026)
% Example: JD2026 = juliandate(datetime(2026,1,1));
%          [ObjC,CX,CY,Info] = ml.detrend.correctSourceDCR(IFsys, Obj, ...
%              'FitJDRange',[-Inf JD2026], 'Verbosity',1);
%
% <Unknown Bugs>: fitting a free chromatic term per source can absorb a real
%   astrometric signal, the parallactic angle being correlated with time within
%   a season. Keep the signal window out of FitJDRange.

    arguments
        IF;
        Obj;
        Args.FitJDRange       = [];
        Args.MinNfit          = 200;
        Args.MaxCorrection    = 0.5;
        Args.UseWeights       = true;
        Args.Verbosity        = 0;
    end

    if ~isfield(IF.Data,'pa') || ~isfield(IF.Data,'secz')
        error('correctSourceDCR:NoGeometry','IF.Data must hold pa and secz');
    end

    [Rx, Ry] = IF.calculateResiduals;
    Pa   = IF.Data.pa(:,1);
    Secz = IF.Data.secz(:,1);
    H    = [ones(IF.Nepoch,1), sin(Pa).*Secz, cos(Pa).*Secz];

    % epochs the fit is allowed to see
    if isempty(Args.FitJDRange)
        FitFlag = true(IF.Nepoch,1);
        fprintf(['correctSourceDCR: fitting on every epoch. A free chromatic term can absorb a real ', ...
                 'signal, so keep the signal window out of FitJDRange when one is expected.\n']);
    else
        FitFlag = IF.JD >= Args.FitJDRange(1) & IF.JD <= Args.FitJDRange(2);
    end

    if Args.UseWeights
        Wes = calculateWes(IF, 'NormalizeWeights', false);
    else
        Wes = ones(size(Rx));
    end

    CorrX = zeros(IF.Nepoch, IF.Nsrc);
    CorrY = zeros(IF.Nepoch, IF.Nsrc);
    Par   = nan(2, IF.Nsrc, 2);          % [sin cos] x source x [X Y]
    Nfit  = zeros(1, IF.Nsrc);
    Capped  = false(1, IF.Nsrc);
    TooFew  = false(1, IF.Nsrc);

    for Isrc = 1:IF.Nsrc
        Ok = FitFlag & isfinite(Rx(:,Isrc)) & isfinite(Ry(:,Isrc)) & isfinite(Wes(:,Isrc)) & Wes(:,Isrc)>0;
        Nfit(Isrc) = sum(Ok);
        if Nfit(Isrc) < Args.MinNfit
            TooFew(Isrc) = true;
            continue
        end
        Wsrc = Wes(Ok,Isrc);
        Px = lscov(H(Ok,:), Rx(Ok,Isrc), Wsrc);
        Py = lscov(H(Ok,:), Ry(Ok,Isrc), Wsrc);
        % the constant is fitted but not applied, so the mean position stands
        Cx = H(:,2:3)*Px(2:3);
        Cy = H(:,2:3)*Py(2:3);
        if max(max(abs(Cx)), max(abs(Cy))) > Args.MaxCorrection
            Capped(Isrc) = true;
            continue
        end
        CorrX(:,Isrc)  = Cx;
        CorrY(:,Isrc)  = Cy;
        Par(:,Isrc,1)  = Px(2:3);
        Par(:,Isrc,2)  = Py(2:3);
    end

    Info = struct('Par',Par, 'Nfit',Nfit, 'TooFew',TooFew, 'Capped',Capped, ...
                  'NcorrectedSrc',sum(~TooFew & ~Capped), 'NepochFit',sum(FitFlag), ...
                  'MaxCorrectionApplied',max(max(abs(CorrX(:))), max(abs(CorrY(:)))));

    report(Args.Verbosity, ['correctSourceDCR: fitted on %d of %d epochs; corrected %d of %d sources ', ...
                            '(%d with too few epochs, %d capped); largest correction %.4g pix\n'], ...
           Info.NepochFit, IF.Nepoch, Info.NcorrectedSrc, IF.Nsrc, sum(TooFew), sum(Capped), ...
           Info.MaxCorrectionApplied);

    ObjOut = Obj.copy();
    ObjOut.Data.X = ObjOut.Data.X - CorrX;
    ObjOut.Data.Y = ObjOut.Data.Y - CorrY;
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
