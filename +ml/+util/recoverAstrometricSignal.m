function [AmpX, AmpY, Info] = recoverAstrometricSignal(IF, Template, Args)
% Matched-filter amplitude of a template in the residuals of a fit.
% Measures how much of an excursion of a given shape is present in each
% source's residuals, in mas, which is what recovers an injection made by
% ml.util.injectAstrometricSignal and what sets an honest upper limit on a real
% signal.
% The template must first be made orthogonal to the model the fit has already
% removed. A position, a proper motion and a per-epoch transformation are taken
% out before the residuals are formed, and a slow excursion lies mostly inside
% that space: for a 90 day profile in a 142 day season the raw template
% recovers 9% of what was injected and the orthogonalised one 90%. Projecting
% onto the raw template does not bias the answer, since the noise is measured
% the same way, but it throws away most of the sensitivity and makes the limit
% about three times worse than it needs to be.
% Input  : - A fitted IterFit object.
%          - The template, either a function handle of JD or a vector with one
%            entry per epoch of the fit. Its scale is irrelevant; the amplitude
%            is returned in the units of the template's peak.
%          * ...,key,val,...
%            'Orthogonalise' - Take the model out of the template before
%                   projecting. Default is true. False reproduces the naive
%                   estimator and is there only for comparison.
%            'ModelColumns' - Design columns to orthogonalise against, as a
%                   matrix with one row per epoch. Default is [], meaning a
%                   constant and a linear term in time, which is what the
%                   source model holds.
%            'RejectOutliers' - Apply the same outlier mask calculateRstd uses.
%                   Default is true.
%            'PixScale' - mas per pixel. Default is 400.
% Output : - Amplitude along X, one entry per source, in mas.
%          - The same along Y.
%          - A structure with the template actually used, its norm relative to
%            the raw one, and the number of epochs each source contributed.
% Author : ULTRASAT team (Aug 2026)
% <Unknown Bugs>: orthogonalising against a constant and a line does not cover
%   the per-epoch transformation or the reweighting, which between them absorb
%   about a tenth more. Recovery is therefore near 90% rather than exactly 100%,
%   and an injection is the only way to know the figure for a given run.
% Example: [Ax,Ay] = ml.util.recoverAstrometricSignal(IFsys, Inf.Profile);
%          Surv = (Ax(Ind) - AxNull(Ind))./20;

    arguments
        IF
        Template
        Args.Orthogonalise  = true;
        Args.ModelColumns   = [];
        Args.RejectOutliers = true;
        Args.PixScale       = 400;
    end

    JD = IF.JD(:);
    if isa(Template,'function_handle')
        G = Template(JD);
    else
        G = Template(:);
    end
    G = G(:);
    if numel(G)~=IF.Nepoch
        error('ml:util:recoverAstrometricSignal:BadTemplate', ...
              'The template holds %d entries for a fit of %d epochs', numel(G), IF.Nepoch);
    end

    if isempty(Args.ModelColumns)
        H = [ones(IF.Nepoch,1), (JD-IF.JD0)./365.25];
    else
        H = Args.ModelColumns;
    end

    [Rx, Ry] = IF.calculateResiduals;
    W        = IF.calculateWes;
    if Args.RejectOutliers
        Bad = isoutlier(Rx,'movmedian',30,"ThresholdFactor",1.5,'SamplePoints',JD) ...
            | isoutlier(Ry,'movmedian',30,"ThresholdFactor",1.5,'SamplePoints',JD) ...
            | isoutlier(IF.Data.MAG_PSF,'movmedian',30,"ThresholdFactor",1.5,'SamplePoints',JD);
        W(Bad) = 0;
    end
    W(~isfinite(Rx) | ~isfinite(Ry)) = 0;
    Rx(~isfinite(Rx)) = 0;
    Ry(~isfinite(Ry)) = 0;

    Nsrc = IF.Nsrc;
    AmpX = nan(Nsrc,1);
    AmpY = nan(Nsrc,1);
    Nep  = zeros(Nsrc,1);
    Rel  = nan(Nsrc,1);
    % Each source carries its own weights, so the template has to be made
    % orthogonal to the model under each source's own weighting.
    for Isrc = 1:Nsrc
        Wc = W(:,Isrc);
        Nep(Isrc) = sum(Wc>0);
        if sum(Wc) <= 0
            continue
        end
        if Args.Orthogonalise
            Sw = sqrt(Wc);
            Gp = G - H*((H.*Sw)\(G.*Sw));
        else
            Gp = G;
        end
        Den = sum(Wc.*Gp.^2);
        if Den <= 0
            continue
        end
        Rel(Isrc)  = sqrt(Den./max(eps, sum(Wc.*G.^2)));
        AmpX(Isrc) = Args.PixScale.*sum(Wc.*Gp.*Rx(:,Isrc))./Den;
        AmpY(Isrc) = Args.PixScale.*sum(Wc.*Gp.*Ry(:,Isrc))./Den;
    end

    Info = struct('Template',G, 'ModelColumns',H, 'Orthogonalised',Args.Orthogonalise, ...
                  'RelativeNorm',Rel, 'MedianRelativeNorm',median(Rel,'omitnan'), ...
                  'Nepoch',Nep, 'PixScale',Args.PixScale);
end
