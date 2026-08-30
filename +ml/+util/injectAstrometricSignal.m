function [Result, Info] = injectAstrometricSignal(MS, SrcInd, Args)
% Add an astrometric excursion of known shape and size to chosen sources.
% Every upper limit on an astrometric signal rests on how much of that signal
% the detrending leaves behind, and that is not something to be argued from the
% residual scatter: a fit removes a position, a motion and a transformation per
% epoch, and how much of a given excursion falls into that space depends on the
% shape of the excursion and on the length of the run. Injecting one of known
% size and recovering it with ml.util.recoverAstrometricSignal measures it.
% Input  : - A MatchedSources or MMS object.
%          - Indices of the sources to inject into, in that object's own source
%            numbering.
%          * ...,key,val,...
%            'Amplitude' - Peak size of the excursion [mas]. Default is 20.
%            'Shape' - 'gaussian', 'tophat', 'ramp', or a function handle of JD
%                   returning the profile normalised to a peak of one.
%                   Default is 'gaussian'.
%            'FWHM' - Width of the profile [days], the full width for 'tophat'.
%                   Default is 90.
%            'PeakJD' - Where the profile peaks. Default is the median JD.
%            'Axis' - 'x', 'y' or 'both', or one entry per source with 1 for x
%                   and 2 for y. Default is 'both'.
%            'Sign' - Multiplies the amplitude, scalar or one per source.
%                   Injecting several sources with one shape and alternating
%                   signs keeps the ensemble free of a common mode. Default is 1.
%            'PixScale' - mas per pixel, since X and Y are in pixels.
%                   Default is 400.
%            'CreateNewObj' - Work on a copy. Default is true.
% Output : - The object with the excursion added to X and Y.
%          - A structure holding the profile as a function of JD, its values at
%            the object's own epochs, and a record of what went where.
% Author : ULTRASAT team (Aug 2026)
% <Unknown Bugs>: the excursion is added to the registered pixel coordinates,
%   so it is injected downstream of the registration and of anything else
%   KMT_pipelineI did. A real signal would have passed through those steps and
%   might have been attenuated by them; this measures only what the detrending
%   chain does.
% Example: [MSi,Inf] = ml.util.injectAstrometricSignal(MS, 420, 'Amplitude',20);
%          % then fit MSi and call ml.util.recoverAstrometricSignal

    arguments
        MS
        SrcInd
        Args.Amplitude     = 20;
        Args.Shape         = 'gaussian';
        Args.FWHM          = 90;
        Args.PeakJD        = [];
        Args.Axis          = 'both';
        Args.Sign          = 1;
        Args.PixScale      = 400;
        Args.CreateNewObj  = true;
    end

    if Args.CreateNewObj
        Result = MS.copy();
    else
        Result = MS;
    end

    JD = Result.JD(:);
    if isempty(Args.PeakJD)
        Args.PeakJD = median(JD, 'omitnan');
    end

    if isa(Args.Shape,'function_handle')
        Profile = Args.Shape;
    else
        switch lower(Args.Shape)
            case 'gaussian'
                Sigma   = Args.FWHM./2.3548;
                Profile = @(T) exp(-0.5.*((T-Args.PeakJD)./Sigma).^2);
            case 'tophat'
                Profile = @(T) double(abs(T-Args.PeakJD) <= 0.5.*Args.FWHM);
            case 'ramp'
                Profile = @(T) (T-min(JD))./max(eps, max(JD)-min(JD));
            otherwise
                error('ml:util:injectAstrometricSignal:BadShape', ...
                      'Shape must be gaussian, tophat, ramp or a function handle, not %s', Args.Shape);
        end
    end

    SrcInd = SrcInd(:);
    Nsrc   = numel(SrcInd);
    if any(SrcInd<1 | SrcInd>Result.Nsrc)
        error('ml:util:injectAstrometricSignal:BadSrcInd', ...
              'Source indices must lie in 1..%d', Result.Nsrc);
    end
    Sgn = Args.Sign(:);
    if isscalar(Sgn)
        Sgn = repmat(Sgn, Nsrc, 1);
    elseif numel(Sgn)~=Nsrc
        error('ml:util:injectAstrometricSignal:BadSign', ...
              'Sign must be scalar or hold one entry per source (%d)', Nsrc);
    end
    if ischar(Args.Axis) || isstring(Args.Axis)
        switch lower(string(Args.Axis))
            case "x",    DoX = true(Nsrc,1);  DoY = false(Nsrc,1);
            case "y",    DoX = false(Nsrc,1); DoY = true(Nsrc,1);
            case "both", DoX = true(Nsrc,1);  DoY = true(Nsrc,1);
            otherwise
                error('ml:util:injectAstrometricSignal:BadAxis', ...
                      'Axis must be x, y, both, or one entry per source');
        end
    else
        Ax  = Args.Axis(:);
        DoX = Ax==1;
        DoY = Ax==2;
    end

    G   = Profile(JD);
    Amp = Args.Amplitude./Args.PixScale;      % X and Y are in pixels
    for I = 1:Nsrc
        if DoX(I)
            Result.Data.X(:,SrcInd(I)) = Result.Data.X(:,SrcInd(I)) + Sgn(I).*Amp.*G;
        end
        if DoY(I)
            Result.Data.Y(:,SrcInd(I)) = Result.Data.Y(:,SrcInd(I)) + Sgn(I).*Amp.*G;
        end
    end

    Info = struct('Profile',Profile, 'Template',G, 'JD',JD, 'PeakJD',Args.PeakJD, ...
                  'Amplitude',Args.Amplitude, 'AmplitudePix',Amp, 'Shape',Args.Shape, ...
                  'FWHM',Args.FWHM, 'SrcInd',SrcInd, 'Sign',Sgn, ...
                  'InjectedX',SrcInd(DoX), 'InjectedY',SrcInd(DoY), ...
                  'PixScale',Args.PixScale);
end
