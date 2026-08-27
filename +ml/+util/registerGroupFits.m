function [Pos, Info] = registerGroupFits(IFs, Args)
% Bring separately fitted observing groups onto one common astrometric frame.
% Each group is solved with its own per-epoch transformation, so its fitted
% positions carry an arbitrary origin, orientation and scale. Positions from
% different groups therefore cannot be compared until the frames are tied
% together. This fits one affine per group, from the stars it shares with a
% reference group, and applies it to every source.
% The registration is deliberately fitted from well-measured stars only, and
% the target must be kept out of it: a source included in the registration
% helps define the frame it is then measured against, which suppresses exactly
% the motion one is looking for.
% Input  : - Cell array of IterFit objects, one per group, each already
%            detrended on its own.
%          * ...,key,val,...
%            'SrcInd' - Cell array, one entry per group, giving for each of
%                   that group's sources its index in a common source list.
%                   Info.SrcInd from ml.scripts.runIterDetrendMSc serves.
%                   Without it the groups are matched by position instead.
%                   Default is {}.
%            'RefGroup' - Group whose frame the others are brought onto.
%                   Default is the one holding the most sources.
%            'RegisterFlag' - Cell array of logicals, one per group, marking
%                   the sources the registration may use. Without it the
%                   choice is made by MaxRegisterMag and MinRegisterN.
%                   Default is {}.
%            'MaxRegisterMag' - Only sources brighter than this take part in
%                   the registration. Default is 17.
%            'ExcludeXY' - Position [pix] of a source to keep out of the
%                   registration, typically the target. Empty for none.
%                   Default is [150 150].
%            'ExcludeRadius' - Sources within this distance [pix] of ExcludeXY
%                   are kept out. Default is 3.
%            'MinRegisterN' - Fewest shared stars a group needs to be
%                   registered. Default is 20.
%            'Order' - 'shift', 'affine' or 'quadratic'. Default is 'affine'.
%            'Niter' - Rejection iterations on the registration fit.
%                   Default is 3.
%            'SigmaClip' - Rejection threshold, in robust sigma. Default is 3.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - A structure with, for each group, the registered position of every
%            source in the common list, the group's mean JD, and which sources
%            it actually measured.
%          - A structure recording the transformation fitted for each group,
%            how many stars it used, and the residual of the registration,
%            which is the floor below which cross-group positions cannot be
%            trusted.
% Author : ULTRASAT team (Aug 2026)
% <Unknown Bugs>: an affine registration absorbs any motion shared by the
%   registration stars, so a real signal common to the whole field would be
%   removed with it. That is the intended behaviour for a frame tie, but it
%   means this cannot detect motion of the field as a whole.
% Example: [Pos,Info] = ml.util.registerGroupFits(IFs, 'SrcInd',SrcInds, 'Verbosity',1);

    arguments
        IFs
        Args.SrcInd          = {};
        Args.RefGroup        = [];
        Args.RegisterFlag    = {};
        Args.MaxRegisterMag  = 17;
        Args.ExcludeXY       = [150 150];
        Args.ExcludeRadius   = 3;
        Args.MinRegisterN    = 20;
        Args.Order           = 'affine';
        Args.Niter           = 3;
        Args.SigmaClip       = 3;
        Args.Verbosity       = 0;
    end

    Ngroup = numel(IFs);
    if Ngroup < 2
        error('registerGroupFits:TooFewGroups','At least two groups are needed');
    end

    % --- a common source list ---------------------------------------------
    UseInd = ~isempty(Args.SrcInd);
    if UseInd
        AllInd = unique(cell2mat(cellfun(@(V) V(:), Args.SrcInd(:), 'UniformOutput',false)));
    else
        % fall back on matching the reference group's positions
        AllInd = (1:IFs{1}.Nsrc).';
    end
    Ncommon = numel(AllInd);

    % --- gather each group's positions on the common list ------------------
    X = nan(Ncommon, Ngroup); Y = X; Mag = X; Rstd = X;
    MeanJD = nan(1,Ngroup);
    for K = 1:Ngroup
        IF = IFs{K};
        MeanJD(K) = mean(IF.JD);
        Dt = (MeanJD(K) - IF.JD0)./365.25;
        Px = IF.ParS(1,:).' + IF.ParS(3,:).'.*Dt;      % at the group's own epoch
        Py = IF.ParS(2,:).' + IF.ParS(4,:).'.*Dt;
        Mg = IF.medianFieldSource({'MAG_PSF'});
        [Rx,Ry] = IF.calculateRstd;
        if UseInd
            [Tf, Loc] = ismember(Args.SrcInd{K}(:), AllInd);
            X(Loc(Tf),K)    = Px(Tf);
            Y(Loc(Tf),K)    = Py(Tf);
            Mag(Loc(Tf),K)  = Mg(Tf);
            Rstd(Loc(Tf),K) = sqrt(Rx(Tf).^2 + Ry(Tf).^2);
        else
            N = min(numel(Px), Ncommon);
            X(1:N,K)=Px(1:N); Y(1:N,K)=Py(1:N); Mag(1:N,K)=Mg(1:N);
            Rstd(1:N,K)=sqrt(Rx(1:N).^2+Ry(1:N).^2);
        end
    end

    % --- which group defines the frame ------------------------------------
    if isempty(Args.RefGroup)
        [~, Iref] = max(sum(isfinite(X),1));
    else
        Iref = Args.RefGroup;
    end

    % --- which sources may define the registration ------------------------
    CanUse = false(Ncommon, Ngroup);
    for K = 1:Ngroup
        if ~isempty(Args.RegisterFlag)
            Ok = false(Ncommon,1);
            Ok(1:min(numel(Args.RegisterFlag{K}),Ncommon)) = logical(Args.RegisterFlag{K}(1:min(numel(Args.RegisterFlag{K}),Ncommon)));
        else
            Ok = isfinite(Mag(:,K)) & Mag(:,K) < Args.MaxRegisterMag;
        end
        CanUse(:,K) = Ok & isfinite(X(:,K)) & isfinite(Y(:,K));
    end
    % the target must not help define the frame it is measured against
    Excluded = false(Ncommon,1);
    if ~isempty(Args.ExcludeXY)
        Xr = X(:,Iref); Yr = Y(:,Iref);
        Excluded = hypot(Xr-Args.ExcludeXY(1), Yr-Args.ExcludeXY(2)) < Args.ExcludeRadius;
        Excluded(~isfinite(Xr)) = false;
        CanUse(Excluded,:) = false;
    end

    % --- register each group onto the reference ---------------------------
    Xr = nan(Ncommon, Ngroup); Yr = Xr;
    Par = cell(1,Ngroup); Nused = zeros(1,Ngroup); Resid = nan(1,Ngroup);
    for K = 1:Ngroup
        if K == Iref
            Xr(:,K) = X(:,K); Yr(:,K) = Y(:,K);
            Par{K} = 'reference'; Nused(K) = sum(CanUse(:,K)); Resid(K) = 0;
            continue
        end
        Sel = CanUse(:,K) & CanUse(:,Iref);
        if sum(Sel) < Args.MinRegisterN
            report(Args.Verbosity, 'group %d: only %d shared stars, not registered\n', K, sum(Sel));
            continue
        end
        [Px, Py, Keep, Rms] = fitRegistration(X(:,K), Y(:,K), X(:,Iref), Y(:,Iref), Sel, ...
                                              Args.Order, Args.Niter, Args.SigmaClip);
        H = designMat(X(:,K), Y(:,K), Args.Order);
        Xr(:,K) = H*Px;  Yr(:,K) = H*Py;
        Par{K} = [Px, Py]; Nused(K) = sum(Keep); Resid(K) = Rms;
        report(Args.Verbosity, 'group %d: registered on %d stars, residual %.3f pix (%.1f mas)\n', ...
               K, sum(Keep), Rms, 400*Rms);
    end

    Pos = struct('X',Xr, 'Y',Yr, 'XRaw',X, 'YRaw',Y, 'Mag',Mag, 'Rstd',Rstd, ...
                 'MeanJD',MeanJD, 'SrcInd',AllInd, 'Measured',isfinite(X));
    Info = struct('RefGroup',Iref, 'Par',{Par}, 'Nused',Nused, 'ResidPix',Resid, ...
                  'ResidMas',400.*Resid, 'Order',Args.Order, 'Ncommon',Ncommon, ...
                  'NExcluded',sum(Excluded), 'MaxRegisterMag',Args.MaxRegisterMag);
end


% -------------------------------------------------------------------------
function H = designMat(X, Y, Order)
    % Design matrix of the requested transformation
    One = ones(size(X));
    switch lower(Order)
        case 'shift',     H = One;
        case 'affine',    H = [One, X, Y];
        case 'quadratic', H = [One, X, Y, X.^2, Y.^2, X.*Y];
        otherwise, error('registerGroupFits:BadOrder','Unknown Order "%s"', Order);
    end
end


function [Px, Py, Keep, Rms] = fitRegistration(X, Y, Xref, Yref, Sel, Order, Niter, SigmaClip)
    % Least squares with iterative rejection, both coordinates sharing the flag
    Keep = Sel;
    Px = []; Py = []; Rms = NaN;
    for Iter = 1:Niter
        H  = designMat(X(Keep), Y(Keep), Order);
        Px = H\Xref(Keep);
        Py = H\Yref(Keep);
        Dx = Xref(Keep) - H*Px;
        Dy = Yref(Keep) - H*Py;
        D  = hypot(Dx, Dy);
        Rms = tools.math.stat.rstd(D);
        if Iter == Niter || Rms == 0 || ~isfinite(Rms)
            break
        end
        Bad = D > SigmaClip.*Rms;
        if ~any(Bad)
            break
        end
        Idx = find(Keep);
        Keep(Idx(Bad)) = false;
    end
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
