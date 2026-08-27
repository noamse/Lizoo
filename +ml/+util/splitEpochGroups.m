function [GroupId, Info] = splitEpochGroups(JD, Args)
% Split epochs into observing groups, by the gaps between them.
% A target is visible for part of the year, so its epochs fall into blocks
% separated by long gaps. Rather than assuming those blocks line up with
% calendar years, which holds for a bulge field observed from one site and not
% in general, the split is made wherever the gap between consecutive epochs
% exceeds MinGap. Groups too small to be useful are dropped.
% Input  : - Vector of JD. Need not be sorted.
%          * ...,key,val,...
%            'MinGap' - A gap longer than this [days] starts a new group.
%                   Default is 30.
%            'MinEpochs' - Groups holding fewer epochs than this are dropped,
%                   receiving group 0. Default is 500.
%            'MinNights' - Groups spanning fewer distinct nights than this are
%                   dropped. Default is 50.
%            'MinSpan' - Groups spanning fewer than this many days are dropped.
%                   Default is 60.
%            'Groups' - An explicit group number per epoch, used instead of
%                   splitting on gaps. Zero or NaN marks an epoch as excluded.
%                   Default is [].
% Output : - Group number per epoch, in the order of the input JD, numbered
%            from 1 in time order. Zero marks an excluded epoch.
%          - A structure with one row per surviving group: its number, epoch
%            count, night count, first and last JD, span and mean JD, plus a
%            record of the groups that were dropped and why.
% Author : ULTRASAT team (Aug 2026)
% Example: [G,Info] = ml.util.splitEpochGroups(IFsys.JD, 'Verbosity',1);
%          % KMT-2026-BLG-0058 gives 10 usable groups: 2016-2019, 2021-2026,
%          % 2020 being dropped as a 110 epoch stub cut short in March.

    arguments
        JD
        Args.MinGap      = 30;
        Args.MinEpochs   = 500;
        Args.MinNights   = 50;
        Args.MinSpan     = 60;
        Args.Groups      = [];
    end

    JD  = JD(:);
    Nep = numel(JD);

    if ~isempty(Args.Groups)
        Raw = Args.Groups(:);
        Raw(isnan(Raw)) = 0;
    else
        [Sorted, Isort] = sort(JD);
        RawSorted = cumsum([1; diff(Sorted) > Args.MinGap]);
        Raw = zeros(Nep,1);
        Raw(Isort) = RawSorted;
    end

    % --- drop the groups too small to stand on their own ------------------
    Uraw    = unique(Raw(Raw>0));
    Keep    = false(numel(Uraw),1);
    Dropped = struct('Group',{},'Nepoch',{},'Nnights',{},'Span',{},'Reason',{});
    for K = 1:numel(Uraw)
        Sel   = Raw==Uraw(K);
        Nn    = numel(unique(floor(JD(Sel)-0.5)));
        Span  = max(JD(Sel)) - min(JD(Sel));
        Why   = '';
        if sum(Sel) < Args.MinEpochs
            Why = sprintf('%d epochs < MinEpochs %d', sum(Sel), Args.MinEpochs);
        elseif Nn < Args.MinNights
            Why = sprintf('%d nights < MinNights %d', Nn, Args.MinNights);
        elseif Span < Args.MinSpan
            Why = sprintf('%.0f day span < MinSpan %g', Span, Args.MinSpan);
        end
        if isempty(Why)
            Keep(K) = true;
        else
            Dropped(end+1) = struct('Group',Uraw(K), 'Nepoch',sum(Sel), ...
                                    'Nnights',Nn, 'Span',Span, 'Reason',Why); %#ok<AGROW>
        end
    end

    % --- renumber the survivors from 1, in time order ---------------------
    GroupId  = zeros(Nep,1);
    Survivor = Uraw(Keep);
    MeanJD   = arrayfun(@(G) mean(JD(Raw==G)), Survivor);
    [~, Iord] = sort(MeanJD);
    Survivor  = Survivor(Iord);
    for K = 1:numel(Survivor)
        GroupId(Raw==Survivor(K)) = K;
    end

    Ngroup = numel(Survivor);
    Info = struct('Ngroup',Ngroup, 'Ndropped',numel(Dropped), 'Dropped',Dropped, ...
                  'MinGap',Args.MinGap, 'Nepoch',zeros(Ngroup,1), 'Nnights',zeros(Ngroup,1), ...
                  'FirstJD',zeros(Ngroup,1), 'LastJD',zeros(Ngroup,1), ...
                  'MeanJD',zeros(Ngroup,1), 'Span',zeros(Ngroup,1));
    for K = 1:Ngroup
        Sel = GroupId==K;
        Info.Nepoch(K)  = sum(Sel);
        Info.Nnights(K) = numel(unique(floor(JD(Sel)-0.5)));
        Info.FirstJD(K) = min(JD(Sel));
        Info.LastJD(K)  = max(JD(Sel));
        Info.MeanJD(K)  = mean(JD(Sel));
        Info.Span(K)    = Info.LastJD(K) - Info.FirstJD(K);
    end
    Info.NepochExcluded = sum(GroupId==0);
end
