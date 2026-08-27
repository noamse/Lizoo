function [Res, Info] = detrendGroupScaling(MS, Args)
% Detrend blocks of consecutive observing groups, to see how the fit improves
% as more of them are combined.
% The epochs are split into groups by ml.util.splitEpochGroups, which for a
% bulge field observed from one site gives one group per observing season. For
% every requested block length the groups are taken in consecutive,
% non-overlapping runs, each block is detrended on its own, and the resulting
% precision is recorded. Comparing across block lengths shows what combining
% two, three or ten seasons is actually worth.
% Input  : - A MatchedSources object as produced by KMT_pipelineI, or the path
%            of a .mat file holding one.
%          * ...,key,val,...
%            'NGroups' - Block lengths to try, in groups. Default is
%                   [1 2 3 5 10].
%            'MaxBlocksPerN' - At most this many blocks are fitted at each
%                   length, taken from the start. Default is 3.
%            'splitEpochGroupsArgs' - Cell array for ml.util.splitEpochGroups.
%                   Default is {}.
%            'runIterDetrendMScArgs' - Cell array for
%                   ml.scripts.runIterDetrendMSc. Default is {}.
%            'TargetXY' - Position of the source to follow, in pixels.
%                   Default is [150 150].
%            'SaveDir' - Directory for one .mat per block. Nothing is saved
%                   when empty. Default is ''.
%            'Verbosity' - 0 is silent. Default is 0.
% Output : - A table with one row per block: the block length, which groups it
%            covers, its epoch count and time span, the residual scatter of the
%            bright sources and of the target, and the target's fitted position
%            evaluated at the block's own mean epoch.
%          - A structure with the grouping, and for each block length the
%            scatter of the target position about a straight line in time,
%            which is the precision a signal would have to exceed.
% Author : ULTRASAT team (Aug 2026)
% <Unknown Bugs>: a single group spans well under a year, so the proper motion
%   is barely constrained within it and absorbs part of any real motion.
%   Positions are therefore reported at each block's own mean epoch, where they
%   are well determined, rather than extrapolated to a common one.
% Example: Res = ml.scripts.detrendGroupScaling( ...
%              '~/KMTdata/Results16_26_v2/KMT_260058_BLG41_MSc.mat', ...
%              'NGroups',[1 2 5 10], 'Verbosity',1);

    arguments
        MS
        Args.NGroups                = [1 2 3 5 10];
        Args.MaxBlocksPerN          = 3;
        Args.splitEpochGroupsArgs   = {};
        Args.runIterDetrendMScArgs  = {};
        Args.TargetXY               = [150 150];
        Args.SaveDir                = '';
        Args.Verbosity              = 0;
    end

    if ischar(MS) || isstring(MS)
        Loaded = load(char(MS));
        Vars   = fieldnames(Loaded);
        IsMS   = cellfun(@(V) isa(Loaded.(V),'MatchedSources'), Vars);
        MS     = Loaded.(Vars{find(IsMS,1)});
    end
    if ~isempty(Args.SaveDir) && ~isfolder(Args.SaveDir)
        mkdir(Args.SaveDir);
    end

    [GroupId, GInfo] = ml.util.splitEpochGroups(MS.JD, Args.splitEpochGroupsArgs{:});
    Info.Grouping = GInfo;
    report(Args.Verbosity, '%d groups, %d epochs excluded\n', GInfo.Ngroup, GInfo.NepochExcluded);
    for K = 1:GInfo.Ngroup
        report(Args.Verbosity, '  group %2d: %5d epochs, %3d nights, %s to %s\n', K, ...
            GInfo.Nepoch(K), GInfo.Nnights(K), ...
            string(datetime(GInfo.FirstJD(K),'convertfrom','juliandate'),'yyyy-MM-dd'), ...
            string(datetime(GInfo.LastJD(K),'convertfrom','juliandate'),'yyyy-MM-dd'));
    end

    Rows = {};
    for N = Args.NGroups(:).'
        if N > GInfo.Ngroup
            report(Args.Verbosity, 'skipping N=%d, only %d groups exist\n', N, GInfo.Ngroup);
            continue
        end
        Starts = 1:N:(GInfo.Ngroup-N+1);
        Starts = Starts(1:min(numel(Starts), Args.MaxBlocksPerN));
        for Ib = 1:numel(Starts)
            First = Starts(Ib);
            Sel   = ismember(GroupId, First:(First+N-1));
            report(Args.Verbosity, '\n--- N=%d block %d: groups %d-%d, %d epochs ---\n', ...
                   N, Ib, First, First+N-1, sum(Sel));
            T0 = tic;
            try
                Sub = MS.selectByEpoch(find(Sel));
                [IF, ~, ~, DInfo] = ml.scripts.runIterDetrendMSc(Sub, ...
                    'Verbosity',Args.Verbosity, Args.runIterDetrendMScArgs{:});
                Ie   = IF.findClosestSource(Args.TargetXY);
                Mg   = IF.medianFieldSource({'MAG_PSF'});
                Br   = Mg<17 & isfinite(Mg);
                [RX,RY] = IF.calculateRstd;
                MeanJD  = mean(IF.JD);
                % the position where it is measured, not extrapolated to JD0
                Dt   = (MeanJD - IF.JD0)./365.25;
                PosX = IF.ParS(1,Ie) + IF.ParS(3,Ie).*Dt;
                PosY = IF.ParS(2,Ie) + IF.ParS(4,Ie).*Dt;
                Rows{end+1} = {N, Ib, First, First+N-1, sum(Sel), IF.Nsrc, ...
                               max(IF.JD)-min(IF.JD), MeanJD, ...
                               median(RX(Br),'omitnan'), median(RY(Br),'omitnan'), ...
                               RX(Ie), RY(Ie), 400*PosX, 400*PosY, toc(T0)/60}; %#ok<AGROW>
                report(Args.Verbosity, 'N=%d block %d: rstd bright %.2f/%.2f, target %.1f/%.1f mas, %.1f min\n', ...
                    N, Ib, median(RX(Br),'omitnan'), median(RY(Br),'omitnan'), RX(Ie), RY(Ie), toc(T0)/60);
                if ~isempty(Args.SaveDir)
                    IFsys = IF; BlockInfo = DInfo;
                    save(fullfile(Args.SaveDir, sprintf('IFsys_N%02d_block%02d.mat',N,Ib)), ...
                         'IFsys','BlockInfo','-v7.3');
                end
                clear IF Sub
            catch ME
                fprintf('N=%d block %d FAILED after %.1f min: %s | %s\n', N, Ib, toc(T0)/60, ME.identifier, ME.message);
            end
        end
    end

    if isempty(Rows)
        Res = table();
        Info.Scaling = table();
        return
    end
    Res = cell2table(vertcat(Rows{:}), 'VariableNames', ...
        {'N','Block','FirstGroup','LastGroup','Nepoch','Nsrc','SpanDays','MeanJD', ...
         'RstdBrightX','RstdBrightY','TargetRstdX','TargetRstdY','TargetPosX','TargetPosY','Minutes'});

    % --- how the target position scatters at each block length ------------
    Un = unique(Res.N);
    Sc = nan(numel(Un),4);
    for K = 1:numel(Un)
        Sub = Res(Res.N==Un(K),:);
        Sc(K,1) = Un(K);
        Sc(K,2) = height(Sub);
        if height(Sub) >= 3
            % remove a straight line in time, which is the proper motion
            T = (Sub.MeanJD - mean(Sub.MeanJD))./365.25;
            Sc(K,3) = std(Sub.TargetPosX - polyval(polyfit(T,Sub.TargetPosX,1),T));
            Sc(K,4) = std(Sub.TargetPosY - polyval(polyfit(T,Sub.TargetPosY,1),T));
        elseif height(Sub) == 2
            Sc(K,3) = abs(diff(Sub.TargetPosX))/sqrt(2);
            Sc(K,4) = abs(diff(Sub.TargetPosY))/sqrt(2);
        end
    end
    Info.Scaling = array2table(Sc, 'VariableNames', {'N','Nblocks','TargetScatterX','TargetScatterY'});
end


function report(Verbosity, varargin)
    % Print only when verbosity has been switched on
    if Verbosity > 0
        fprintf(varargin{:});
    end
end
