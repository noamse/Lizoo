addpath('/home/ocs/matlab/Lizoo');
addpath('/home/ocs/matlab/Lizoo/OGLEdata/OB260058');   % KMT_pipelineI does load('OB160058.mat')
OutDir = '/home/ocs/KMTdata/Results16_26_v2';
if ~isfolder(OutDir), mkdir(OutDir); end
Fid = fopen(fullfile(OutDir,'.wtest'),'w'); assert(Fid>0,'output dir not writable'); fclose(Fid);
delete(fullfile(OutDir,'.wtest'));

Fields = {'BLG41','BLG01'};
for K = 1:numel(Fields)
    Fl   = Fields{K};
    Temp = sprintf('/bigdata3/projects/KMTdata/Images/260058/KB260058_20*_CTIO_I_%s/RAW/*I*.fits', Fl);
    fprintf('\n================ %s start %s ================\n', Fl, datestr(now,'yyyy-mm-dd HH:MM:SS'));
    T0 = tic;
    try
        [~, JD, MSc, ~] = KMT_pipelineI([], 'TempName', Temp);
        Elapsed = toc(T0);
        fprintf('%s: KMT_pipelineI done in %.0f s (%.2f h)\n', Fl, Elapsed, Elapsed/3600);
        fprintf('%s: MSc Nepoch=%d Nsrc=%d fracNaN=%.3f\n', Fl, MSc.Nepoch, MSc.Nsrc, mean(isnan(MSc.Data.X(:))));
        UD = MSc.UserData;
        if isstruct(UD)
            fprintf('%s: UserData present - ShiftX %d, FlagGoodEpoch %d/%d kept, FlagGoodSrc %d/%d kept\n', Fl, ...
                numel(UD.ShiftX), sum(UD.FlagGoodEpoch), numel(UD.FlagGoodEpoch), ...
                sum(UD.FlagGoodSrc), numel(UD.FlagGoodSrc));
        else
            fprintf('%s: WARNING UserData absent\n', Fl);
        end
        NI = sum(~isnan(MSc.SrcData.I_ogle));
        fprintf('%s: I_ogle coverage %d/%d (%.1f%%)\n', Fl, NI, MSc.Nsrc, 100*NI/MSc.Nsrc);
        Out = fullfile(OutDir, sprintf('KMT_260058_%s_MSc.mat', Fl));
        Ts = tic; save(Out, 'MSc', 'JD', '-v7.3');
        D = dir(Out);
        fprintf('%s: saved %s (%.2f GB) in %.0f s\n', Fl, Out, D.bytes/1e9, toc(Ts));
    catch ME
        fprintf('%s: FAILED after %.0f s: %s | %s\n', Fl, toc(T0), ME.identifier, ME.message);
        fprintf('%s\n', getReport(ME,'extended','hyperlinks','off'));
    end
    clear MSc JD
end
fprintf('\n=== all fields done %s ===\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
