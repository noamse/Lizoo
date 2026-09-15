% Full chain on the re-registered v3 MSc files.
addpath('~/matlab/Lizoo');
F = getenv('FIELD');
In  = sprintf('/home/ocs/KMTdata/Results16_26_v3/KMT_260058_%s_MSc.mat', F);
PM  = sprintf('/home/ocs/KMTdata/Results/v3_JointPM1625/IFsys1625_260058_CTIO_%s.mat', F);
Sys = sprintf('/home/ocs/KMTdata/Results/v3_DecadeSysRem/SysCor_260058_CTIO_%s.mat', F);
Fin = sprintf('/home/ocs/KMTdata/Results/v3_Final/IFfinal_260058_CTIO_%s.mat', F);

L=load(In,'MSc'); MS=L.MSc; clear L
[G,GI]=ml.util.splitEpochGroups(MS.JD);
fprintf('%s: %d epochs, %d sources, %d seasons | I_ogle coverage %.1f%%\n', F, MS.Nepoch, MS.Nsrc, ...
    GI.Ngroup, 100*mean(isfinite(MS.SrcData.I_ogle(:))));

% --- 1. proper motions, event season excluded ---------------------------
T0=tic;
Sub=MS.selectByEpoch(find(ismember(G,1:GI.Ngroup-1)));
[IFsys,~,~,Info]=ml.scripts.runIterDetrendMSc(Sub,'Verbosity',0);
save(PM,'IFsys','Info','-v7.3');
M=IFsys.medianFieldSource({'MAG_PSF'}); M=M(:); Br=M<17&isfinite(M);
[rx,ry]=IFsys.calculateRstd;
fprintf('%s STEP1 joint 2016-2025: %.1f min | bright %.3f/%.3f | colour %s %.1f%%\n', F, toc(T0)/60, ...
    median(rx(Br),'omitnan'), median(ry(Br),'omitnan'), Info.Colour.Mode, 100*Info.Colour.FracValid);
clear IFsys Sub

% --- 2. shared decade SysRem correction ---------------------------------
% There is deliberately no per-season fitting step. One was run here until
% 2026-09-09, writing ~/KMTdata/Results/v3_Seasons, but the final solution below
% takes the full MS and never read it, so it cost ~40 min per field and fed
% nothing. Solving once over the decade is also what cut the target's
% season-to-season scatter from 5.32 to 2.28 mas in the first place.
T0=tic;
[~,~,~,I2]=ml.scripts.runIterDetrendMSc(MS,'FixedPM',PM,'Verbosity',0);
SysCorX=I2.SysCorX; SysCorY=I2.SysCorY; SrcInd=I2.SrcInd(:); JD=MS.JD(I2.EpochInd); JD=JD(:);
save(Sys,'SysCorX','SysCorY','JD','SrcInd','-v7.3');
fprintf('%s STEP2 decade SysRem: %.1f min | correction %dx%d\n', F, toc(T0)/60, size(SysCorX,1), size(SysCorX,2));
clear I2 SysCorX SysCorY

% --- 3. final global solution, motions fitted ---------------------------
T0=tic;
[IFsys,~,IFsysB,Info]=ml.scripts.runIterDetrendMSc(MS, ...
    'PixPhase',true,'NiterWeightsBeforeSys',15,'NiterWeightsAfterSys',6,'Verbosity',0);
M=IFsys.medianFieldSource({'MAG_PSF'}); M=M(:); Br=M<17&isfinite(M);
[bx,by]=IFsysB.calculateRstd; [fx,fy]=IFsys.calculateRstd;
Ie=IFsys.findClosestSource([150 150]);
Track=cellfun(@(v) median(v(Br),'omitnan'), IFsys.RMSTrack); n=numel(Track);
Config=struct('Input',In,'Field',F,'Note','v3: OGLE re-registered per field','Created',string(datetime('now')));
save(Fin,'IFsys','Info','Config','-v7.3');
fprintf('%s STEP3 final: %.1f min | pass1 %.3f/%.3f -> %.3f/%.3f | target %.2f/%.2f (I_ogle %.2f) | conv %.4f\n', ...
    F, toc(T0)/60, median(bx(Br),'omitnan'),median(by(Br),'omitnan'), ...
    median(fx(Br),'omitnan'),median(fy(Br),'omitnan'), fx(Ie),fy(Ie), ...
    Info.SrcData.I_ogle(Ie), max(Track(n-2:n))-min(Track(n-2:n)));
fprintf('V3CHAIN DONE %s\n', F);
