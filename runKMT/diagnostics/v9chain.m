% v9: as v8, but SysRem is refitted SEPARATELY IN EACH SEASON and the resulting
% per-season corrections are stitched into one decade-long correction, which the
% final global fit then applies. The source set, selection and every other
% setting are identical to v8, so v8 vs v9 isolates the SysRem cadence.
addpath('~/matlab/Lizoo');
F = getenv('FIELD');
In  = sprintf('/home/ocs/KMTdata/Results16_26_v3/KMT_260058_%s_MSc.mat', F);
Sys = sprintf('/home/ocs/KMTdata/Results/v9/SysCor_%s.mat', F);
Fin = sprintf('/home/ocs/KMTdata/Results/v9/IFfinal_%s.mat', F);
P=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',F));
L=load(In,'MSc'); MS0=L.MSc; clear L
MS=MS0.selectBySrcIndex(P.SrcIdx); clear MS0
[G,GI]=ml.util.splitEpochGroups(MS.JD);
fprintf('%s: %d sources, %d epochs, %d seasons\n', F, MS.Nsrc, MS.Nepoch, GI.Ngroup);

% --- per-season SysRem ----------------------------------------------------
T0=tic;
SysCorX=zeros(MS.Nepoch,MS.Nsrc); SysCorY=zeros(MS.Nepoch,MS.Nsrc);
Hit=false(MS.Nepoch,1);
for b=1:GI.Ngroup
  ei=find(G==b);
  if numel(ei)<200, fprintf('  season %2d: only %d epochs, skipped\n', b, numel(ei)); continue; end
  Sub=MS.selectByEpoch(ei);
  try
    [~,~,~,I2]=ml.scripts.runIterDetrendMSc(Sub,'Verbosity',0);
  catch E
    fprintf('  season %2d FAILED: %s\n', b, E.message); continue;
  end
  if ~isfield(I2,'SysCorX') || isempty(I2.SysCorX)
    fprintf('  season %2d: no SysRem correction returned\n', b); continue;
  end
  % I2 indexes the season subset; map back onto the full arrays
  eg = ei(I2.EpochInd(:));
  sg = I2.SrcInd(:);
  ok = sg>=1 & sg<=MS.Nsrc;
  SysCorX(eg, sg(ok)) = I2.SysCorX(:, ok);
  SysCorY(eg, sg(ok)) = I2.SysCorY(:, ok);
  Hit(eg)=true;
  fprintf('  season %2d: %d epochs, correction %dx%d, rms %.2f/%.2f mas\n', b, numel(ei), ...
      size(I2.SysCorX,1), size(I2.SysCorX,2), ...
      400*std(I2.SysCorX(isfinite(I2.SysCorX))), 400*std(I2.SysCorY(isfinite(I2.SysCorY))));
  clear I2 Sub
end
JD=MS.JD(:); SrcInd=(1:MS.Nsrc)';
save(Sys,'SysCorX','SysCorY','JD','SrcInd','-v7.3');
fprintf('%s STEP A per-season SysRem: %.1f min | %d of %d epochs covered | overall rms %.2f/%.2f mas\n', ...
    F, toc(T0)/60, sum(Hit), MS.Nepoch, 400*std(SysCorX(:)), 400*std(SysCorY(:)));

% --- final global fit, applying the stitched correction --------------------
T0=tic; Pre=load(Sys);
[IFsys,~,IFsysB,Info]=ml.scripts.runIterDetrendMSc(MS, ...
    'SysRemCorrection',Pre,'PixPhase',true, ...
    'NiterWeightsBeforeSys',15,'NiterWeightsAfterSys',6,'Verbosity',0);
[fx,fy]=IFsys.calculateRstd;
Ie=find(P.SrcIdx(Info.SrcInd(:))==P.TargetMsc);
if isempty(Ie), error('v9chain:TargetMissing','target MSc %d absent', P.TargetMsc); end
M=IFsys.medianFieldSource({'MAG_PSF'}); M=M(:); Br=M<17&isfinite(M);
Tr=cellfun(@(v) median(v(Br),'omitnan'), IFsys.RMSTrack); n=numel(Tr);
Config=struct('Input',In,'Field',F,'Note','v9: per-season SysRem, otherwise identical to v8', ...
              'Created',string(datetime('now')));
save(Fin,'IFsys','Info','Config','-v7.3');
fprintf('%s STEP B: %.1f min | bright(<17) %.3f/%.3f | target %.2f/%.2f at (%.1f,%.1f) | conv %.4f | Nsrc %d\n', ...
    F, toc(T0)/60, median(fx(Br),'omitnan'), median(fy(Br),'omitnan'), fx(Ie), fy(Ie), ...
    IFsys.ParS(1,Ie), IFsys.ParS(2,Ie), max(Tr(n-2:n))-min(Tr(n-2:n)), IFsys.Nsrc);
fprintf('V9 DONE %s\n', F);
