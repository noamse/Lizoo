% v10: as v8, second clipping iteration on the season scale (prep10).
addpath('~/matlab/Lizoo');
F = getenv('FIELD');
In  = sprintf('/home/ocs/KMTdata/Results16_26_v3/KMT_260058_%s_MSc.mat', F);
Sys = sprintf('/home/ocs/KMTdata/Results/v10/SysCor_%s.mat', F);
Fin = sprintf('/home/ocs/KMTdata/Results/v10/IFfinal_%s.mat', F);
P=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep10_%s.mat',F));
L=load(In,'MSc'); MS0=L.MSc; clear L
MS=MS0.selectBySrcIndex(P.SrcIdx); clear MS0
fprintf('%s: %d sources, %d epochs, nothing pinned, no RUWE cut at selection\n', F, MS.Nsrc, MS.Nepoch);
T0=tic;
[~,~,~,I2]=ml.scripts.runIterDetrendMSc(MS,'Verbosity',0);
SysCorX=I2.SysCorX; SysCorY=I2.SysCorY; SrcInd=I2.SrcInd(:); JD=MS.JD(I2.EpochInd); JD=JD(:);
save(Sys,'SysCorX','SysCorY','JD','SrcInd','-v7.3');
fprintf('%s STEP A: %.1f min | correction %dx%d\n', F, toc(T0)/60, size(SysCorX,1), size(SysCorX,2));
clear I2 SysCorX SysCorY
T0=tic; Pre=load(Sys);
[IFsys,~,IFsysB,Info]=ml.scripts.runIterDetrendMSc(MS, ...
    'SysRemCorrection',Pre,'PixPhase',true, ...
    'NiterWeightsBeforeSys',15,'NiterWeightsAfterSys',6,'Verbosity',0);
[fx,fy]=IFsys.calculateRstd;
% Identify the target by its MSc index. findClosestSource assumes the fitted
% positions sit near the pixel grid, which is not guaranteed, and it silently
% returns a neighbour when the target is absent or the gauge has moved.
Ie=find(P.SrcIdx(Info.SrcInd(:))==P.TargetMsc);
if isempty(Ie)
    error('v10chain:TargetMissing','the target (MSc %d) is not in the solution', P.TargetMsc);
end
Ichk=IFsys.findClosestSource([150 150]);
fprintf('%s target: index %d at (%.1f,%.1f); findClosestSource would give %d at (%.1f,%.1f)\n', ...
    F, Ie, IFsys.ParS(1,Ie), IFsys.ParS(2,Ie), Ichk, IFsys.ParS(1,Ichk), IFsys.ParS(2,Ichk));
M=IFsys.medianFieldSource({'MAG_PSF'}); M=M(:); Br=M<17&isfinite(M);
Tr=cellfun(@(v) median(v(Br),'omitnan'), IFsys.RMSTrack); n=numel(Tr);
Config=struct('Input',In,'Field',F,'Note','v10: v8 minus 2-sigma season-offset chi2 outliers', ...
              'Created',string(datetime('now')));
save(Fin,'IFsys','Info','Config','-v7.3');
fprintf('%s STEP B: %.1f min | bright(<17) %.3f/%.3f | target %.2f/%.2f | conv %.4f mas | Nsrc %d\n', ...
    F, toc(T0)/60, median(fx(Br),'omitnan'), median(fy(Br),'omitnan'), fx(Ie), fy(Ie), ...
    max(Tr(n-2:n))-min(Tr(n-2:n)), IFsys.Nsrc);
fprintf('V10 DONE %s\n', F);
