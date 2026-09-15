% v14: as v13 (prep8 calibration set, 6 SysRem components) but with the annual
% term switched off in both steps and both passes. Output goes to Results/v14.
addpath('~/matlab/Lizoo');
F = getenv('FIELD'); NSYS = 6; Out = 'v14';
In  = sprintf('/home/ocs/KMTdata/Results16_26_v3/KMT_260058_%s_MSc.mat', F);
Sys = sprintf('/home/ocs/KMTdata/Results/%s/SysCor_%s.mat', Out, F);
Fin = sprintf('/home/ocs/KMTdata/Results/%s/IFfinal_%s.mat', Out, F);
P=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',F));
L=load(In,'MSc'); MS0=L.MSc; clear L
MS=MS0.selectBySrcIndex(P.SrcIdx); clear MS0
fprintf('%s: %d sources, %d epochs, SysRem with %d components, ANNUAL TERM OFF\n', F, MS.Nsrc, MS.Nepoch, NSYS);
T0=tic;
[~,~,~,I2]=ml.scripts.runIterDetrendMSc(MS,'NIterSysRem',NSYS,'AnnualEffect',false,'Verbosity',0);
SysCorX=I2.SysCorX; SysCorY=I2.SysCorY; SrcInd=I2.SrcInd(:); JD=MS.JD(I2.EpochInd); JD=JD(:);
save(Sys,'SysCorX','SysCorY','JD','SrcInd','-v7.3');
fprintf('%s STEP A: %.1f min | correction %dx%d\n', F, toc(T0)/60, size(SysCorX,1), size(SysCorX,2));
clear I2 SysCorX SysCorY
T0=tic; Pre=load(Sys);
[IFsys,~,IFsysB,Info]=ml.scripts.runIterDetrendMSc(MS, ...
    'SysRemCorrection',Pre,'PixPhase',true,'AnnualEffect',false, ...
    'NiterWeightsBeforeSys',15,'NiterWeightsAfterSys',6,'Verbosity',0);
[fx,fy]=IFsys.calculateRstd;
% Identify the target by its MSc index. findClosestSource assumes the fitted
% positions sit near the pixel grid, which is not guaranteed, and it silently
% returns a neighbour when the target is absent or the gauge has moved.
Ie=find(P.SrcIdx(Info.SrcInd(:))==P.TargetMsc);
if isempty(Ie)
    error('v14chain:TargetMissing','the target (MSc %d) is not in the solution', P.TargetMsc);
end
Ichk=IFsys.findClosestSource([150 150]);
fprintf('%s target: index %d at (%.1f,%.1f); findClosestSource would give %d at (%.1f,%.1f)\n', ...
    F, Ie, IFsys.ParS(1,Ie), IFsys.ParS(2,Ie), Ichk, IFsys.ParS(1,Ichk), IFsys.ParS(2,Ichk));
M=IFsys.medianFieldSource({'MAG_PSF'}); M=M(:); Br=M<17&isfinite(M);
Tr=cellfun(@(v) median(v(Br),'omitnan'), IFsys.RMSTrack); n=numel(Tr);
Config=struct('Input',In,'Field',F,'Note','v14: v13 without the annual term', ...
              'Created',string(datetime('now')));
save(Fin,'IFsys','Info','Config','-v7.3');
fprintf('%s STEP B: %.1f min | bright(<17) %.3f/%.3f | target %.2f/%.2f | conv %.4f mas | Nsrc %d\n', ...
    F, toc(T0)/60, median(fx(Br),'omitnan'), median(fy(Br),'omitnan'), fx(Ie), fy(Ie), ...
    max(Tr(n-2:n))-min(Tr(n-2:n)), IFsys.Nsrc);
fprintf('SysRem applied: %d (%s)\nV14 DONE %s\n', Info.SysRemApplied, Info.SysRemSource, F);
