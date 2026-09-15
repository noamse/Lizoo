% Do single-detection epochs carry the blend centroid? and are there duplicate assignments?
addpath('/home/ocs/matlab/Lizoo'); C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys'); IF=V.IFsys;
  Ie=C.Sel.(['Ie_' f]); Ic=C.Sel.(f)(strcmp(C.Sel.tags,'m18_d04'));
  W=IF.calculateWes; [Rx,Ry]=IF.calculateResiduals; ok=isfinite(Rx)&isfinite(Ry)&W>0;
  X=IF.Data.X; Y=IF.Data.Y; detT=isfinite(X(:,Ie)); detN=isfinite(X(:,Ic));
  dp=IF.ParS(1:2,Ic)-IF.ParS(1:2,Ie); u=dp/norm(dp);
  same=detT&detN&abs(X(:,Ie)-X(:,Ic))<1e-6&abs(Y(:,Ie)-Y(:,Ic))<1e-6;
  fprintf('\n%s: raw detections target %d, neighbour %d, both %d, identical X,Y in both: %d\n', f, sum(detT), sum(detN), sum(detT&detN), sum(same));
  cls={detT&detN,'both detected'; detT&~detN,'target only'};
  for k=1:2
    s=cls{k,1}&ok(:,Ie); r=400*[Rx(s,Ie) Ry(s,Ie)]*u;
    fprintf('   target  in %-14s (n=%5d): mean residual along pair axis %+7.1f mas (s.e. %.1f), rms %.0f\n', cls{k,2}, sum(s), mean(r), std(r)/sqrt(numel(r)), std(r));
  end
  cls={detT&detN,'both detected'; ~detT&detN,'neighbour only'};
  for k=1:2
    s=cls{k,1}&ok(:,Ic); r=400*[Rx(s,Ic) Ry(s,Ic)]*u;
    fprintf('   neighb. in %-14s (n=%5d): mean residual along pair axis %+7.1f mas (s.e. %.1f), rms %.0f\n', cls{k,2}, sum(s), mean(r), std(r)/sqrt(numel(r)), std(r));
  end
  % magnitudes per class
  M=IF.Data.MAG_PSF;
  fprintf('   target KMT mag: both-detected %.2f, target-only %.2f | neighbour: both %.2f, neighbour-only %.2f\n', median(M(detT&detN,Ie),'omitnan'), median(M(detT&~detN,Ie),'omitnan'), median(M(detT&detN,Ic),'omitnan'), median(M(~detT&detN,Ic),'omitnan'));
  % seeing per class
  fw=median(IF.Data.fwhm,2,'omitnan');
  fprintf('   median fwhm: both %.2f, target-only %.2f, neighbour-only %.2f, neither %.2f px\n', median(fw(detT&detN),'omitnan'), median(fw(detT&~detN),'omitnan'), median(fw(~detT&detN),'omitnan'), median(fw(~detT&~detN),'omitnan'));
  % relative PM of the pair from both-detected epochs only, raw positions, straight line fit in time
  t=(IF.JD-IF.JD0)/365.25; s=detT&detN&ok(:,Ie)&ok(:,Ic);
  dxr=400*(X(s,Ic)-X(s,Ie)); dyr=400*(Y(s,Ic)-Y(s,Ie)); A=[ones(sum(s),1) t(s)];
  cx=A\dxr; cy=A\dyr; ex=dxr-A*cx; ey=dyr-A*cy; Sxx=sum((t(s)-mean(t(s))).^2);
  fprintf('   raw pair separation drift (both-detected, n=%d): %+.2f+-%.2f / %+.2f+-%.2f mas/yr (pixel X,Y; X opposite to RA) -> sky %+.2f / %+.2f ; Gaia nb-target +4.22 / +9.40\n', sum(s), cx(2), std(ex)/sqrt(Sxx), cy(2), std(ey)/sqrt(Sxx), -cx(2), cy(2));
end
