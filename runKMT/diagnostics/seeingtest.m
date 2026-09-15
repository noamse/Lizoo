% Hypothesis: each star has its own seeing-dependent centroid bias (set by its
% blend geometry). Seeing drifts within a season, so each star gets its own
% within-season slope, reproducing between fields because the geometry is fixed.
%   predicted slope(i,b) = sensitivity_i * dFWHM_b
addpath('/home/ocs/matlab/Lizoo');
for F={'BLG41','BLG01'}
  f=F{1};
  V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f));
  IF=V.IFsys; NF=V.Info; I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  fw=median(IF.Data.fwhm,2,'omitnan');   % per-epoch seeing, median over sources
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup;
  fprintf('\n===== %s =====  seeing %.2f..%.2f pix, median %.2f\n', f, ...
      min(fw,[],'omitnan'), max(fw,[],'omitnan'), median(fw,'omitnan'));
  % within-season drift of seeing
  dfw=nan(nb,1);
  for b=1:nb
    s=G==b & isfinite(fw); if sum(s)<50, continue; end
    t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t));
    c=[ones(sum(s),1) t]\fw(s); dfw(b)=c(2);
  end
  fprintf('  within-season seeing drift by season: %s\n', mat2str(round(dfw',3)));
  % per-star seeing sensitivity, and per-star-season slopes
  sens=nan(IF.Nsrc,2); SL=nan(IF.Nsrc,nb,2);
  for i=1:IF.Nsrc
    ok=W(:,i)>0 & isfinite(fw);
    if sum(ok)<300, continue; end
    q=fw(ok)-mean(fw(ok));
    A=[ones(sum(ok),1) q];
    cx=A\(400*Rx(ok,i)); cy=A\(400*Ry(ok,i));
    sens(i,:)=[cx(2) cy(2)];               % mas per pix of seeing
    for b=1:nb
      s=W(:,i)>0 & G==b; if sum(s)<50, continue; end
      t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t)); B=[ones(sum(s),1) t];
      a1=B\(400*Rx(s,i)); a2=B\(400*Ry(s,i));
      SL(i,b,1)=a1(2); SL(i,b,2)=a2(2);
    end
  end
  fprintf('  per-star seeing sensitivity: median |dX/dFWHM| %.2f mas/pix, spread %.2f\n', ...
      median(abs(sens(:,1)),'omitnan'), tools.math.stat.rstd(sens(:,1)));
  % does sensitivity x seeing drift predict the slopes?
  for ax=1:2
    pred=[]; obs=[];
    for i=1:IF.Nsrc
      if ~isfinite(sens(i,ax)), continue; end
      for b=1:nb
        if ~isfinite(SL(i,b,ax)) || ~isfinite(dfw(b)), continue; end
        pred(end+1,1)=sens(i,ax)*dfw(b); obs(end+1,1)=SL(i,b,ax);
      end
    end
    c=corr(pred,obs);
    r=obs-pred;
    fprintf('  axis %d: corr(predicted, observed slope) = %+.3f over %d star-seasons\n', ax, c, numel(obs));
    fprintf('           slope rms  observed %.2f -> after removing the prediction %.2f mas (%.0f%% of variance)\n', ...
        std(obs), std(r), 100*(1-var(r)/var(obs)));
  end
  clear IF V
end
