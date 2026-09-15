% Seeing-driven blend shift as a PM bias: refit position = p0 + mu*t + k*(fwhm - mean) per epoch
addpath('/home/ocs/matlab/Lizoo');
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120); g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec');
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info'); T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f)); IF=V.IFsys;
  I=V.Info.SrcData.I_ogle(:); I(I>=99)=NaN; Ie=C.Sel.(['Ie_' f]); Ic=C.Sel.(f)(strcmp(C.Sel.tags,'m18_d04'));
  W=IF.calculateWes; [Rx,Ry]=IF.calculateResiduals;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD)|isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0; t=(IF.JD-IF.JD0)/365.25;
  fw=median(IF.Data.fwhm,2,'omitnan'); [G,GI]=ml.util.splitEpochGroups(IF.JD);
  % seeing trend
  fs=arrayfun(@(b) median(fw(G==b),'omitnan'), 1:GI.Ngroup); ts=arrayfun(@(b) mean(t(G==b)), 1:GI.Ngroup);
  A=[ones(GI.Ngroup,1) ts(:)]; c=A\fs(:);
  fprintf('\n%s: median fwhm per season: %s ; trend %+.4f px/yr\n', f, mat2str(round(fs,2)), c(2));
  % joint refit for the target, the neighbour and all stars
  N=IF.Nsrc; dmu=nan(N,2); kk=nan(N,2); sig=nan(N,2);
  for i=1:N
    ok=W(:,i)>0&isfinite(fw); if sum(ok)<500, continue; end
    px=400*(IF.ParS(3,i)*t(ok))+400*Rx(ok,i); py=400*(IF.ParS(4,i)*t(ok))+400*Ry(ok,i);   % PM kept, mas
    A=[ones(sum(ok),1) t(ok) fw(ok)-median(fw(ok))];
    cx=A\px; cy=A\py;
    dmu(i,:)=[cx(2)-400*IF.ParS(3,i), cy(2)-400*IF.ParS(4,i)]; kk(i,:)=[cx(3) cy(3)];
    Ainv=inv(A'*A); sig(i,:)=[std(px-A*cx) std(py-A*cy)]*sqrt(Ainv(2,2));
  end
  % to the sky
  S=@(v) (1000*(T.Amap\(v'/1000)))';
  for k=1:2
    if k==1, i=Ie; nm='target'; else, i=Ic; nm='neighbour'; end
    ds=S(dmu(i,:)); ks=S(kk(i,:)); mu0=T.PMabs(:,i)'; mu1=mu0+ds;
    fprintf('   %-9s k = %+6.1f / %+6.1f mas per px of fwhm (sky) | PM change %+.2f / %+.2f mas/yr | PM %+.2f / %+.2f -> %+.2f / %+.2f | Gaia %+.2f / %+.2f\n', ...
       nm, ks, ds, mu0, mu1, PMRA(T.J(i)), PMDec(T.J(i)));
  end
  gg=isfinite(dmu(:,1))&(1:N)'~=Ie&(1:N)'~=Ic;
  fprintf('   field stars: |k| median %.1f / %.1f mas/px ; PM change rstd %.3f / %.3f mas/yr (pixel axes), 18<I<19: %.3f / %.3f\n', ...
     median(abs(kk(gg,1))), median(abs(kk(gg,2))), tools.math.stat.rstd(dmu(gg,1)), tools.math.stat.rstd(dmu(gg,2)), tools.math.stat.rstd(dmu(gg&I>18,1)), tools.math.stat.rstd(dmu(gg&I>18,2)));
end
