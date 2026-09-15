% Is the target/neighbour pair statically blended? positions, relative PM, seeing dependence, detection overlap
addpath('/home/ocs/matlab/Lizoo');
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120);
g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec'); Gm=g('phot_g_mean_mag');
Xi=(g('RA')*180/pi-RAd).*cosd(Decd).*3600; Eta=(g('Dec')*180/pi-Decd).*3600;
d=sqrt(Xi.^2+Eta.^2); [~,o]=sort(d); it=o(1); in=o(3); i3=o(2);   % target, the G=19.69 neighbour, the G=20.33 one
fprintf('Gaia: target G %.2f | nb G %.2f at (%+.2f,%+.2f) | faint G %.2f at (%+.2f,%+.2f); nb-faint sep %.2f arcsec\n', Gm(it),Gm(in),Xi(in)-Xi(it),Eta(in)-Eta(it),Gm(i3),Xi(i3)-Xi(it),Eta(i3)-Eta(it), hypot(Xi(in)-Xi(i3),Eta(in)-Eta(i3)));
fprintf('Gaia quality cols: %s\n', strjoin(GCol(contains(GCol,'astrometric')|contains(GCol,'non_single')),', '));
for nm={'astrometric_excess_noise','astrometric_excess_noise_sig','astrometric_gof_al','astrometric_n_good_obs_al','astrometric_params_solved','non_single_star'}
  c=find(strcmp(GCol,nm{1})); if ~isempty(c), fprintf('   %-32s target %8.3f  nb %8.3f  faint %8.3f\n', nm{1}, GC(it,c), GC(in,c), GC(i3,c)); end
end
dG=[Xi(in)-Xi(it); Eta(in)-Eta(it)];  % Gaia separation vector, arcsec, target->nb
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info'); T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f));
  IF=V.IFsys; Ie=C.Sel.(['Ie_' f]); Ic=C.Sel.(f)(strcmp(C.Sel.tags,'m18_d04'));
  % measured separation on the sky (Amap: arcsec -> pix)
  dp=IF.ParS(1:2,Ic)-IF.ParS(1:2,Ie); ds=T.Amap\dp;
  fprintf('\n%s: measured nb-target separation %.2f pix = (%+.2f,%+.2f) arcsec, |%.2f| ; Gaia (%+.2f,%+.2f) |%.2f| -> ratio %.3f\n', f, norm(dp), ds, norm(ds), dG, norm(dG), norm(ds)/norm(dG));
  % relative PM
  rel=T.PMabs(:,Ic)-T.PMabs(:,Ie); relG=[PMRA(in)-PMRA(it); PMDec(in)-PMDec(it)];
  fprintf('   relative PM nb-target: ours %+.2f / %+.2f ; Gaia %+.2f / %+.2f  (ratio along Gaia direction %.2f)\n', rel, relG, dot(rel,relG)/dot(relG,relG));
  % detection overlap
  W=IF.calculateWes; [Rx,Ry]=IF.calculateResiduals; ok=isfinite(Rx)&isfinite(Ry)&W>0;
  fprintf('   epochs with target %d, neighbour %d, both %d, target only %d, neighbour only %d (of %d)\n', sum(ok(:,Ie)), sum(ok(:,Ic)), sum(ok(:,Ie)&ok(:,Ic)), sum(ok(:,Ie)&~ok(:,Ic)), sum(~ok(:,Ie)&ok(:,Ic)), IF.Nepoch);
  % seeing dependence of the residual along the separation direction
  fw=median(IF.Data.fwhm,2,'omitnan'); u=dp/norm(dp);
  for k=1:2
    if k==1, i=Ie; nm='target'; else, i=Ic; nm='neighbour'; end
    s=ok(:,i)&isfinite(fw); r=400*[Rx(s,i) Ry(s,i)]*u; A=[ones(sum(s),1) fw(s)-median(fw(s))]; c=A\r;
    fprintf('   %s: residual along the pair axis vs fwhm: slope %+.2f mas per pixel of fwhm (fwhm range %.1f-%.1f px, 10-90%%: %.1f-%.1f)\n', nm, c(2), min(fw(s)), max(fw(s)), quantile(fw(s),.1), quantile(fw(s),.9));
    % per-epoch measured separation when both detected
  end
  s=ok(:,Ie)&ok(:,Ic)&isfinite(fw); sep=400*(([Rx(s,Ic) Ry(s,Ic)]-[Rx(s,Ie) Ry(s,Ie)])*u); A=[ones(sum(s),1) fw(s)-median(fw(s))]; c=A\sep;
  fprintf('   measured separation (both detected, n=%d) vs fwhm: slope %+.2f mas/px ; separation rms %.1f mas\n', sum(s), c(2), std(sep));
  % relative position by season (does the separation drift as Gaia predicts, 9.4 mas/yr?)
  [G,GI]=ml.util.splitEpochGroups(IF.JD); t=(IF.JD-IF.JD0)/365.25; tb=[]; sb=[];
  for b=1:GI.Ngroup, q=s&G==b; if sum(q)<30, continue; end; tb(end+1)=mean(t(q)); v=400*(([Rx(q,Ic) Ry(q,Ic)]-[Rx(q,Ie) Ry(q,Ie)])); sb(end+1,:)=mean(v,1); end
  A=[ones(numel(tb),1) tb(:)]; cx=A\sb(:,1); cy=A\sb(:,2);
  fprintf('   separation drift from residuals of both-detected epochs, per season: %+.2f / %+.2f mas/yr (pixel axes), %d seasons\n', cx(2), cy(2), numel(tb));
end
