% Control: for other close pairs in the calibration set, does the raw separation drift agree with Gaia's relative PM?
addpath('/home/ocs/matlab/Lizoo');
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120); g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec'); Gm=g('phot_g_mean_mag');
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info'); T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f)); IF=V.IFsys;
  I=V.Info.SrcData.I_ogle(:); I(I>=99)=NaN; Ie=C.Sel.(['Ie_' f]); Ic=C.Sel.(f)(strcmp(C.Sel.tags,'m18_d04'));
  X=IF.Data.X; Y=IF.Data.Y; W=IF.calculateWes; [Rx,Ry]=IF.calculateResiduals; ok=isfinite(Rx)&isfinite(Ry)&W>0; t=(IF.JD-IF.JD0)/365.25;
  P=IF.ParS(1:2,:)'; N=IF.Nsrc; J=T.J; gm=T.Tie.Matched>0&isfinite(J);
  D=squareform(pdist(P)); [a,b]=find(triu(D>0&D<20,1));   % pairs closer than 8 arcsec
  fprintf('\n%s: %d pairs within 8 arcsec among the %d stars\n', f, numel(a), N);
  fprintf('   %-5s %-5s %5s %5s %5s | %14s | %14s | %14s\n','i','j','sep"','Ii','Ij','ours dRA/dDec','Gaia dRA/dDec','diff (sigma)');
  R=[];
  for k=1:numel(a)
    i=a(k); j=b(k); if ~(gm(i)&&gm(j)), continue; end
    s=ok(:,i)&ok(:,j); if sum(s)<500, continue; end
    dx=400*(X(s,j)-X(s,i)); dy=400*(Y(s,j)-Y(s,i)); A=[ones(sum(s),1) t(s)]; cx=A\dx; cy=A\dy;
    Sxx=sum((t(s)-mean(t(s))).^2); ex=std(dx-A*cx)/sqrt(Sxx); ey=std(dy-A*cy)/sqrt(Sxx);
    sky=T.Amap\[cx(2);cy(2)]*1;  % pixel/yr -> arcsec/yr? no: cx is mas/yr in pixel axes; Amap maps arcsec->pix (scale ~2.52), so sky mas/yr = Amap\ (mas/yr pix)
    ours=1000*(T.Amap\([cx(2);cy(2)]/1000/0.4*0.4));   % keep units: cx in mas(pixel-scaled)/yr -> pix/yr = cx/400 ; sky arcsec/yr = Amap\pix ; mas/yr = *1000
    ours=1000*(T.Amap\([cx(2);cy(2)]/400));
    err=1000*abs(T.Amap\([ex;ey]/400));
    gaia=[PMRA(J(j))-PMRA(J(i)); PMDec(J(j))-PMDec(J(i))];
    R(end+1,:)=[D(i,j)*0.4 ours' gaia' err' I(i) I(j) sum(s)];
    tag=''; if (i==Ie&&j==Ic)||(i==Ic&&j==Ie), tag='  <-- target pair'; end
    fprintf('   %-5d %-5d %5.2f %5.1f %5.1f | %+6.2f %+6.2f  | %+6.2f %+6.2f  | %+5.1f %+5.1f%s\n', i,j,D(i,j)*0.4,I(i),I(j),ours,gaia,(ours-gaia)./err',tag);
  end
  d=R(:,2:3)-R(:,4:5); e=R(:,6:7);
  fprintf('   %d pairs: rstd(ours-Gaia) %.2f / %.2f mas/yr ; median |diff/err| %.1f / %.1f ; sep<3": n=%d rstd %.2f / %.2f\n', size(R,1), tools.math.stat.rstd(d(:,1)), tools.math.stat.rstd(d(:,2)), median(abs(d(:,1)./e(:,1))), median(abs(d(:,2)./e(:,2))), sum(R(:,1)<3), tools.math.stat.rstd(d(R(:,1)<3,1)), tools.math.stat.rstd(d(R(:,1)<3,2)));
end
