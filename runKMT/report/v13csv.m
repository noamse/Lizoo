addpath('/home/ocs/matlab/Lizoo');
Doc='/home/ocs/matlab/Lizoo/doc/';
RAc=celestial.coo.convertdms('17:52:38.09','gH','d');
Decc=celestial.coo.convertdms('-31:47:36.1','gD','d');
C=load('/home/ocs/KMTdata/Results/v13/cmp.mat'); Sel=C.Sel;
for K=1:2
  Fs={'BLG41','BLG01'}; f=Fs{K};
  V=load(sprintf('/home/ocs/KMTdata/Results/v13/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v13/Tie_%s.mat',f));
  IF=V.IFsys; NF=V.Info; PMabs=T.PMabs; Amap=T.Amap;
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN; Ie=C.Sel.(['Ie_' f]);
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(IF.Data.MAG_PSF,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  Good=~Bad & isfinite(Rx) & isfinite(Ry) & W>0;
  [G,GI]=ml.util.splitEpochGroups(IF.JD);
  List={Ie,''};
  for q=1:numel(Sel.tags)
     if isfinite(Sel.(f)(q)), List(end+1,:)={Sel.(f)(q), Sel.tags{q}}; end
  end
  for L=1:size(List,1)
    Ii=List{L,1}; tag=List{L,2}; g0=Good(:,Ii);
    if sum(g0)<200, continue; end
    eAx=400*rms(Rx(g0,Ii),'omitnan'); eAy=400*rms(Ry(g0,Ii),'omitnan');
    eBx=nan(IF.Nepoch,1); eBy=nan(IF.Nepoch,1); Seas=zeros(IF.Nepoch,1);
    for b=1:GI.Ngroup
      s=(G==b); gg=s&g0; if sum(gg)<50, continue; end
      eBx(s)=400*rms(Rx(gg,Ii),'omitnan'); eBy(s)=400*rms(Ry(gg,Ii),'omitnan'); Seas(s)=b;
    end
    ok=g0&isfinite(eBx);
    JD=IF.JD(ok); t=(JD-IF.JD0)/365.25;
    px=400*(IF.ParS(1,Ii)+IF.ParS(3,Ii).*t)+400*Rx(ok,Ii);
    py=400*(IF.ParS(2,Ii)+IF.ParS(4,Ii).*t)+400*Ry(ok,Ii);
    px=px-mean(px); py=py-mean(py);
    vv=Amap\([IF.ParS(1,Ii);IF.ParS(2,Ii)]-[T.ax(1);T.ay(1)]);
    RAs=RAc+vv(1)/3600/cosd(Decc); Decs=Decc+vv(2)/3600;
    dist=sqrt((IF.ParS(1,Ii)-IF.ParS(1,Ie))^2+(IF.ParS(2,Ii)-IF.ParS(2,Ie))^2);
    if isempty(tag), nm=sprintf('%s',f); ttl='the target';
    else, nm=sprintf('%s_%s',f,tag); ttl=sprintf('calibration star "%s"',tag); end
    Csv=sprintf('%ssource_motion_v13_%s.csv',Doc,nm);
    fid=fopen(Csv,'w');
    fprintf(fid,'# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058, field %s, v13 reduction\n',f);
    fprintf(fid,'# %s\n', ttl);
    fprintf(fid,'# Solution: ~/KMTdata/Results/v13/IFfinal_%s.mat\n',f);
    fprintf(fid,'# v13 fits the calibration set 14<I<19, isolated to 2 arcsec, minus 2-sigma RMS outliers: %d stars, with 6 SysRem components.\n', IF.Nsrc);
    fprintf(fid,'# All are freely fitted; Gaia RUWE<1.4 is applied only in the post-hoc tie.\n');
    fprintf(fid,'# I_OGLE = %.2f, at pixel (%.1f, %.1f), %.1f pix from the target\n', ...
        I(Ii), IF.ParS(1,Ii), IF.ParS(2,Ii), dist);
    fprintf(fid,'# RA, Dec (J2000) = %s , %s  =  %.6f , %.6f deg\n', ...
        celestial.coo.convertdms(RAs,'d','SH'), celestial.coo.convertdms(Decs,'d','SD'), RAs, Decs);
    fprintf(fid,'# absolute proper motion (mu_alpha cos d, mu_d) = %+.3f, %+.3f mas/yr\n', PMabs(1,Ii), PMabs(2,Ii));
    fprintf(fid,'#   full six-parameter gauge removed (constant AND linear gradient)\n');
    fprintf(fid,'# X_mas, Y_mas: position in mas, proper motion NOT removed, relative to the mean\n');
    fprintf(fid,'# errX_decade, errY_decade: this star own residual RMS over the run. Constant.\n');
    fprintf(fid,'# errX_season, errY_season: this star own residual RMS within each season.\n');
    fprintf(fid,'# season: 1=2016 ... 10=2026 (2020 is dropped as a short stub)\n');
    fprintf(fid,'JD,X_mas,Y_mas,errX_decade,errY_decade,errX_season,errY_season,season\n');
    Ex=eBx(ok); Ey=eBy(ok); Sn=Seas(ok);
    for q=1:numel(JD)
      fprintf(fid,'%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n',JD(q),px(q),py(q),eAx,eAy,Ex(q),Ey(q),Sn(q));
    end
    fclose(fid);
    fprintf('  %-22s I %.2f  decade err %6.2f/%6.2f  %d rows\n', nm, I(Ii), eAx, eAy, numel(JD));
  end
  clear IF V
end
fprintf('V6CSV DONE\n');
