% v8 vs v11/v11z on the SAME epochs: is the RMS gain a cleaner solution or just cleaner epochs?
addpath('/home/ocs/matlab/Lizoo'); C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for R={'v11','v11z'}
  r=R{1};
  for F={'BLG41','BLG01'}
    f=F{1}; V8=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info'); Vn=load(sprintf('/home/ocs/KMTdata/Results/%s/IFfinal_%s.mat',r,f),'IFsys','Info');
    P8=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',f));
    m8=P8.SrcIdx(V8.Info.SrcInd(:)); mn=P8.SrcIdx(Vn.Info.SrcInd(:)); [~,i8,in]=intersect(m8,mn);
    I=V8.Info.SrcData.I_ogle(:); I(I>=99)=NaN; I=I(i8); Ie8=C.Sel.(['Ie_' f]); ie=find(i8==Ie8);
    [Rx8,Ry8]=V8.IFsys.calculateResiduals; [Rxn,Ryn]=Vn.IFsys.calculateResiduals;
    W8=V8.IFsys.calculateWes; Wn=Vn.IFsys.calculateWes;
    [tf,loc]=ismember(V8.IFsys.JD, Vn.IFsys.JD);   % v8 epochs that survive in the variant
    fx8=nan(numel(i8),1); fy8=fx8; fxn=fx8; fyn=fx8;
    for k=1:numel(i8)
      a=i8(k); b=in(k);
      s=tf&W8(:,a)>0; fx8(k)=400*tools.math.stat.rstd(Rx8(s,a)); fy8(k)=400*tools.math.stat.rstd(Ry8(s,a));
      s=Wn(:,b)>0;     fxn(k)=400*tools.math.stat.rstd(Rxn(s,b)); fyn(k)=400*tools.math.stat.rstd(Ryn(s,b));
    end
    fprintf('%s %s: %d common epochs of %d (v8) ; %d stars\n', r, f, sum(tf), numel(tf), numel(i8));
    for e=[14 17;17 18;18 19]'
      s=I>=e(1)&I<e(2);
      fprintf('   I %2d-%2d (n=%3d): v8 on same epochs %.3f/%.3f -> %s %.3f/%.3f  (%+.1f%% / %+.1f%%)\n', e, sum(s), median(fx8(s)),median(fy8(s)),r,median(fxn(s)),median(fyn(s)), 100*(median(fxn(s))/median(fx8(s))-1), 100*(median(fyn(s))/median(fy8(s))-1));
    end
    fprintf('   target: v8 on same epochs %.2f/%.2f -> %s %.2f/%.2f ; v8 all epochs %.2f/%.2f\n', fx8(ie),fy8(ie),r,fxn(ie),fyn(ie), 400*tools.math.stat.rstd(Rx8(W8(:,Ie8)>0,Ie8)), 400*tools.math.stat.rstd(Ry8(W8(:,Ie8)>0,Ie8)));
    s=I<17; fprintf('   stars improved (I<17): X %.0f%%, Y %.0f%%\n', 100*mean(fxn(s)<fx8(s)), 100*mean(fyn(s)<fy8(s)));
  end
end
