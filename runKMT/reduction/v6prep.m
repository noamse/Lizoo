addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
OGf=[Doc(1:end-4) 'OGLEdata/OB260058/OB160058.mat'];
for F={'BLG41','BLG01'}
  f=F{1};
  S=load(sprintf('/home/ocs/KMTdata/Results/v3_Final/IFfinal_260058_CTIO_%s.mat',f));
  IF=S.IFsys; NF=S.Info;
  Iog=NF.SrcData.I_ogle(:); Iog(Iog>=99)=NaN; Ie=IF.findClosestSource([150 150]);
  Cat=ml.util.ogleCompanionCat(OGf,'Field',f);
  Cal=logical(ml.util.selectRefSources(IF,'RefMag',Iog,'MagRange',[14 19], ...
      'CompanionRadius',5,'CompanionMaxMag',18,'CompanionCat',Cat)); Cal=Cal(:);
  Cal(Ie)=true;                       % belt and braces; it already passes
  SrcIdx=NF.SrcInd(Cal);
  save(sprintf('/home/ocs/KMTdata/GaiaRef/prep6_%s.mat',f),'SrcIdx','-v7.3');
  fprintf('%s: %d sources selected (target included)\n', f, numel(SrcIdx));
  clear IF S
end
