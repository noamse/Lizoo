addpath('~/matlab/Lizoo');
P=load('/home/ocs/KMTdata/GaiaRef/prep8_BLG41.mat');
L=load('/home/ocs/KMTdata/Results16_26_v3/KMT_260058_BLG41_MSc.mat','MSc');
MS=L.MSc.selectBySrcIndex(P.SrcIdx); clear L
MS=MS.selectByEpoch(1:500);
for nb=[6 12 20]
  [IF,~,~,~]=ml.scripts.runIterDetrendMSc(MS,'NColourBins',nb,'Verbosity',0, ...
      'NiterNoWeightsBeforeSys',1,'NiterWeightsBeforeSys',2,'NiterWeightsAfterSys',1);
  [NC,ed,bc]=IF.generateBins;
  fprintf('requested %2d bins -> IF.NColourBins %2d, bins used %2d, stars per bin %s\n', ...
      nb, IF.NColourBins, numel(NC), mat2str(NC(:)'));
  fprintf('    ParHalat %s (18 params x bins)\n', mat2str(size(IF.ParHalat)));
  clear IF
end
