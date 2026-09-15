% Re-match the saved MSc files to OGLE with the corrected per-field offset.
% Only SrcData.I_ogle and SrcData.V_ogle depend on it: FlagGoodSrc is computed
% before the match and the source and epoch lists are untouched, so this is
% exactly what rerunning KMT_pipelineI would produce.
addpath('~/matlab/Lizoo');
Off = struct('BLG41',[229.202 229.283], 'BLG01',[227.542 229.689]);
L0=load('/home/ocs/matlab/Lizoo/OGLEdata/OB260058/OB160058.mat');
V=fieldnames(L0); OB=L0.(V{1});
for F=["BLG41","BLG01"]
  In =sprintf('/home/ocs/KMTdata/Results16_26_v2/KMT_260058_%s_MSc.mat',F);
  Out=sprintf('/home/ocs/KMTdata/Results16_26_v3/KMT_260058_%s_MSc.mat',F);
  S=load(In); MSc=S.MSc; JD=S.JD; clear S
  Old=MSc.SrcData.I_ogle(:);
  T=OB.CatData.Table;
  D=Off.(F);
  Xo=(double(T.corrX)-D(1)).*0.26./0.4+150;
  Yo=(double(T.corrY)-D(2)).*0.26./0.4+150;
  Io=double(T.I); Vo=double(T.V);
  k=Io<19 & isfinite(Xo) & isfinite(Yo);            % the pipeline's own cut
  Xo=Xo(k); Yo=Yo(k); Io=Io(k); Vo=Vo(k);
  [Yo,so]=sort(Yo); Xo=Xo(so); Io=Io(so); Vo=Vo(so); % matchCatalogsXY needs Y-sorted
  [Ind1,~,~,~,~,~] = imUtil.match.mex.matchCatalogsXY(MSc.SrcData.X(:), MSc.SrcData.Y(:), Xo, Yo, 1.5);
  NN=~isnan(Ind1);
  MSc.SrcData.I_ogle = nan(1,MSc.Nsrc); MSc.SrcData.V_ogle = nan(1,MSc.Nsrc);
  MSc.SrcData.I_ogle(NN)=Io(Ind1(NN));  MSc.SrcData.V_ogle(NN)=Vo(Ind1(NN));
  New=MSc.SrcData.I_ogle(:);
  fprintf('%s: I_ogle %d -> %d of %d (%.1f%% -> %.1f%%) | V_ogle %d\n', F, ...
      sum(isfinite(Old)), sum(isfinite(New)), MSc.Nsrc, ...
      100*mean(isfinite(Old)), 100*mean(isfinite(New)), ...
      sum(isfinite(MSc.SrcData.V_ogle(:))));
  both=isfinite(Old)&isfinite(New);
  fprintf('   of the %d matched both ways, %d changed star (|dI|>0.05): median |dI| %.3f\n', ...
      sum(both), sum(abs(Old(both)-New(both))>0.05), median(abs(Old(both)-New(both))));
  save(Out,'MSc','JD','-v7.3');
  fprintf('   wrote %s\n', Out);
  clear MSc
end
fprintf('REMATCH DONE\n');
function v=S0(~), v=[]; end
