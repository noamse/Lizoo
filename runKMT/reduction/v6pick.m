addpath('/home/ocs/matlab/Lizoo');
Fs={'BLG41','BLG01'}; D=struct();
for k=1:2
  f=Fs{k};
  V=load(sprintf('/home/ocs/KMTdata/Results/v6/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v6/Tie_%s.mat',f));
  IF=V.IFsys; NF=V.Info;
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  Ie=IF.findClosestSource([150 150]);
  d=sqrt((IF.ParS(1,:)'-IF.ParS(1,Ie)).^2+(IF.ParS(2,:)'-IF.ParS(2,Ie)).^2);
  D.(f)=struct('I',I,'d',d,'Ie',Ie,'gid',T.J,'ok',T.Tie.Matched>0&isfinite(T.J), ...
               'X',IF.ParS(1,:)','Y',IF.ParS(2,:)');
  clear IF V
end
% stars present in BOTH fields, by Gaia id
g1=D.BLG41.gid; g2=D.BLG01.gid;
pick=@(lo,hi,n) deal([]);
fprintf('%-10s %-6s %-7s %-22s %-22s\n','tag','I','dist','BLG41 idx / pixel','BLG01 idx / pixel');
Sel=struct('BLG41',[],'BLG01',[],'tags',{{}});
for band={[17.8 18.4],'f'; [16.0 17.0],'b'}'
end
Bands={[17.8 18.4],'m18'; [16.0 17.0],'m16'};
for q=1:2
  lo=Bands{q,1}(1); hi=Bands{q,1}(2); pre=Bands{q,2};
  c1=find(D.BLG41.ok & D.BLG41.I>=lo & D.BLG41.I<hi & (1:numel(D.BLG41.I))'~=D.BLG41.Ie);
  c2=find(D.BLG01.ok & D.BLG01.I>=lo & D.BLG01.I<hi & (1:numel(D.BLG01.I))'~=D.BLG01.Ie);
  [~,i1,i2]=intersect(g1(c1),g2(c2));
  a=c1(i1); b=c2(i2);
  dd=(D.BLG41.d(a)+D.BLG01.d(b))/2;
  [~,o]=sort(dd); o=o(1:min(3,numel(o)));
  for z=1:numel(o)
    tag=sprintf('%s_d%02.0f',pre,D.BLG41.d(a(o(z))));
    Sel.BLG41(end+1)=a(o(z)); Sel.BLG01(end+1)=b(o(z)); Sel.tags{end+1}=tag;
    fprintf('%-10s %6.2f %6.1f   %4d (%5.1f,%5.1f)   %4d (%5.1f,%5.1f)\n', tag, ...
       D.BLG41.I(a(o(z))), D.BLG41.d(a(o(z))), a(o(z)), D.BLG41.X(a(o(z))), D.BLG41.Y(a(o(z))), ...
       b(o(z)), D.BLG01.X(b(o(z))), D.BLG01.Y(b(o(z))));
  end
end
save('/home/ocs/KMTdata/Results/v6/cmp.mat','Sel','-v7.3');
fprintf('\nsaved %d comparison stars\n', numel(Sel.tags));
