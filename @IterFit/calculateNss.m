function Nss = calculateNss(IF)


Bss = zeros([numel(IF.ParS(:,1)),numel(IF.ParS(:,1)),IF.Nsrc]);

%[Ax,Ay] = IF.generateSourceDesignMat;
Ax = IF.AsX;
Ay = IF.AsY;

Wes = calculateWes(IF);

% A parameter that is held rather than fitted is taken out of its source's
% normal matrix: the row and the column are cleared and the diagonal put back at
% the scale the block already carries. calculateBs clears the matching right
% hand side, so the increment of a held parameter comes out exactly zero and the
% free ones solve the constrained system, rather than being taken from the
% unconstrained one and overwritten afterwards. The diagonal is matched to the
% block instead of being set to one because the weights are normalised to sum to
% unity, which leaves the real entries orders of magnitude below one and would
% stop bicg from converging.
Held = [];
if ~isempty(IF.ParSFixed)
    Held = isfinite(IF.ParSFixed);
end

for Isrc = 1:IF.Nsrc
     W   = Wes(:,Isrc);
     Blk = Bss(:,:,Isrc) + (Ax'*(Ax.*W) + Ay'*(Ay.*W));
     if ~isempty(Held) && any(Held(:,Isrc))
         Ih  = find(Held(:,Isrc));
         Dg  = diag(Blk);
         Rep = Dg(Ih);
         Rep(~isfinite(Rep) | Rep<=0) = 1;    % a source carrying no weight at all
         Blk(Ih,:) = 0;
         Blk(:,Ih) = 0;
         Blk(sub2ind(size(Blk),Ih,Ih)) = Rep;
     end
     Bss(:,:,Isrc)= Blk;

end

Nss = sparse(Bss(:,:,1));

for Iblk = 2:numel(Bss(1,1,:)); Nss = blkdiag(Nss,Bss(:,:,Iblk));end

end