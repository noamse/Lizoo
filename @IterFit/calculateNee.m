function Nee = calculateNee(IF,Args)

arguments
    IF;
    Args.Chromatic = false;
end


Bee = zeros([numel(IF.ParE(:,1)),numel(IF.ParE(:,1)),IF.Nepoch]);

W = calculateWes(IF);

%[Aex,Aey] = IF.generateEpochDesignMat;
Aex = IF.AeX;
Aey = IF.AeY;
W = median(W,1,'omitnan');
% Restrict the frame to the reference sources, when one has been chosen. The
% weight here is per source, so zeroing it keeps a source out of the epoch
% solution while leaving its own parameters to be solved as usual.
if ~isempty(IF.RefSrcFlag)
    W(~IF.RefSrcFlag) = 0;
end
% if IF.Chromatic
%     pa = IF.getTimeSeriesField(1,{'pa'});
%     
%     for Iep=1:IF.Nepoch
%         AexC=Aex;
%         AeyC=Aey;
%         
%         if isnan(pa(Iep))
%             AexC(:,7) = Aex(:,7).*0;
%             AeyC(:,8) = Aey(:,8).*0;    
%         else
%             %AexC(:,7) = -Aex(:,7).*sin(pa(Iep)) + Aey(:,8).*cos(pa(Iep));
%             %AexC(:,8) = zeros(size(AexC(:,8)));
%             AexC(:,7) = Aex(:,7).*sin(pa(Iep)) ;
%             AeyC(:,8) = Aey(:,8).*cos(pa(Iep));
%         end
%         Bee(:,:,Iep)= Bee(:,:,Iep) + (AexC'.*W)*AexC + (AeyC'.*W)*AeyC;
%     end
% else
    for Iep=1:IF.Nepoch
        Bee(:,:,Iep)= Bee(:,:,Iep) + (Aex'.*W)*Aex + (Aey'.*W)*Aey;
    end
% end


Nee = sparse(Bee(:,:,1));

for Iblk = 2:numel(Bee(1,1,:)); Nee = blkdiag(Nee,Bee(:,:,Iblk));end

end
