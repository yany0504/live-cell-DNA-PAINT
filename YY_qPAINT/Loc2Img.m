function PALMLoc = Loc2Img(Orte,Mag,height,width)

%Pixel=1/Pixel;
PALMLoc=uint16(zeros(floor(Mag.*height),floor(Mag.*width)));
%PALMLoc=uint32(zeros(floor(Mag.*height),floor(Mag.*width)));
% erstellt Bildmatrix mit der gewünschten effektive Pixelgröße
for o=1:size(Orte,1)
PALMLoc(ceil(Orte(o,2).*Mag),ceil(Orte(o,1).*Mag))=PALMLoc(ceil(Orte(o,2).*Mag),ceil(Orte(o,1).*Mag))+1;
end
% alle Elemente durchgehen und den Ortsvektor im Bild
% inkrementieren
%Mean=mean(mean(Orte(:,7:8)'))
end
