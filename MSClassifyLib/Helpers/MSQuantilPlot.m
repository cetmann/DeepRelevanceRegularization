function [ ] = MSQuantilPlot( X,Y,color )
%MSQuantilPlot plotet die Minima, untere Quartil, Median, obere
%Quartil und Maxima der Daten.
% X sind die Daten für die xAchse
% Y hat die Dimension 5xN wobei N die Länge von X ist.
%   In der ersten Zeile von Y stehen die Maxima.
%   In der zweiten Zeile von Y stehen die oberen Quartil.
%   In der dritten Zeile von Y stehen die Medians.
%   In der vierten Zeile von Y stehen die unteren Quartil.
%   In der fünften Zeile von Y stehen die Minima.
% color ist ein string und gibt optional die Farbe an.

% Die Daten können einfach mit der Funktion quantile erzeugt werden:
% Y = quantile(Data,[1,0.75,0.5,0.25,0]);
% Die Funktion quantile erhebt spaltenweise Minimum, unteres Quartil, Median, oberes
% Quartil und Maximum der Daten.

if nargin<3
    color = 'r';
end
   
cplot = fill([X,flip(X)],[Y(4,:),flip(Y(2,:))],color,'LineStyle','none');
set(cplot, 'FaceAlpha', 0.15);
hold on
plot(X,Y(3,:),color);
plot(X,Y(5,:),[color '--']);
plot(X,Y(1,:),[color '--']);
hold off

end

