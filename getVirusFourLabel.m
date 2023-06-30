function [oneVirus] = getVirusFourLabel(oneVirusSeq)
%GETVIRUSFOURLABEL  
% Charge, pH, Mass, Hydro indice

oneVirusSeqE = oneVirusSeq(:,1);
oneVirusSeqPH = oneVirusSeq(:,2);
oneVirusSeqM = oneVirusSeq(:,3);
oneVirusSeqHydro = oneVirusSeq(:,4);

oneVirusSeqE(isnan(oneVirusSeqE)) = [];
oneVirusSeqPH(isnan(oneVirusSeqPH)) = [];
oneVirusSeqM(isnan(oneVirusSeqM)) = [];
oneVirusSeqHydro(isnan(oneVirusSeqHydro)) = [];

oneVirus(:,1)=oneVirusSeqE;
oneVirus(:,2)=oneVirusSeqPH;
oneVirus(:,3)=oneVirusSeqM;
oneVirus(:,4)=oneVirusSeqHydro;

end

