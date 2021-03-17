clear

v = load('fcch_rad.mat');
ring_rad = v.ring;

%Get 6D orbit no bb / w. radiations
co = findorbit6(ring_rad,1:length(ring_rad)+1);
cox_i = co(1,:);
cop_i = co(5,:);
spos_i = findspos(ring_rad,1:length(ring_rad)+1);

bends = findcells(ring_rad,'Class','Bend');
quads = findcells(ring_rad,'Class','Quadrupole');
sexts = findcells(ring_rad,'Class','Sextupole');

k0 = atgetfieldvalues(ring_rad(bends),'BendingAngle')./atgetfieldvalues(ring_rad(bends),'Length');
%k1 = atgetfieldvalues(ring_rad(quads),'PolynomB',{2});
%k2 = atgetfieldvalues(ring_rad(sexts),'PolynomB',{3});

ring_rad(bends) = atsetfieldvalues(ring_rad(bends),'PolynomB',{1},k0.*cop_i(bends)');
%ring_rad(quads) = atsetfieldvalues(ring_rad(quads),'PolynomB',{2},k1.*(1+cop_i(quads)'));
%ring_rad(sexts) = atsetfieldvalues(ring_rad(sexts),'PolynomB',{3},k2.*(1+cop_i(sexts)'));

%Get 6D orbit no bb / w. radiations
co_t = findorbit6(ring_rad,1:length(ring_rad)+1);
cox_t = co_t(1,:);
cop_t = co_t(5,:);
spos_t = findspos(ring_rad,1:length(ring_rad)+1);

figure(1)
plot(spos_i,cox_i,spos_t,cox_t);
figure(2)
atplot(ring_rad)

ring = ring_rad;
save fcch_rad_tapered.mat ring


