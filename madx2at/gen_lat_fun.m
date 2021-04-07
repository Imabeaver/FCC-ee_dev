function gen_lat_fun(fname)

global NUMDIFPARAMS
NUMDIFPARAMS.XYStep = 1e-9;
NUMDIFPARAMS.DPStep = 1e-9;

v = load(strcat(fname,'.mat'));
ring = v.ring;

cav = findcells(ring,'Class','RFCavity');
for i=length(cav):-1:1
    if ring{cav(i)}.Voltage==0
        ring{cav(i)}=atmarker('CAV_marker');
    end
end

cav = findcells(ring,'Class','RFCavity');
dip = findcells(ring,'Class','Bend');
quad = findcells(ring,'Class','Quadrupole');
sext = findcells(ring,'Class','Sextupole');
ring(dip) = atsetfieldvalues(ring(dip),'NumIntSteps',20);
ring(quad) = atsetfieldvalues(ring(quad),'NumIntSteps',20);
ring(sext) = atsetfieldvalues(ring(sext),'NumIntSteps',20);

ring(dip) = atsetfieldvalues(ring(dip),'EntranceAngle',0);
ring(dip) = atsetfieldvalues(ring(dip),'ExitAngle',0);

volt = sum(atgetfieldvalues(ring(cav),'Voltage'));
freq = mean(atgetfieldvalues(ring(cav),'Frequency'));
len = findspos(ring,length(ring)+1);
harm = floor(freq/2.99792e8*len);

idpass = findcells(ring,'PassMethod','IdentityPass');
ring(idpass) = atsetfieldvalues(ring(idpass),'Length',0);

%ring w.o. radiations
ring = atsetcavity(ring,volt,0,harm);
ring = atradoff(ring,'IdentityPass','auto','auto');
eval(['save ',fname,'_norad.mat ring']);
end