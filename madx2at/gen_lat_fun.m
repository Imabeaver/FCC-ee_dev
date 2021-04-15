function gen_lat_fun(fname)

global NUMDIFPARAMS
NUMDIFPARAMS.XYStep = 1e-9;
NUMDIFPARAMS.DPStep = 1e-9;
clight = 2.99792458e8;

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

%volt = sum(atgetfieldvalues(ring(cav),'Voltage'));
freq = atgetfieldvalues(ring(cav),'Frequency');
[fc,~,ib] = unique(freq);
len = findspos(ring,length(ring)+1);
harm = floor(fc/clight*len);
for i = 1: length(cav)
    ring{cav(i)}.HarmNumber = harm(ib(i));
    ring{cav(i)}.Frequency = ring{cav(i)}.HarmNumber/len*clight;
end

idpass = findcells(ring,'PassMethod','IdentityPass');
ring(idpass) = atsetfieldvalues(ring(idpass),'Length',0);

%ring w.o. radiations
ring = atradoff(ring,'IdentityPass','auto','auto');
eval(['save ',fname,'_norad.mat ring']);
end