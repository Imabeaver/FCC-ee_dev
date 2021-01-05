clear

mod = py.importlib.import_module('Integrated_kick');

v = load('fcch_norad.mat');
ring_norad = v.ring;

v = load('fcch_rad.mat');
ring_rad = v.ring;

%Get 6D orbit no bb / zero w.o radiations
co = findorbit6(ring_rad,1:length(ring_rad)+1);
cox = co(1,:);
cos = co(5,:);
spos = findspos(ring_rad,1:length(ring_rad)+1);

figure(1)
plot(spos,cos,spos,cox);
legend('dp/p', 'x')

%Compare tune
%track on particle through the ring
n_turns = 1024;
rin = [0.0 0 1.0e-6 0 0 0];
rout_norad = ringpass(ring_norad, rin', n_turns);
rout_rad = ringpass(ring_rad, rin', n_turns);

%Compare tune bb

%generate at element
bb=atbaselem('BeambeamKick','bb_kick');
bb.lib = mod;

ring_norad_bb = [ring_norad; {bb}];
ring_rad_bb = [ring_rad; {bb}];

rout_norad_bb = ringpass(ring_norad_bb, rin', n_turns);
rout_rad_bb = ringpass(ring_rad_bb, rin', n_turns);

figure(2)
plot((0:n_turns-1)/(n_turns-1),abs(fft(rout_norad(3,:))))
hold on
plot((0:n_turns-1)/(n_turns-1),abs(fft(rout_rad(3,:))))
plot((0:n_turns-1)/(n_turns-1),abs(fft(rout_norad_bb(3,:))))
plot((0:n_turns-1)/(n_turns-1),abs(fft(rout_rad_bb(3,:))))
legend('rad off / bb off','rad on / bb off',...
       'rad off / bb on','rad on / bb on')
hold off



