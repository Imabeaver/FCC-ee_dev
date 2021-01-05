% real implementation
clear

mod = py.importlib.import_module('Integrated_kick');

%import fcc ring
v = load('fcch_norad.mat');
ring = v.ring;

%generate at element
bb=atbaselem('BeambeamKick','bb_kick');
bb.lib = mod;

%track on particle through the ring
n_turns = 100;
rin = [1.0e-6 0 0.0e-6 0 0 0];
ring_new = [ring; {bb}];
rout = ringpass(ring_new, rin', n_turns);

figure(1)
plot(rout(1,:), rout(2,:), '.');
xlabel('x [m]');
ylabel('x^{prime} [rad]');

figure(2)
plot(rout(3,:), rout(4,:), '.');
xlabel('y [m]');
ylabel('y^{prime} [rad]');

%now test kick
x = (-30:30)*2.0e-6;
xp = zeros(length(x),1);
for i = 1: length(x)
    rin = [x(i) 0 0 0 0 0];
    kick = bb_kick(bb,rin');
    xp(i) = kick(2);
end

figure(3)
plot(x,xp)


%clear classes at the end to avoid matlab crash
%clear classes

