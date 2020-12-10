% real implementation
clear

%import fcc ring
v = load('fcch_norad.mat');
ring = v.ring;

%generate at element
bb=atbaselem('BeambeamKick','bb_kick');

n_turns = 100;
rin = [1.0e-6 0 0 0 0 0];

fin_x = zeros(1, n_turns);
fin_xp = zeros(1, n_turns);

xp = zeros(1, n_turns);

figure; hold on

for i = 1:n_turns
    x = rin(1);
    y = rin(3); 
    
    fin_x(i) = rin(1);
    fin_xp(i) = rin(2);
    
    a = py.test.kick_calc(x,y);
    dx = a{1};
    dx = double(dx);
    dy = a{2};
    dy = double(dy);
    
    amp = [0 dx*1.0e-6 0 dy*1.0e-6 0 0]';
    bb.Amplitude = amp; 
    ring_new = [ring; {bb}];
    
    rout = ringpass(ring_new, rin', 1);
    rin = rout'; 
    
    xp(i) = dx;    
    
end
plot(fin_x, fin_xp, '.');
xlabel('x [m]');
ylabel('x^{prime} [rad]');
hold off

figure; 
[fin_x, indx] = sort(fin_x, 'descend');
xp = xp(indx);
plot(fin_x/6.4*1.0e6, xp, '.b');
xlabel('\sigma_x[m]');
ylabel('dx[\mu rad]');

