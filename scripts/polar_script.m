% polar plot script

theta = linspace(0, 2*pi, 9);
rho = randi(8,[8,1]);
rho = [rho;rho(1)];

figure
polar(theta,rho')