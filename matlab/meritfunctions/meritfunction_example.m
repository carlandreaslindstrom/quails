% MERIT FUNCTION AUTO-GENERATED BY MATHEMATICA 
function f = meritfunction_example(x)
	I1 = x(1);
	I2 = x(2);
	d2 = x(3);

	f = ((-0.434783E2).*(0.114929E-3.*I1.^(-1).*I2.^(-1/2).*(0.2E1.*I1.^( ...
  3/2).*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+ ...
  0.1E1.*I1.^(1/2).*I2.*cos(0.316328E-1.*I1.^(1/2)).*cos( ...
  0.158164E-1.*I2.^(1/2))+0.1E1.*d2.*I1.^(1/2).*I2.*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.65312E2.* ...
  I1.*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+ ...
  0.104322E1.*I2.*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.* ...
  I1.^(1/2))+(-0.958569E0).*d2.*I1.*I2.*cos(0.158164E-1.*I2.^(1/2)) ...
  .*sin(0.316328E-1.*I1.^(1/2))+0.663552E2.*I1.^(1/2).*I2.^(1/2).* ...
  cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+ ...
  0.632256E2.*d2.*I1.^(1/2).*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)) ...
  .*sin(0.158164E-1.*I2.^(1/2))+(-0.191714E1).*d2.*I1.^(3/2).*I2.^( ...
  1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+( ...
  -0.3E1).*I1.*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+(-0.123212E3).*d2.*I1.*I2.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))).*(( ...
  -0.163247E1).*I1.^(1/2).*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).* ...
  cos(0.158164E-1.*I2.^(1/2))+0.1E1.*d2.*I1.^(3/2).*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.344046E2.* ...
  I2.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2) ...
  )+0.156483E1.*I1.*I2.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+0.32656E2.*d2.*I1.*I2.^(1/2).*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.344046E2.* ...
  I1.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2) ...
  )+0.104322E1.*I1.^(3/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.521611E0.*I1.^(1/2).*I2.*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.521611E0.* ...
  d2.*I1.^(1/2).*I2.*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.* ...
  I2.^(1/2))+0.108831E1.*I1.*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.544156E0.*I2.*sin(0.316328E-1.*I1.^(1/2) ...
  ).*sin(0.158164E-1.*I2.^(1/2))+(-0.5E0).*d2.*I1.*I2.*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2)))+ ...
  0.299742E-4.*I1.^(-1/2).*((-0.636429E3).*I1.^(3/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+(-0.321344E3) ...
  .*I1.^(1/2).*I2.*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.* ...
  I2.^(1/2))+(-0.316128E3).*d2.*I1.^(1/2).*I2.*cos(0.316328E-1.* ...
  I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.383428E1.*d2.*I1.^(3/2) ...
  .*I2.*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.123882E5).*I1.*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.* ...
  I1.^(1/2))+0.4E1.*I1.^2.*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+(-0.197875E3).*I2.*cos(0.158164E-1.*I2.^( ...
  1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.5E1.*I1.*I2.*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.428242E3.* ...
  d2.*I1.*I2.*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^( ...
  1/2))+(-0.12586E5).*I1.^(1/2).*I2.^(1/2).*cos(0.316328E-1.*I1.^( ...
  1/2)).*sin(0.158164E-1.*I2.^(1/2))+(-0.119924E5).*d2.*I1.^(1/2).* ...
  I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2) ...
  )+0.8E1.*I1.^(3/2).*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.852485E3.*d2.*I1.^(3/2).*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.1E1.*I1.^( ...
  1/2).*I2.^(3/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.* ...
  I2.^(1/2))+0.1E1.*d2.*I1.^(1/2).*I2.^(3/2).*cos(0.316328E-1.*I1.^( ...
  1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.956729E3.*I1.*I2.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.312871E5.* ...
  d2.*I1.*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.* ...
  I2.^(1/2))+(-0.383428E1).*d2.*I1.^2.*I2.^(1/2).*sin(0.316328E-1.* ...
  I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.104322E1.*I2.^(3/2).* ...
  sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+( ...
  -0.958569E0).*d2.*I1.*I2.^(3/2).*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))).*(0.E-307+cos(0.316328E-1.*I1.^(1/2)).* ...
  cos(0.158164E-1.*I2.^(1/2))+0.104322E1.*I1.^(-1/2).*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+cos( ...
  0.316328E-1.*I1.^(1/2)).*(d2.*cos(0.158164E-1.*I2.^(1/2))+ ...
  0.104322E1.*I2.^(-1/2).*sin(0.158164E-1.*I2.^(1/2)))+(-0.958569E0) ...
  .*I1.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*(d2.*cos(0.158164E-1.* ...
  I2.^(1/2))+0.104322E1.*I2.^(-1/2).*sin(0.158164E-1.*I2.^(1/2))))+ ...
  0.119897E-3.*I1.^(-1/2).*I2.^(-1/2).*(0.516068E2.*I1.^(1/2).*I2.^( ...
  1/2).*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.208644E1).*I1.^(3/2).*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).* ...
  cos(0.158164E-1.*I2.^(1/2))+(-0.159107E3).*d2.*I1.^(3/2).*I2.^( ...
  1/2).*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.260805E0).*I1.^(1/2).*I2.^(3/2).*cos(0.316328E-1.*I1.^(1/2)).* ...
  cos(0.158164E-1.*I2.^(1/2))+(-0.260805E0).*d2.*I1.^(1/2).*I2.^( ...
  3/2).*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.108762E4).*I2.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+(-0.150583E3).*I1.*I2.^(1/2).*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+(-0.309704E4) ...
  .*d2.*I1.*I2.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+0.1E1.*d2.*I1.^2.*I2.^(1/2).*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+(-0.272078E0) ...
  .*I2.^(3/2).*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^( ...
  1/2))+0.25E0.*d2.*I1.*I2.^(3/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+(-0.108762E4).*I1.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+(-0.100026E3) ...
  .*I1.^(3/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^( ...
  1/2))+(-0.508291E2).*I1.^(1/2).*I2.*cos(0.316328E-1.*I1.^(1/2)).* ...
  sin(0.158164E-1.*I2.^(1/2))+(-0.494687E2).*d2.*I1.^(1/2).*I2.*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.1E1.*d2.* ...
  I1.^(3/2).*I2.*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^( ...
  1/2))+(-0.344046E2).*I1.*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.104322E1.*I1.^2.*sin(0.316328E-1.*I1.^( ...
  1/2)).*sin(0.158164E-1.*I2.^(1/2))+(-0.172023E2).*I2.*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.130403E1.* ...
  I1.*I2.*sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+ ...
  0.800752E2.*d2.*I1.*I2.*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))).*(0.E-307+(-0.958569E0).*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+(-0.1E1).* ...
  I1.^(-1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+cos(0.316328E-1.*I1.^(1/2)).*(cos( ...
  0.158164E-1.*I2.^(1/2))+(-0.958569E0).*d2.*I2.^(1/2).*sin( ...
  0.158164E-1.*I2.^(1/2)))+(-0.958569E0).*I1.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*(cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.958569E0).*d2.*I2.^(1/2).*sin(0.158164E-1.*I2.^(1/2)))))+( ...
  -0.23E-1).*(0.299742E-4.*((-0.632256E3).*I1.*cos(0.316328E-1.* ...
  I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+(-0.316128E3).*I2.*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.383428E1.* ...
  d2.*I1.*I2.*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^( ...
  1/2))+(-0.119924E5).*I1.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+0.4E1.*I1.^(3/2).*cos(0.158164E-1.*I2.^( ...
  1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.5E1.*I1.^(1/2).*I2.*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.424242E3.* ...
  d2.*I1.^(1/2).*I2.*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.* ...
  I1.^(1/2))+(-0.119924E5).*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).* ...
  sin(0.158164E-1.*I2.^(1/2))+0.8E1.*I1.*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.848485E3.* ...
  d2.*I1.*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.* ...
  I2.^(1/2))+0.1E1.*I2.^(3/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.948383E3.*I1.^(1/2).*I2.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.306548E5.* ...
  d2.*I1.^(1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+(-0.383428E1).*d2.*I1.^(3/2).*I2.^(1/2).* ...
  sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+( ...
  -0.958569E0).*d2.*I1.^(1/2).*I2.^(3/2).*sin(0.316328E-1.*I1.^(1/2) ...
  ).*sin(0.158164E-1.*I2.^(1/2))).*(0.E-307+cos(0.316328E-1.*I1.^( ...
  1/2)).*cos(0.158164E-1.*I2.^(1/2))+(-0.958569E0).*I1.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*(d2.*cos(0.158164E-1.*I2.^(1/2))+ ...
  0.104322E1.*I2.^(-1/2).*sin(0.158164E-1.*I2.^(1/2))))+ ...
  0.119897E-3.*I2.^(-1/2).*((-0.208644E1).*I1.*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+(-0.158064E3) ...
  .*d2.*I1.*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*cos( ...
  0.158164E-1.*I2.^(1/2))+(-0.260805E0).*I2.^(3/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+(-0.148406E3) ...
  .*I1.^(1/2).*I2.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+(-0.29981E4).*d2.*I1.^(1/2).*I2.^(1/2).* ...
  cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.1E1.* ...
  d2.*I1.^(3/2).*I2.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+0.25E0.*d2.*I1.^(1/2).*I2.^(3/2).*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+(-0.989374E2) ...
  .*I1.*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+( ...
  -0.494687E2).*I2.*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.* ...
  I2.^(1/2))+0.1E1.*d2.*I1.*I2.*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+(-0.443809E-12).*I1.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.104322E1.* ...
  I1.^(3/2).*sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2) ...
  )+0.130403E1.*I1.^(1/2).*I2.*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.79032E2.*d2.*I1.^(1/2).*I2.*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))).*(0.E-307+( ...
  -0.958569E0).*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+(-0.958569E0).*I1.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*(cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.958569E0).*d2.*I2.^(1/2).*sin(0.158164E-1.*I2.^(1/2))))+ ...
  0.151611E-1.*I2.^(-1/2).*(0.1E1.*d2.*I1.*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.156483E1.* ...
  I1.^(1/2).*I2.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+0.316128E2.*d2.*I1.^(1/2).*I2.^(1/2).*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.104322E1.* ...
  I1.*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+ ...
  0.521611E0.*I2.*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.* ...
  I2.^(1/2))+(-0.5E0).*d2.*I1.^(1/2).*I2.*sin(0.316328E-1.*I1.^(1/2) ...
  ).*sin(0.158164E-1.*I2.^(1/2))).*(0.758055E-2.*I2.*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.479284E0.* ...
  I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2) ...
  )+(-0.151611E-1).*I1.^(1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2) ...
  ).*sin(0.158164E-1.*I2.^(1/2))+(-0.958569E0).*(((-0.158164E-1).* ...
  I1.*cos(0.316328E-1.*I1.^(1/2))+(-1/2).*I1.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2))).*(cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.958569E0).*d2.*I2.^(1/2).*sin(0.158164E-1.*I2.^(1/2)))+I1.^( ...
  1/2).*sin(0.316328E-1.*I1.^(1/2)).*(0.758055E-2.*d2.*I2.*cos( ...
  0.158164E-1.*I2.^(1/2))+0.790819E-2.*I2.^(1/2).*sin(0.158164E-1.* ...
  I2.^(1/2))+0.479284E0.*d2.*I2.^(1/2).*sin(0.158164E-1.*I2.^(1/2))) ...
  )))).^2;

end
