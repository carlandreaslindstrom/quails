% CONSTRAINTS AUTO-GENERATED BY MATHEMATICA 
function [c, ceq] = constraints_example(x)
	I1 = x(1);
	I2 = x(2);
	d2 = x(3);

	ceq1 = (-0.23E-1).*(0.E-307+cos(0.316328E-1.*I1.^(1/2)).*(0.E-307+cos( ...
  0.158164E-1.*I2.^(1/2)))+(-0.958569E0).*I1.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*(d2.*cos(0.158164E-1.*I2.^(1/2))+ ...
  0.104322E1.*I2.^(-1/2).*sin(0.158164E-1.*I2.^(1/2)))).*(0.E-307+( ...
  -0.958569E0).*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+(-0.958569E0).*I1.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*(cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.958569E0).*d2.*I2.^(1/2).*sin(0.158164E-1.*I2.^(1/2))))+( ...
  -0.434783E2).*(0.E-307+cos(0.316328E-1.*I1.^(1/2)).*(0.E-307+cos( ...
  0.158164E-1.*I2.^(1/2)))+0.104322E1.*I1.^(-1/2).*(0.E-307+cos( ...
  0.158164E-1.*I2.^(1/2))).*sin(0.316328E-1.*I1.^(1/2))+cos( ...
  0.316328E-1.*I1.^(1/2)).*(d2.*cos(0.158164E-1.*I2.^(1/2))+ ...
  0.104322E1.*I2.^(-1/2).*sin(0.158164E-1.*I2.^(1/2)))+(-0.958569E0) ...
  .*I1.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*(d2.*cos(0.158164E-1.* ...
  I2.^(1/2))+0.104322E1.*I2.^(-1/2).*sin(0.158164E-1.*I2.^(1/2)))).* ...
  (0.E-307+(-0.958569E0).*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).* ...
  sin(0.158164E-1.*I2.^(1/2))+(-0.1E1).*I1.^(-1/2).*I2.^(1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+cos( ...
  0.316328E-1.*I1.^(1/2)).*(cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.958569E0).*d2.*I2.^(1/2).*sin(0.158164E-1.*I2.^(1/2)))+( ...
  -0.958569E0).*I1.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*(cos( ...
  0.158164E-1.*I2.^(1/2))+(-0.958569E0).*d2.*I2.^(1/2).*sin( ...
  0.158164E-1.*I2.^(1/2))));

	ceq2 = (-0.23E-1).*((cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^( ...
  1/2))+(-0.958569E0).*d2.*I1.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).* ...
  sin(0.316328E-1.*I1.^(1/2))+(-0.1E1).*I1.^(1/2).*I2.^(-1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))).*( ...
  0.151611E-1.*I1.*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.* ...
  I2.^(1/2))+0.758055E-2.*I2.*cos(0.316328E-1.*I1.^(1/2)).*cos( ...
  0.158164E-1.*I2.^(1/2))+0.479284E0.*I1.^(1/2).*cos(0.158164E-1.* ...
  I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+(-0.726648E-2).*d2.*I1.^( ...
  1/2).*I2.*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2)) ...
  +0.479284E0.*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+(-0.14533E-1).*d2.*I1.*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+( ...
  -0.227416E-1).*I1.^(1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).* ...
  sin(0.158164E-1.*I2.^(1/2))+(-0.918854E0).*d2.*I1.^(1/2).*I2.^( ...
  1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2)))+( ...
  0.151611E-1.*d2.*I1.*cos(0.316328E-1.*I1.^(1/2)).*cos( ...
  0.158164E-1.*I2.^(1/2))+0.237246E-1.*I1.^(1/2).*cos(0.158164E-1.* ...
  I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.479284E0.*d2.*I1.^(1/2) ...
  .*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+ ...
  0.158164E-1.*I1.*I2.^(-1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.790819E-2.*I2.^(1/2).*cos(0.316328E-1.* ...
  I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+(-0.758055E-2).*d2.*I1.^( ...
  1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.* ...
  I2.^(1/2))).*((-0.958569E0).*I1.^(1/2).*cos(0.158164E-1.*I2.^(1/2) ...
  ).*sin(0.316328E-1.*I1.^(1/2))+(-0.958569E0).*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.918854E0.* ...
  d2.*I1.^(1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))))+(-0.434783E2).*((0.E-307+cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+d2.*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.104322E1.* ...
  I1.^(-1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^( ...
  1/2))+(-0.958569E0).*d2.*I1.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).* ...
  sin(0.316328E-1.*I1.^(1/2))+0.104322E1.*I2.^(-1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+(-0.1E1).* ...
  I1.^(1/2).*I2.^(-1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))).*(0.151611E-1.*I1.*cos(0.316328E-1.*I1.^( ...
  1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.758055E-2.*I2.*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.758055E-2.* ...
  d2.*I2.*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+ ...
  0.495101E0.*I1.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+0.790819E-2.*I1.^(-1/2).*I2.*cos( ...
  0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+( ...
  -0.726648E-2).*d2.*I1.^(1/2).*I2.*cos(0.158164E-1.*I2.^(1/2)).* ...
  sin(0.316328E-1.*I1.^(1/2))+0.503009E0.*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.479284E0.* ...
  d2.*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^( ...
  1/2))+(-0.14533E-1).*d2.*I1.*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2) ...
  ).*sin(0.158164E-1.*I2.^(1/2))+(-0.227416E-1).*I1.^(1/2).*I2.^( ...
  1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+( ...
  -0.934016E0).*d2.*I1.^(1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2) ...
  ).*sin(0.158164E-1.*I2.^(1/2)))+(0.E-307+(-0.2475E-1).*cos( ...
  0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+0.151611E-1.* ...
  d2.*I1.*cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+ ...
  0.521611E0.*I1.^(-1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+0.237246E-1.*I1.^(1/2).*cos(0.158164E-1.* ...
  I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+0.495101E0.*d2.*I1.^(1/2) ...
  .*cos(0.158164E-1.*I2.^(1/2)).*sin(0.316328E-1.*I1.^(1/2))+ ...
  0.521611E0.*I2.^(-1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.158164E-1.*I1.*I2.^(-1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.790819E-2.* ...
  I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2) ...
  )+0.790819E-2.*d2.*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+0.165E-1.*I1.^(1/2).*I2.^(-1/2).*sin( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.825E-2.* ...
  I1.^(-1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin( ...
  0.158164E-1.*I2.^(1/2))+(-0.758055E-2).*d2.*I1.^(1/2).*I2.^(1/2).* ...
  sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))).*( ...
  0.E-307+cos(0.316328E-1.*I1.^(1/2)).*cos(0.158164E-1.*I2.^(1/2))+( ...
  -0.958569E0).*I1.^(1/2).*cos(0.158164E-1.*I2.^(1/2)).*sin( ...
  0.316328E-1.*I1.^(1/2))+(-0.958569E0).*I2.^(1/2).*cos( ...
  0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2))+(-0.958569E0) ...
  .*d2.*I2.^(1/2).*cos(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.* ...
  I2.^(1/2))+(-0.1E1).*I1.^(-1/2).*I2.^(1/2).*sin(0.316328E-1.*I1.^( ...
  1/2)).*sin(0.158164E-1.*I2.^(1/2))+0.918854E0.*d2.*I1.^(1/2).* ...
  I2.^(1/2).*sin(0.316328E-1.*I1.^(1/2)).*sin(0.158164E-1.*I2.^(1/2) ...
  )));

	c = [];
	ceq = [ceq1 ceq2 ];
end