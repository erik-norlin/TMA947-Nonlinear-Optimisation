function df = holmg(x)
%df = [6*x(1)^5+8*x(1)^3-30*x(1)^2*x(2)-8*x(2)^3+24*x(1)*x(2)^2; 8*x(2)-10*x(1)^3-24*x(1)*x(2)^2+24*x(1)^2*x(2)+2*x(2)];
df = [6*x(1)^2*(x(1)^3-x(2)) + 8*(x(1)-x(2))^3; -2*(x(1)^3-x(2)) - 8*(x(1)-x(2))^3];
