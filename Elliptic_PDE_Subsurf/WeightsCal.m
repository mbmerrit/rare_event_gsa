function w = WeightsCal(nx, ny)

dx = 2/nx;
dy = 1/ny;

x = [0 : dx : 2]';
y = [0 : dy : 1]';

wx = 2*ones(size(x));
wx(1) = 1;
wx(end) = 1;
wx = 0.5 * dx * wx;

wy = 2*ones(size(y));
wy(1) = 1;
wy(end) = 1;
wy = 0.5 * dy * wy;

w = zeros(nx+1,ny+1);
for i = 1 : nx+1
   for j = 1 : ny+1
      w(i, j) = wx(i) * wy(j);
   end
end

