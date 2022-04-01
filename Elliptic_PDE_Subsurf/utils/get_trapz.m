function w = get_trapz(xmesh)
dx = xmesh(2)-xmesh(1);
w = ones(size(xmesh))*2;
w(1) = 1;
w(end) = 1;
w = 0.5 * dx * w;
