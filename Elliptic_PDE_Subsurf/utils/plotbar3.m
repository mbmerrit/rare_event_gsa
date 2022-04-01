function plotbar3(x1, y1, x2, y2, nx, ny)

if nargin < 8
   bdry = 0;
end

x = linspace(x1, x2, nx);
y = linspace(y1, y2, ny);
[X Y] = meshgrid(x, y);
xx = X(:);
yy = Y(:);

ww = ones(size(xx));

W = reshape(ww, ny, nx);
W = W(end : -1 : 1, :);
h = bar3(xx, yy, W);

for i = 1:numel(h)
  index = logical(kron(isnan(W(:,i)),ones(6,1)));
  zData = get(h(i),'ZData');
  zData(index,:) = nan;
  set(h(i),'ZData',zData);
end
