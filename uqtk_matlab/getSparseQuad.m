function quadrature = getSparseQuad(accuracyLevel, ndim);

type = 'KPU';
[nodes, weights] = nwspgr(type, ndim, accuracyLevel);

quadrature.nodes = -1 + 2 * nodes;         % need to map from [0, 1] to [-1, 1]
quadrature.weights = weights;
quadrature.nquad = length(weights);
