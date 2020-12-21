% Vectorize a matrix of variances and covariances

function f =vec_var(M)
  n = size(M);
  G = M(1,:);
  for i = 2:n(2)
    g =  M(i,i:n(1));
    G = [G g];
  end
  f = G';
  end