% We define this function to transform cubic to rectangular matrices, when
% this is possible

function f = rect_D2(D2F)
  n = size(D2F);
  aux = zeros(n(1),n(1));
  for i = 1:n(1)
    aux(:,:) = D2F(i,:,:);
    f(i,:) = vec_var(aux)';
  end
  end