% We use this function to make the MXM derivative with respect to vectorized X.

function f = Gamma2(H)
  n = size(H);
  for i1 = 1:n(2)
    for i2 = i1:n(2)
      f(:,(i1-1)*n(2)+i2-i1*(i1-1)/2) = vec_var(H(:,i1)*H(:,i2)')...
      +vec_var(H(:,i2)*H(:,i1)');
      if i1 == i2
        f(:,(i1-1)*n(2)+i2-i1*(i1-1)/2) = f(:,(i1-1)*n(2)+i2-i1*(i1-1)/2)/2;
      end
    end
  end
      end