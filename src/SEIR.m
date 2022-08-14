function dy = SEIR(t,y,B,N,p,tk)
beta = interp1(0:1:tk,B,t);
dy = zeros(4,1);
dy(1) = (-beta*y(1)*y(3))/N;
dy(2) = (beta*y(1)*y(3))/N - p(1)*y(2);
dy(3) = p(1)*y(2) - p(2)*y(3);
dy(4) = p(2)*y(3);
end