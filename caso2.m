function z=caso2(y,r)
  f=@(z)norm(z-y,2)^2;
  h=@(z) norm(z,2)-r;
  n=length(y);
  z0=zeros(1,n);
  z0(1)=1;
  [z, obj, info, iter, nf, lambda] =sqp(z0,f,h,[],[],[], 300, 10*e^-17);
  endfunction