function S = MSFitTiling (N, ratio)
%MSFITTILING Summary of this function goes here
%   Detailed explanation goes here

p = sqrt(N/ratio);
q = N/p;
r1 = ceil(p);
c1 = ceil(N/r1);
c2 = ceil(q);
r2 = ceil(N/c2);

if r1*c1 < r2*c2
  select = 1;
elseif r1*c1 > r2*c2
  select = 2;
elseif abs(c1/r1-ratio) <= abs(c2/r2-ratio)
  select = 1;
else
  select = 2;
end

if select == 1
  S = [r1 c1];
else
  S = [r2 c2];
end

end

