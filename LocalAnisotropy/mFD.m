function FD = mFD(xi, a, r, d)

% get needed Tensor U

u13fu23f = U13fU23f(xi);

u13f = U13f(xi);
u23f = U23f(xi);
u13tu23t = U13tU23t(xi);
u13t = U13t(xi);
u23t = U23t(xi);
u13fu23t = U13fU23t(xi);
u13tu23f = U13tU23f(xi);



u33f = U33f(xi);

u13tu33t = U13tU33t(xi);
u23tu33t = u13tu33t;

u13tu23t = U13tU23t(xi);
u13t = U13t(xi);
u23t = u13t;
u33t = U33t(xi);

u13tu23tu33f = U13tU23tU33f(xi);
u13tu23tu33t = U13tU23tU33t(xi);

u13fu33t = U13fU33t(xi);
u23fu33t = u13fu33t;


FD = a^2 + r^2*u13tu23tu33f + d^2*(u33f + u13tu23t + u13tu33t + u23tu33t) + ...
    a*r*(u13tu33t + u23tu33t) + a*d*(u23t + u13t + 2*u33t) + ...
    r*d*(u13fu33t+ 2*u13tu23tu33t + u23fu33t);
end