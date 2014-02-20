amat
a=mat;
aat = a*a';

umatsvd
usvd=mat;

dmatsvd
dsvd=diag(mat);

vmatsvd
vsvd=mat;

umatschur
uschur=mat;

dmatschur
dschur=diag(mat);

zerotol = 1e-10;
mzerotol = -zerotol;

udiff = -uschur-usvd;
[r c] = size(udiff);
for i=1:r
    for j=1:c
        if udiff(i,j) > mzerotol && udiff(i,j) < zerotol
            udiff(i,j) = 0;
        end
    end
end



max(max(abs(usvd*dsvd*vsvd-a)))
