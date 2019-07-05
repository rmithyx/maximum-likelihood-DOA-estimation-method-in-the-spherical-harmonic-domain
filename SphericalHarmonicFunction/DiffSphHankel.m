%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%球汉克尔函数导数%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = DiffSphHankel(n,k,x)

if(k==2)
    I = DiffSphBesselj(n,x)-1j*DiffSphBessely(n,x);
else
    I = DiffSphBesselj(n,x)+1j*DiffSphBessely(n,x);
end

end