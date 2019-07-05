%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Çòºº¿Ë¶ûº¯Êý%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = SphHankel(n,k,x)

if(k==2)
    I = SphBesselj(n,x)-1j*SphBessely(n,x);
else
    I = SphBesselj(n,x)+1j*SphBessely(n,x);
end

end