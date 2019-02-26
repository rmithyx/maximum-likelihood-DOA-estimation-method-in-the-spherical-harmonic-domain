%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%¼ÆËãYnm£¬m´Ó-nµ½n%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = SphHarmonic(n,theta,phi)
    
    P=legendre(n,cos(theta));
    y=zeros(2*n+1,1);
    for m=-n:n
        temp=sqrt((2*n+1)/(4*pi) * factorial(n-abs(m))/factorial(n+abs(m))) * exp(1i*m*phi);
        if (m>=0)
            I(m+n+1,1)=temp * P(m+1);
        else
            fm=-m;
            I(m+n+1,1)=temp * P(fm+1) * (-1)^fm;
        end    
    end
    
end