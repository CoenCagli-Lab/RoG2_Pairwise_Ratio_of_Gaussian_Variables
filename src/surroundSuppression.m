function [num,den]=surroundSuppression(Rmax,eps,contrasts)
    %% surroundSuppresion Tuning curve for a pseudo-surround suppression model
%   Used for parametrizing the mean in the Ratio of Gaussians model
%   INPUT -
%       Rmax (2,1) - Maximum firing rate
%       eps (2,1) - Semi-saturation constant
%       contrasts (1,number of contrasts) - vector of contrast values
%   OUTPUT:
%       num  (2,number of contrasts) 
%       den (2,number of contrasts) 
    arguments
        Rmax (2,1);
        eps (2,1);
        contrasts (1,:);
    end   
    CONST_WC = 25; 
    c2 = contrasts.^2;
    c2m = min(c2,CONST_WC.^2);
    num= Rmax.*c2m;
    den = (eps.*eps+c2m+max(c2-CONST_WC.^2,0));
end

