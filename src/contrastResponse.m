function [num,den]=contrastResponse(Rmax,eps,contrasts)
%% contrastResponse Tuning curve for standard contrast-response hyperbolic ratio function
%   Used for parametrizing the mean in the Ratio of Gaussians model
%   R = Rmax*c^2/(eps^2+c^2)
%   INPUT -
%       Rmax (2,1) - Maximum firing rate
%       eps (2,1) - Semi-saturation constant
%       contrasts (1,number of contrasts) - vector of contrast values
%   OUTPUT:
%       num  (2,number of contrasts) - numerator, Rmax*c^2
%       den (2,number of contrasts) - denominator, eps^2+c^2
%%
    arguments
        Rmax (2,1);
        eps (2,1);
        contrasts (1,:);
    end
    c2 = contrasts.*contrasts;
    num= Rmax.*c2;
    den = (eps.*eps+c2);
end
