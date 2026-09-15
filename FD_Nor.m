function FD = FD_Nor(x, R, R_ext)
% Dimensionless pressure for the Nordbotten solution.
% x can be scalar or vector. Returns 0 outside R_ext.
% If R > R_ext, uses the confined (shape-factor) branch.

    FD   = zeros(size(x));
    mask = x < R_ext;

    if any(mask)
        if R <= R_ext
            FD(mask) = log(R ./ x(mask));
        else
            FD(mask) = log(R_ext ./ x(mask)) + 2/2.25*(R/R_ext)^2 - 3/4;
        end
    end
end