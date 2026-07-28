function [coeff, lambda_ritz] = upper_matrix_trial(model, x, y)
%UPPER_MATRIX_TRIAL Generate the first Ritz vector in a fixed high-order space.
K0 = (y^2+x^2)*I_mid(model.Kxx) ...
    -2*x*I_mid(model.Kxy)+I_mid(model.Kyy);
M0 = y^2*I_mid(model.M);
try
    opts.tol = 1e-12;
    opts.maxit = 3000;
    [v,d] = eigs(K0,M0,1,'smallestreal',opts);
    lambda_ritz = real(d(1,1));
catch
    try
        [v,d] = eigs(K0,M0,1,'sm');
        lambda_ritz = real(d(1,1));
    catch
        [V,D] = eig(full(K0),full(M0));
        values = real(diag(D));
        values(values <= 0 | ~isfinite(values)) = Inf;
        [lambda_ritz,idx] = min(values);
        v = V(:,idx);
    end
end
coeff = real(v);
coeff = coeff/max(abs(coeff));
end
