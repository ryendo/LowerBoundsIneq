function invs = upper_matrix_trial_invariants(model, coeff, mode)
%UPPER_MATRIX_TRIAL_INVARIANTS Outward evaluation of four fixed-trial forms.
if strcmpi(mode,'interval')
    q = I_intval(coeff);
else
    q = coeff;
end
invs.Exx = q'*model.Kxx*q;
invs.Exy = q'*model.Kxy*q;
invs.Eyy = q'*model.Kyy*q;
invs.M0 = q'*model.M*q;
end
