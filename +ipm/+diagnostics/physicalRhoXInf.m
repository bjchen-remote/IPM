function value = physicalRhoXInf(flow,scale)
%IPM.DIAGNOSTICS.PHYSICALRHOXINF Global physical max |rho_x|.
%   FLOW.SOURCE is the stored-coordinate derivative R_X. The isotropic and
%   anisotropic physical reconstructions both multiply it by C_x/C_omega.

value = exp(scale.logC_l-scale.logC_omega) * ...
    max(abs(flow.source),[],'all');
end
