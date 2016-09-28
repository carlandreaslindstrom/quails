% CALCULATES NORMALIZED EMITTANCE GIVEN A BEAM
function [ emit_Nx, emit_Ny, beta_x, beta_y, alpha_x, alpha_y, E, sigma_E] = beamStatistics( beam )
    
    % calculate geometric emittances
    cx = cov(beam(:,1:2));
    cy = cov(beam(:,3:4));
    emit_gx = sqrt(det(cx));
    emit_gy = sqrt(det(cy));
    
    % energy and energy spread
    E = mean(beam(:,6));
    sigma_E = std(beam(:,6))/E;
    
    % normalized emittances
    gamma = E/0.511e-3;
    emit_Nx = emit_gx*gamma;
    emit_Ny = emit_gy*gamma;
    
    % alphas and betas
    beta_x = cx(1,1)/emit_gx;
    beta_y = cy(1,1)/emit_gy;
    alpha_x = -cx(2,1)/emit_gx;
    alpha_y = -cy(2,1)/emit_gy;
    
end

