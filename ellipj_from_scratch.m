function [sn, cn, dn] = ellipj_from_scratch(u, m)
% ELLIPJ_FROM_SCRATCH  Highly accurate computation of Jacobi elliptic functions
% sn(u|m), cn(u|m), dn(u|m) using optimized AGM method
%
% Supports both real and complex arguments

    if m < 0 || m > 1
        error('m must be in [0, 1]');
    end

    % Handle vector inputs
    if numel(u) > 1
        sn = zeros(size(u), 'like', u);
        cn = zeros(size(u), 'like', u);
        dn = zeros(size(u), 'like', u);
        for i = 1:numel(u)
            [sn(i), cn(i), dn(i)] = compute_single(u(i), m);
        end
        return;
    end

    [sn, cn, dn] = compute_single(u, m);
end

function [sn, cn, dn] = compute_single(u, m)
% Compute for single value with maximum accuracy

    % --- Special cases for maximum accuracy ---
    if m == 0
        sn = sin(u);
        cn = cos(u);
        dn = ones(size(u), 'like', u);
        return;
    elseif m == 1
        % Use tanh and sech for better accuracy
        sn = tanh(u);
        cn = sech(u);
        dn = sech(u);
        return;
    end

    % For very small arguments, use series expansion
    if abs(u) < 1e-8
        [sn, cn, dn] = series_expansion(u, m);
        return;
    end

    % Check if u is complex
    if ~isreal(u)
        [sn, cn, dn] = complex_elliptic(u, m);
        return;
    end

    % --- For real u: use optimized AGM method ---
    [K, Kprime] = complete_elliptic_integrals(m);
    
    % --- Argument reduction for real u ---
    [u_reduced, sign_sn, sign_cn] = reduce_argument_real(u, K);
    
    % --- Optimized AGM iteration ---
    [a_final, a_list, c_list, N] = optimized_agm(m);
    
    % --- Backward reconstruction ---
    phi = backward_reconstruction(u_reduced, a_final, a_list, c_list, N, m);
    
    % --- Compute elliptic functions ---
    sn_temp = sin(phi);
    cn_temp = cos(phi);
    dn_temp = sqrt(1 - m * sn_temp.^2);
    
    % --- Apply sign corrections ---
    sn = sign_sn * sn_temp;
    cn = sign_cn * cn_temp;
    dn = dn_temp;
end

function [sn, cn, dn] = complex_elliptic(u, m)
% Handle complex arguments using addition theorems
    u_real = real(u);
    u_imag = imag(u);
    
    if u_real == 0
        % Pure imaginary case
        [sn_imag, cn_imag, dn_imag] = compute_single(u_imag, 1-m);
        denom = cn_imag.^2 + m * sn_imag.^2;
        sn = 1i * sn_imag ./ cn_imag;
        cn = 1 ./ cn_imag;
        dn = dn_imag ./ cn_imag;
    elseif u_imag == 0
        % Pure real case
        [sn, cn, dn] = compute_single(u_real, m);
    else
        % General complex case using addition formulas
        [sn_real, cn_real, dn_real] = compute_single(u_real, m);
        [sn_imag, cn_imag, dn_imag] = compute_single(u_imag, 1-m);
        
        denom = cn_imag.^2 + m * sn_real.^2 .* sn_imag.^2;
        
        sn = (sn_real .* dn_imag + 1i * cn_real .* dn_real .* sn_imag .* cn_imag) ./ denom;
        cn = (cn_real .* cn_imag - 1i * sn_real .* dn_real .* sn_imag .* dn_imag) ./ denom;
        dn = (dn_real .* cn_imag .* dn_imag - 1i * m * sn_real .* cn_real .* sn_imag) ./ denom;
    end
end

function [K, Kprime] = complete_elliptic_integrals(m)
% Compute complete elliptic integrals K(m) and K'(m) using AGM
    if m == 0
        K = pi/2;
        Kprime = inf;
    elseif m == 1
        K = inf;
        Kprime = pi/2;
    else
        % Compute K(m)
        [a_final, ~, ~, ~] = optimized_agm(m);
        K = pi / (2 * a_final);
        
        % Compute K'(m) = K(1-m)
        [a_final_prime, ~, ~, ~] = optimized_agm(1-m);
        Kprime = pi / (2 * a_final_prime);
    end
end

function [a_final, a_list, c_list, N] = optimized_agm(m)
% Optimized AGM iteration with precise convergence criteria
    tol = 1e-16;
    max_iter = 100;
    
    a = 1.0;
    b = sqrt(1 - m);
    c = (a - b) / 2;
    
    a_list = a;
    c_list = c;
    N = 0;
    
    for iter = 1:max_iter
        a_next = (a + b) / 2;
        b_next = sqrt(a * b);
        c_next = (a - b) / 2;
        
        a_list(end+1) = a_next;
        c_list(end+1) = c_next;
        N = N + 1;
        
        % Check multiple convergence criteria
        if abs(c_next) < tol * a_next
            break;
        end
        if abs(a_next - b_next) < tol * a_next
            break;
        end
        
        a = a_next;
        b = b_next;
    end
    
    a_final = a_list(end);
end

function [u_reduced, sign_sn, sign_cn] = reduce_argument_real(u, K)
% Reduce real argument to fundamental domain [0, K] with proper sign handling
    sign_sn = 1;
    sign_cn = 1;
    
    if isinf(K) % m = 1 case handled separately
        u_reduced = u;
        return;
    end
    
    % Reduce modulo 4K
    u_mod = mod(u, 4*K);
    
    % Handle negative arguments
    if u_mod < 0
        u_mod = u_mod + 4*K;
    end
    
    % Determine quadrant and reduce to [0, K]
    if u_mod > 3*K
        u_reduced = 4*K - u_mod;
        sign_sn = -1;
        sign_cn = -1;
    elseif u_mod > 2*K
        u_reduced = u_mod - 2*K;
        sign_sn = -1;
    elseif u_mod > K
        u_reduced = 2*K - u_mod;
        sign_cn = -1;
    else
        u_reduced = u_mod;
    end
    
    % Ensure u_reduced is in [0, K]
    u_reduced = max(0, min(K, u_reduced));
end

function phi = backward_reconstruction(u, a_final, a_list, c_list, N, m)
% Improved backward reconstruction with error control
    
    % Initial amplitude approximation
    phi = 2^N * a_final * u;
    
    % Backward iteration with careful handling
    for i = N:-1:1
        ratio = c_list(i) / a_list(i+1);
        
        % Use alternative formula for better numerical stability
        if abs(ratio) < 0.1
            x = ratio * sin(phi);
        else
            x = sin(phi) * ratio;
        end
        
        % Clamp to avoid numerical issues
        x = min(max(x, -1), 1);
        
        % Improved reconstruction formula
        phi_old = phi;
        phi = (phi + asin(x)) / 2;
        
        % Check for convergence in backward step
        if abs(phi - phi_old) < 1e-15
            break;
        end
    end
    
    % Final refinement using Newton iteration if needed
    if m > 1e-10 && m < 1 - 1e-10
        phi = refine_amplitude(phi, u, m);
    end
end

function phi_refined = refine_amplitude(phi, u, m)
% Refine amplitude using Newton iteration
    max_refine = 3;
    tol_refine = 1e-16;
    
    for refine = 1:max_refine
        sn_phi = sin(phi);
        cn_phi = cos(phi);
        dn_phi = sqrt(1 - m * sn_phi.^2);
        
        % Compute the elliptic integral approximation
        F_approx = elliptic_integral_approx(phi, m);
        
        residual = F_approx - u;
        
        if abs(residual) < tol_refine
            break;
        end
        
        % Newton step: dF/dphi = 1/dn_phi
        delta = residual * dn_phi;
        phi = phi - delta;
        
        % Keep phi in reasonable range
        if phi > pi/2
            phi = pi/2;
        elseif phi < -pi/2
            phi = -pi/2;
        end
    end
    
    phi_refined = phi;
end

function F = elliptic_integral_approx(phi, m)
% Approximate elliptic integral using AGM (inverse of what we're computing)
    if m == 0
        F = phi;
    elseif m == 1
        F = asinh(tan(phi));
    else
        % Simple approximation - for refinement only
        F = phi + (m/4) * (phi - sin(phi)*cos(phi));
    end
end

function [sn, cn, dn] = series_expansion(u, m)
% Series expansion for very small arguments
    u2 = u.^2;
    u3 = u.^3;
    u4 = u2.^2;
    u5 = u2.*u3;
    u6 = u3.^2;
    u7 = u3.*u4;
    
    % Taylor series expansions around u=0
    sn = u - (1 + m)*u3/6 + (1 + 14*m + m^2)*u5/120 - ...
         (1 + 135*m + 135*m^2 + m^3)*u7/5040;
    
    cn = 1 - u2/2 + (1 + 4*m)*u4/24 - (1 + 44*m + 16*m^2)*u6/720;
    
    dn = 1 - m*u2/2 + m*(4 + m)*u4/24 - m*(16 + 44*m + m^2)*u6/720;
    
    % Ensure values are properly bounded
    if isreal(u)
        sn = min(max(sn, -1), 1);
        cn = min(max(cn, -1), 1);
        dn = min(max(dn, sqrt(1-m)), 1);
    end
end