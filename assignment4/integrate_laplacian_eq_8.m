function p=integrate_laplacian_eq_8(lambda,k,q)
    % This function is the equivalent of Neelamani2006 eq (8):
    % integral of (lambda/2 * exp(-lambda * abs(tau))) for tau between (k-0.5)q and (k+0.5)q.
    % Because the abs() function is used, this can't be integrated directly
    % Instead, check the ends of the range.
    % If the whole range is above or below 0, abs(tau) = tau or -tau
    % If the range contains zero, split it into two ranges and sum up the results.
    
    if (k-0.5)*q >= 0
        % entire range above 0
        p = integrate_laplacian_positive(lambda, (k-0.5)*q, (k+0.5)*q);
    elseif (k+0.5)*q <= 0
        % entire range below 0
        p = integrate_laplacian_negative(lambda, (k-0.5)*q, (k+0.5)*q);
    else
        % range intersects 0
        % Assume k and q are positive
        p = integrate_laplacian_positive(lambda, 0, (k+0.5)*q);
        p = p + integrate_laplacian_negative(lambda, (k-0.5)*q, 0);
    end
end
function p=integrate_laplacian_positive(lambda, a, b)
    % integrated this one by hand
    p = 0.5*(-exp(-lambda*b)+exp(-lambda*a))*10;
end
function p=integrate_laplacian_negative(lambda, a, b)
    % integrated this one by hand
    p = 0.5*(exp(lambda*b)-exp(lambda*a))*10;
end
