function qs = identify_qs(d, h)
    qs_fft_est = estimate_fft_qs(h)
    
    % Imagine a dirac comb with period T, sampled between T = [0, 1].
    % The FFT of that sample will have a period F = 1/T
    % *If* the histogram values were sampled between T = {0, 1}, 
    % then T = 1/F
    % but they're not, the histogram values are over the range X = {1,
    % length(h)}
    % => T = length(h)
    % e.g. length(h)/qs_fft_est
    qs_est = round(length(h)./qs_fft_est)
    
    REFINE_ESTIMATE = 1;
    
    if ~REFINE_ESTIMATE || isempty(qs_est)
        qs = qs_est;
    else
        % Take Neelamani2006 approach to refining estimated Qs
        % They treat d as the 'set D of the observed coefficients'
        % This would seemingly imply to remove duplicate elements,
        % but I don't think that's correct.
        % More common coefficient values should have a larger effect on the score.
        % By inspection, this produced more sensible values on test2.png
%         d = unique(d);
        N = length(d);
        lambda = N/sum(abs(d));

        qs = [];
        % Refine each q_est in qs_est
        for q_est = qs_est
            % As in [Fan2003], choose possible parameter space from q_est
            % Q, Q+1, Q-1.
            % We don't consider integer factors here: a) those would have shown up in the frequency
            % distribution, b) the selection process is unfairly biased towards smaller factors, c)
            % it makes detecting multiple factors more challenging (e.g. if qs_est = [40 5], we
            % could estimate [5 5]).
            qs_poss = [q_est, q_est-1, q_est+1];

            % Neelamani2006 equations:
            % (5) Xtilde (DCT coeff) = Xbar_q (quantized coeff) + round-off error
            %
            % (6) P(roundoff = t) = gaussian if t in range [-6, 6], 0 otherwise
            % (7) P(Xtilde = t) = laplacian, lambda estimated from observed image
            %
            % (8) P(X_q = k*q | q) = integral of (7) over (k+0.5)q, (k-0.5)q
            % (e.g. the probability of a quantized result being some integer k * q given q)
            %
            % (11) P(Xtilde = t|q) = sum(P(roundoff=t-kq)*P(X_q = kq | q) for each possible offset k)
            %
            % (14) q_selected = argmax(prod(P(Xtilde | q) for Xtilde in d))

            % Working backwards:
            % for each Xtilde in D
            %   multiply up P(Xtilde | q)
            %       = sum((6) * (8) for each integer k where Xtilde - k*q is in range [-6, 6]
            
            % BUT
            % individually the scores for each q get really small
            % because the laplace distribution can return really tiny ~1E-5 value, and mutiplying
            % lots of those together gets *tiny*
            % Instead, use logarithms
            % sum up the log(p_of_xtilde), and because we're adding up length(d) values, divide
            % by length(d)
            
            log_scores = zeros(1, length(qs_poss));
            for Xtilde = d
                for i_q = 1:length(qs_poss)
                    q_poss = qs_poss(i_q);
                    
                    % Neelamani2006 reference Fan2003 for their definition of gamma.
                    % They take gamma=6, but we can just take the value directly from Fan2003:
                    % B(m,n) = D(m)D(n), where D(0) = 2
                    % => their gamma = B(0,0) = 4
                    % Note that this is only for one level of compression - I make the bounds wider
                    % to account for that.
                    gamma = 4;
                    % Xtilde - k*q is in range [-gamma, gamma]
                    % Xtilde - k_max*q = -gamma, k_max = (Xtilde + gamma)/q
                    % k_min = (Xtilde - gamma)/q
                    % Round both towards zero to get integer ks
                    ks = ceil((Xtilde - gamma)/q_poss):floor((Xtilde + gamma)/q_poss);
                    p_of_xtilde = 0;
                    for k = ks
                        expected_err = Xtilde - (k*q_poss);
                        p_roundoff = roundoff_gauss_p_eq_6(gamma,expected_err);
                        p_laplace = integrate_laplacian_eq_8(lambda, k, q_poss);% / (lambda/2);
%                         if p_roundoff <= 0 || p_laplace <= 0 || (p_roundoff * p_laplace) <= 0
%                             p_roundoff, p_laplace, Xtilde, lambda, k, q_poss
%                         end
                        p_of_xtilde = p_of_xtilde + (p_roundoff * p_laplace);
                    end
%                     if p_of_xtilde == 0
%                         Xtilde, lambda, q_poss, ks
%                     end

                    if p_of_xtilde == 0
                        % Due to potential extra rounding factors etc. it may be that we encounter a
                        % number out of range
                        p_of_xtilde = 1E-20;
                    end

                    log_scores(i_q) = log_scores(i_q) + log(p_of_xtilde)/length(d);
                end
            end
            
            qs_poss
            log_scores 
            % OK, all q_poss' have been scored
            % Get the index of the maximum score (= the index of the q with the maximum score)
            [s_max, i_q_max] = max(log_scores);
            qs = [qs qs_poss(i_q_max)];
        end
    end
end
function p=roundoff_gauss_p_eq_6(gamma, expected_err)
    if abs(expected_err) > gamma
        p = 0;
    else
        norm_constant = 1;
        % This is also taken from Fan2003
        sigma_2 = 3;
        p = norm_constant * exp(-expected_err*expected_err/(2*sigma_2));
    end
end
function p=integrate_laplacian_eq_8(lambda,k,q)
    % This function is the equivalent of Neelamani2006 eq (8):
    % integral of (lambda/2 * exp(-lamda * abs(tau))) for tau between (k-0.5)q and (k+0.5)q.
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
% alldivisors taken from 
% https://uk.mathworks.com/matlabcentral/answers/21542-find-divisors-for-a-given-number#answer_419132
function divs = alldivisors(N)
  % compute the set of all integer divisors of the positive integer N
  
  % first, get the list of prime factors of N. 
  facs = factor(N);
  
  divs = [1,facs(1)];
  for fi = facs(2:end)
    % if N is prime, then facs had only one element,
    % and this loop will not execute at all, In that case
    % The set of all divisors is simply 1 and N.
    
    % this outer product will generate all combinations of
    % the divisors found so far, combined with the current
    % divisor fi.
    divs = [1;fi]*divs;
    
    % unique eliminates the replicate divisors, making
    % this an efficient code.
    divs = unique(divs(:)');
  end
end
