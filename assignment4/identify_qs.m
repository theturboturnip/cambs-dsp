function qs = identify_qs(type, d, h)
    if type == "dc"
        h = histcounts(abs(d),'BinMethod','Integer');
        peaks_fft = estimate_dc_fft_qs(h);
        % Imagine a dirac comb with period T, sampled between T = [0, 1].
        % The FFT of that sample will have a period F = 1/T
        % *If* the histogram values were sampled between T = {0, 1}, 
        % then T = 1/F
        % but they're not, the histogram values are over the range X = {1,
        % length(h)}
        % => T = length(h)
        % e.g. length(h)/qs_fft_est
        period = (length(h)./peaks_fft);
        period = round(period);
        if ~isempty(period)
            period = unique(period);
        end

        if length(period) == 2
            % Determine the order of quantizations
            period = [min(period) max(period)];

            % If period(2) is a factor of period(1)
            ratio = period(2) / period(1);
            if abs(ratio - round(ratio)) < 0.1
                % Then period(2) is likely the first point where the first and second
                % quantization levels overlapped, rather than being the
                % actual second quantization level.
                % e.g. for period = [10, 40]
                % Given that there were values at each multiple of 10
                % (because we picked it up in the first place) and not just
                % multiples of 40, this data can't have been quantized with
                % q=40.
                % Instead, it's much more likely that it was quantized with
                % q=8 and q=10. This would leave small peaks at
                % e.g. 10, 20, 30 (which are all close to multiples of 8,
                % but not exact) and large peaks at 40 (which is a multiple
                % of both 10 and 8).

                % Now given [10, 40], we want to find 8.
                % This has to be a number that's
                % a) a divisor of 40
                % b) not a multiple of 10
                % c) not a divisor of 10
                % The numbers that fulfil these conditions for [10, 40] are
                % 4 and 8.
                % If 4 were the correct value, then there would be a larger
                % spike at 20, which would have been picked up instead of
                % 40.
                % => also remove divisors of 10*2, 10*3
                divs = divisors(period(2)); % = divisors of 40
                divs = setdiff(divs, divisors(period(1))); % remove divisors of 10
                divs = setdiff(divs, (1:ratio) .* period(1)); % remove multiples of 10
                % For all multiples of 10, up to and not including 40,
                % remove *those* divisors.
                for mult = (1:(ratio-1)) .* period(1)
                    divs = setdiff(divs, divisors(mult));
                end

                if length(divs) == 1
                    qs = [divs(1) period(1)];
                elseif length(divs) > 1
                    qs = [choose_best_q(divs, d) period(1)];
                else
                    error("period %d is a multiple of period %d, but couldn't find a suitable quantization", period(2), period(1));
                end
          
            else
                qs = period;
            end
        else
            qs = period;
        end


        % Refinement doesn't help here - the Neelamani approach 
        % assumes the 0th peak is highest, and the rest get continually
        % lower as modelled by a laplacian distribution
    elseif type == "ac"
        h = histcounts(abs(d),'BinMethod','Integer');
        qs_est = estimate_ac_histogram_qs(h)
        REFINE_ESTIMATE = 1;
        
        if ~REFINE_ESTIMATE || isempty(qs_est)
            qs = qs_est;
        else
            qs = refine_qs(qs_est, d);
        end
    end
end
function qs=refine_qs(qs_est, d)
    qs = [];
    % Refine each q_est in qs_est
    for q_est = qs_est
        % As in [Fan2003], choose possible parameter space from q_est
        % Q, Q+1, Q-1.
        % We don't consider integer factors here currently.
        % the selection process is unfairly biased towards smaller factors, 
        qs_poss = [q_est, q_est-1, q_est+1];
        if q_est == 1
            % Don't try a quantization of 0
            qs_poss = [1 2];
        end

        qs = [qs choose_best_q(qs_poss, d)];
    end
end
function q=choose_best_q(qs_poss, d)
    % Take Neelamani2006 approach to refining estimated Qs
    % They treat d as the 'set D of the observed coefficients'
    % This would seemingly imply to remove duplicate elements,
    % but I don't think that's correct.
    % More common coefficient values should have a larger effect on the score.
    N = length(d);
    lambda = N/sum(abs(d));

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
            gamma = 6;
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
                p_of_xtilde = 1E-99;
            end

            log_scores(i_q) = log_scores(i_q) + log(p_of_xtilde)/length(d);
        end
    end
    
    % OK, all q_poss' have been scored
    % Get the index of the maximum score (= the index of the q with the maximum score)
    [s_max, i_q_max] = max(log_scores);
    q = qs_poss(i_q_max);

    qs_poss, q
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
