#=
Best possible score of a PWM

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
# Output
    the best possible score this matrix can get
=#
best_score(pwm::Matrix{T}) where T <: Real = @inbounds sum(maximum(pwm[:,i]) for i = 1:size(pwm,2));
best_score(pwm_col::Vector{T}) where T <: Real = maximum(pwm_col);

#=
Worst possible score of a PWM

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
# Output
    the worst possible score this matrix can get
=#
worst_score(pwm::Matrix{T}) where T <: Real = @inbounds sum(minimum(pwm[:,i]) for i = 1:size(pwm,2));
worst_score(pwm_col::Vector{T}) where T <: Real = minimum(pwm_col);

#=
Return a column-permuted PWM that minimize the score range so that δ₁ ≥ δ₂ ≥ … ≥ δₘ
where δᵢ = best_score(pwm[:,i])-worst_score(pwm[:,i]). 

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
# Output
    a column-permuted pwm 
=#
min_score_range(pwm::Matrix{T}) where T <: Real = @inbounds pwm[:,sortperm([best_score(pwm[:,i])-worst_score(pwm[:,i]) 
                                                    for i = 1:size(pwm,2)],rev=true)];

#=
"Round the PWM"

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
    `granularity`: a small positive real number e.g. 0.01 or 0.001, etc.
# Output
    a rounded pwm of the input pwm
=#
round_pwm(pwm::Matrix{T}, granularity::Real) where T <: Real = floor.(pwm ./ granularity) * granularity;

#=
The maximum error induced by the rounded pwm M_ϵ 
(see definition 3 in https://almob.biomedcentral.com/articles/10.1186/1748-7188-2-15; this is the quantity E) 

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
    `granularity`:
# Output
    A positve real number that's the maximum error induced by the rounded pwm M_ϵ 
=#

calc_E(pwm, pwm_rounded) = @inbounds sum(maximum(pwm[:,i]-pwm_rounded[:,i]) for i = size(pwm));

#=
Note: Use a nested dictionary to represent the distribution Q
    Since Q is used to reference the probability of (M[1…i],score),
    the keys in the first layer is i, and the value in the first layer 
    are dictionaries with scores as keys and probability as values

    call create_Q(m) to initialize such a distribution Q
        where m is the "width" of the PWM
=#
create_Q(m::Integer) = Dict{Int16,SortedDict{Double64,Double64}}(i==0 ? i=>SortedDict(0=>1) : i=>SortedDict() for i=0:m);

#=
Input: 
    pwm: a 4 x m matrix
    α, β: score interval [α, β]
    bg: 4 x 1 vector that specifies the multinomial genomic background; default to flat background.    
Output:
    Q: a probability mass table
        e.g. Q[m] shows all the weights of P[pwm_score = η] for α ≤ η ≤ β
=#
function score_distribution(pwm_::Matrix{T}, 
                            alpha::Real, β::Real, 
                            bg=_bg_
                            ) where T <: Real    
    m = size(pwm_,2);
    Q = create_Q(m);
    @inbounds for i = 1:m        
        bs = i+1 > m ? 0 : best_score(pwm_[:,i+1:m]);
        ws = i+1 > m ? 0 : worst_score(pwm_[:,i+1:m]);
        for score in keys(Q[i-1])
            for j = 1:4
                t = score + pwm_[j,i];
                if alpha - bs ≤ t ≤ β - ws
                    if haskey(Q[i], t)
                        Q[i][t] += Q[i-1][score]*bg[j];
                    else
                        Q[i][t] = Q[i-1][score]*bg[j];
                    end
                end
            end
        end
    end
    return Q
end

#=
Return the probability of the background model that has score ≥ α with respect to the input pwm

Input: 
    pwm: a 4 x m matrix
    α: a score threshold
    bg: 4 x 1 vector that specifies the multinomial genomic background; default to flat background.    
Output:
    pval: the probability that the background can have score ≥ α with respect to pwm
=#
function fast_pvalue(pwm::Matrix{T}, 
                     alpha::Real,
                     bg=_bg_
                     ) where T <: Real
    m = size(pwm,2);
    pwm_d64 = Double64.(pwm);
    Q = create_Q(m);
    bg_d64 = Double64.(bg);
    pval = Double64(0.0);
    @inbounds for i = 1:m
        bs = i+1 > m ? 0 : best_score(pwm_d64[:,i+1:m]);
        ws = i+1 > m ? 0 : worst_score(pwm_d64[:,i+1:m]);
        for (score,_) in Q[i-1]
            for j = 1:4
                t = score + pwm_d64[j,i];
                if alpha - ws ≤ t
                    pval = pval + Q[i-1][score]*bg_d64[j];
                elseif alpha - bs ≤ t
                    if haskey(Q[i], t)
                        Q[i][t] += Q[i-1][score]*bg_d64[j];
                    else
                        Q[i][t] = Q[i-1][score]*bg_d64[j];
                    end
                end
            end
        end
    end
    return pval
end

# return the sum of all the weights 
Q_sum(Q_m::SortedDict{Double64,Double64}) = sum(v for v in values(Q_m));

function find_largest_alpha(Q_m::SortedDict{T,T}, pval::T) where T <: Real
    q_sum = Q_sum(Q_m);
    largest_k = nothing;
    for (k,v) in Q_m
        if q_sum ≥ pval
            largest_k = k;
        else
            return k
        end
        q_sum -= v;
    end
    return largest_k
end

function pval_w_Qm(Qm::SortedDict{T,T}, alpha::Real) where T <: Real
    pval = 0;
    for (k,v) in Qm
        if k ≥ alpha
            pval += v;
        end
    end 
    return pval
end

function find_δ(Q_m::SortedDict{T,T}, pval_ϵ::Real, pval::Real) where T <: Real
    q_sum_plus_pval_ϵ = Q_sum(Q_m)+pval_ϵ;
    largest_δ = nothing;
    for (k,v) in Q_m
        if q_sum_plus_pval_ϵ ≥ pval
            largest_δ = k;
        else
            return k
        end
        q_sum_plus_pval_ϵ -= v;
    end
    return largest_δ
end

"""
    pval2score(pwm, pval, ϵ=1e-1, k=10, bg=[.25,.25,.25,.25])
Returns the highest score(M,pval) of a `pwm` such that p-value is greater or equal to `pval`.

Input:
* `pwm`: a 4 x m matrix
* `pval`: a p-value; e.g. pval = 1e-3
* `ϵ`: initial granularity  (optional) 
* `k`: Refinement parameter (optional)
* `bg`: multinomial background (optional)

Output    
* `alpha`: the highest score-threshold
"""
function pvalue2score(pwm::Matrix{T}, 
                      pval::Real, 
                      ϵ=_granularity_, 
                      bg=_bg_
                      ) where T <: Real
    @assert 0 ≤ pval ≤ 1 "pvalue must be in [0,1]"
    @assert size(pwm,1) == 4 "The input matrix must have only 4 rows"
    pwm_double64 = Double64.(pwm);
    pval_double64 = Double64(pval);
    bg_double64 = Double64.(bg);
    mpwm = min_score_range(pwm_double64);
    m = size(pwm, 2);
    pwm_ϵ = round_pwm(mpwm, ϵ);
    # E = calc_E(mpwm, pwm_ϵ);
    Q = create_Q(m);
    Q = score_distribution(pwm_ϵ, worst_score(pwm_ϵ),Inf,bg_double64);
    alpha = find_largest_alpha(Q[m], pval_double64);

    return alpha;
end