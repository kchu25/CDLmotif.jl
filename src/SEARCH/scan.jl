####################################################################################
function motifs_prep(ms::motifs{T,S}) where {T<:Integer,S<:Real}
    positions = [Dict{T,T}() for _ = 1:ms.num_motifs];
    scores = [Dict{T,S}() for _ = 1:ms.num_motifs];
    use_comp = [Dict{T,Bool}() for _ = 1:ms.num_motifs];
    # to_be_rm = Dict{type_int,Bool}(k=>false for k = 1:ms.num_motifs);
    return positions, scores, use_comp
end

function scan_n(data_n, pwm::Matrix{S}, pwm_len::T, L::T) where {S<:Real,T<:Integer}
    # scan the sequence in the orientation as in the dataset 
    max_score_ind = -1; max_score = -Inf;
    for i = 1:L-pwm_len+1 # devectorization
        score = S(0);
        @inbounds for j = 1:4, k = 1:pwm_len
            score += pwm[j,k] * data_n[j,i+k-1];
        end
        score > max_score && (max_score = score; max_score_ind = i);
    end
    return max_score_ind, max_score
end

function overlapping_scan_bg!(g::good_stuff{T,S,Q}) where {T,S,Q}
    p = [Dict{T,T}() for _ = 1:g.ms.num_motifs];
    s = [Dict{T,S}() for _ = 1:g.ms.num_motifs];
    for k = 1:g.ms.num_motifs   
        pwm_comp = pwm_complement(g.ms.pwms[k]);      
        for n = 1:size(g.data.data_matrix_bg,2)
            data_n = reshape(Base.view(g.data.data_matrix_bg,:,n), (4, g.data.L)); 
            argmax_, score = scan_n(data_n,
                                    g.ms.pwms[k], 
                                    g.ms.lens[k], 
                                    g.data.L);
            argmax_c, score_c = scan_n(data_n,
                                    pwm_comp, 
                                    g.ms.lens[k], 
                                    g.data.L);
            if score ≥ score_c
                score > g.ms.thresh[k] && (s[k][n] = score; p[k][n] = argmax_;)
            else
                score_c > g.ms.thresh[k] && (s[k][n] = score_c; p[k][n] = argmax_c;)
            end
        end
    end
    return p, s
end

function non_overlapping_scan_bg!(g::good_stuff{T,S,Q}, positions, scores) where {T,S,Q}
    @inbounds for n = 1:g.data.N
        spans_pos = T[]; 
        spans_len = T[];

        indices_to_keep = Dict(i=>false for i = 1:g.ms.num_motifs);
        scores_n = [haskey(scores[i], n) ? scores[i][n] : 0 for i = 1:g.ms.num_motifs];
        positions_ = [haskey(positions[i], n) ? positions[i][n] : 0 for i = 1:g.ms.num_motifs];   
        max_score_ind = argmax(scores_n);

        while scores_n[max_score_ind] != 0
            intersect_ = false;
            for (p,l) in zip(spans_pos,spans_len)
                # check whether this "segment" that achieves 
                # the best score intersect with any other "segments"
                p_end = p+l-1;
                p_max = positions_[max_score_ind];
                p_max_end = p_max+g.ms.lens[max_score_ind]-1;
                if p ≤ p_max ≤ p_end || p ≤ p_max_end ≤ p_end || p_max ≤ p ≤ p_max_end || p_max ≤ p_end ≤ p_max_end
                    intersect_ = true;
                end
            end
            if !intersect_
                indices_to_keep[max_score_ind] = true;
                push!(spans_pos, positions_[max_score_ind]);
                push!(spans_len, g.ms.lens[max_score_ind]);
            end
            scores_n[max_score_ind] = 0;
            max_score_ind = argmax(scores_n);
        end
        for j = 1:g.ms.num_motifs         
            if !indices_to_keep[j] && haskey(positions[j], n)
                delete!(positions[j], n); 
                delete!(scores[j], n); 
                # so that the set of sequences covered by pwms are disjoint
            end
        end
    end
    return positions, scores
end

function overlapping_scan!(g::good_stuff{T,S,Q},
                           re_evaluate_pfm=false,
                           re_evaluate_pwm=false,
                           re_evaluate_thresh=false,
                           re_evaluate_scores=false) where {T,S,Q}
    p, s, u = motifs_prep(g.ms);

    for k = 1:g.ms.num_motifs   
        pwm_comp = pwm_complement(g.ms.pwms[k]);        
        for n = 1:g.data.N
            data_n = reshape(Base.view(g.data.data_matrix,:,n), (4, g.data.L));                       
            argmax_, score = scan_n(data_n,
                                    g.ms.pwms[k], 
                                    g.ms.lens[k], 
                                    g.data.L);
            argmax_c, score_c = scan_n(data_n,
                                    pwm_comp, 
                                    g.ms.lens[k], 
                                    g.data.L);
            if score ≥ score_c
                (score > g.ms.thresh[k]) && (s[k][n] = score; p[k][n] = argmax_; u[k][n] = false;)
            else
                (score_c > g.ms.thresh[k]) && (s[k][n] = score_c; p[k][n] = argmax_c; u[k][n] = true;)
            end
        end
    end
  
    g.ms.positions = p;
    g.ms.scores = s;
    g.ms.use_comp = u;
    re_evaluations!(g, re_evaluate_pfm, re_evaluate_pwm, 
                        re_evaluate_thresh, re_evaluate_scores);
end

function non_overlap_scan!(g::good_stuff{T,S,Q}; 
                           re_evaluate_pfm=true,
                           re_evaluate_pwm=true,
                           re_evaluate_thresh=true,
                           re_evaluate_scores=false,
                           smoothing=true,
                           less_pseudocount=false
                           ) where {T,S,Q}

    @inbounds for n = 1:g.data.N
        spans_pos = T[]; 
        spans_len = T[];

        indices_to_keep = Dict(i=>false for i = 1:g.ms.num_motifs);
        scores_n = [haskey(g.ms.scores[i], n) ? g.ms.scores[i][n] : 0 for i = 1:g.ms.num_motifs];
        positions = [haskey(g.ms.positions[i], n) ? g.ms.positions[i][n] : 0 for i = 1:g.ms.num_motifs];   
        max_score_ind = argmax(scores_n);

        while scores_n[max_score_ind] != 0
            intersect_ = false;
            for (p,l) in zip(spans_pos,spans_len)
                # check whether this "segment" that achieves 
                # the best score intersect with any other "segments"
                p_end = p+l-1;
                p_max = positions[max_score_ind];
                p_max_end = p_max+g.ms.lens[max_score_ind]-1;
                if p ≤ p_max ≤ p_end || p ≤ p_max_end ≤ p_end || p_max ≤ p ≤ p_max_end || p_max ≤ p_end ≤ p_max_end
                    intersect_ = true;
                end
            end
            if !intersect_
                indices_to_keep[max_score_ind] = true;
                push!(spans_pos, positions[max_score_ind]);
                push!(spans_len, g.ms.lens[max_score_ind]);
            end
            scores_n[max_score_ind] = 0;
            max_score_ind = argmax(scores_n);
        end
        for j = 1:g.ms.num_motifs         
            if !indices_to_keep[j] && haskey(g.ms.positions[j], n)
                delete!(g.ms.positions[j], n); 
                delete!(g.ms.scores[j], n); 
                delete!(g.ms.use_comp[j], n); 
                # so that the set of sequences covered by pwms are disjoint
            end
        end
    end
    re_evaluations!(g, re_evaluate_pfm, re_evaluate_pwm, 
                            re_evaluate_thresh, re_evaluate_scores;
                            smoothing=smoothing,
                            less_pseudocount=less_pseudocount);
end

function scan_w_gpu!(g::good_stuff, data_matrix_gpu; 
                    re_evaluate_pfm=false,
                    re_evaluate_pwm=false,
                    re_evaluate_thresh=false,
                    re_evaluate_scores=false,
                    scan_bg=false
                    )
    # setup 
    type_real = findparam_real(g.ms); type_int = findparam_int(g.ms);
    maxlen = maximum(g.ms.lens); pwms = zeros(type_real, g.ms.num_motifs, 4, maxlen);
    @inbounds for i = 1:g.ms.num_motifs pwms[i,:,1:g.ms.lens[i]] = g.ms.pwms[i]; end

    # execute
    pos_scores = CUDA.zeros(type_real, g.ms.num_motifs, g.data.N, g.data.L);
    @cuda threads=ker_3d blocks=b_size_3d(pos_scores) greedy_search!(cu(pwms), 
                                                                    data_matrix_gpu, 
                                                                    cu(g.ms.lens), 
                                                                    pos_scores,
                                                                    cu(g.ms.thresh));
    argmax_pos_scores = argmax(pos_scores, dims=3);
    argmax_pos_scores_cpu = Array(argmax_pos_scores);
    pos_scores_pos = dropdims(map(x-> x[3]==1 ? 0 : x[3], argmax_pos_scores_cpu),dims=3);
    found = findall(pos_scores_pos .> 0);

    positions = [Dict{type_int, type_int}() for _ = 1:g.ms.num_motifs];
    # use_comp = [Dict{type_int, Bool}() for _ = 1:g.ms.num_motifs];
    for f in found 
        positions[f[1]][f[2]] = pos_scores_pos[f[1],f[2]]; 
        # use_comp = [f[1]][f[2]] = false;
    end

    if scan_bg
        return positions;
    else
        indicator = length.(positions) .> 0;
        g.ms.positions=positions[indicator];
        g.ms.pfms=g.ms.pfms[indicator];
        g.ms.lens=g.ms.lens[indicator];
        g.ms.num_motifs=length(g.ms.pfms);
        # g.ms.use_comp=use_comp;
        re_evaluations!(g, re_evaluate_pfm, re_evaluate_pwm, 
                        re_evaluate_thresh, re_evaluate_scores);
    end

    pwms = nothing; pos_scores = nothing; pos_scores_pos = nothing; argmax_pos_scores = nothing;
end

function greedy_search!(pwms, data_dat_gpu, lens, pos_scores, thresh)
    k = (blockIdx().x - 1) * blockDim().x + threadIdx().x; # kth pair
    n = (blockIdx().y - 1) * blockDim().y + threadIdx().y; # nth sequence
    l = (blockIdx().z - 1) * blockDim().z + threadIdx().z; # lth position

    L, N = size(data_dat_gpu); L_div_4 = CUDA.Int(L/4);
    K, _, _ = size(pwms);
    if k ≤ K && n ≤ N && l ≤ L_div_4-lens[k]+1
        @inbounds for (ind,i) in enumerate(l:l+lens[k]-1)
            for a = 1:4
                pos_scores[k,n,l] += pwms[k,a,ind]*data_dat_gpu[(i-1)*4+a,n];                
            end
        end
        pos_scores[k,n,l] = pos_scores[k,n,l] > thresh[k] ? pos_scores[k,n,l] : 0f0;
    end
    return nothing
end