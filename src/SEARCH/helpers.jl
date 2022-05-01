findparam_int(ar::activation_record{T,S}) where {T,S} = T;
findparam_real(ar::activation_record{T,S}) where {T,S} = S;

findparam_pfms(::motifs{T,S}) where {T, S} = Vector{Matrix{S}};
findparam_lens(::motifs{T,S}) where {T, S} = Vector{T};
findparams_tbrm_key(::motifs{T,S}) where {T, S} = T;

minimum_filter_length(ar::activation_record) = minimum(ar.filter_lengths);
maximum_filter_length(ar::activation_record) = maximum(ar.filter_lengths);
maximum_apart_length(ar::activation_record) = ar.C - minimum_filter_length(ar) - 1; # C-|fₛ|-1

pfm2pwm(pfm) = log2.(pfm ./ .25);
pfm2pwm(pfm, bg) = log2.(pfm ./ bg);

# pfm_complement(pfm::Matrix{T}) where T <: Real = reduce(hcat,[reverse(pfm[:,i]) for i = size(pfm,2):-1:1]);

pfm_complement(pfm::Matrix{T}) where T <: Real = reverse(pfm);
pwm_complement(pwm::Matrix{T}) where T <: Real = reverse(pwm);

const dummy_complementF64 = Dict(Float64.([1,0,0,0])=>Float64.([0,0,0,1]),  
                        Float64.([0,1,0,0])=>Float64.([0,0,1,0]),
                        Float64.([0,0,1,0])=>Float64.([0,1,0,0]),
                        Float64.([0,0,0,1])=>Float64.([1,0,0,0]));                        

const dummy_complementF32 = Dict(Float32.([1,0,0,0])=>Float32.([0,0,0,1]),  
                        Float32.([0,1,0,0])=>Float32.([0,0,1,0]),
                        Float32.([0,0,1,0])=>Float32.([0,1,0,0]),
                        Float32.([0,0,0,1])=>Float32.([1,0,0,0]));

const dummy_complementF16 = Dict(Float16.([1,0,0,0])=>Float16.([0,0,0,1]),  
                           Float16.([0,1,0,0])=>Float16.([0,0,1,0]),
                           Float16.([0,0,1,0])=>Float16.([0,1,0,0]),
                           Float16.([0,0,0,1])=>Float16.([1,0,0,0]));                        

function dummy_comp(vec::Vector{T}) where T <: Real    
    if T == Float16
        return dummy_complementF64[vec]
    elseif T == Float32
        return dummy_complementF32[vec]
    else
        return dummy_complementF16[vec]
    end
end

dummy_complement_seq(v) = reduce(vcat, [dummy_complement[v[i-4+1:i]] for i = length(v):-4:1]);
# reverse complement of the entire data-set
dat_complement(data_dat_cpu) = reduce(hcat, [dummy_complement_seq(data_dat_cpu[:,n]) for n = 1:size(data_dat_cpu,2)]);

function get_pfm(msa::Matrix{T}, pseudo_count::T=1) where T <: Real 
    pfm_vec = (sum(msa, dims=2) .+ pseudo_count) ./ (size(msa,2) .+ 4*pseudo_count);
    pfm = reshape(pfm_vec , (4, Int(size(msa,1)/4)));
    @assert all(sum(pfm, dims=1) .≈ 1.0) "pfm column need to sum to 1 "
    return pfm
end

function submat_comlement(data_matrix, start_, end_, k, len, len4)
    piece = data_matrix[start_:end_,k];
    reshape(reverse(reshape(piece,(4,len))),len4)
end

function posdicts2pfms(positions::Vector{Dict{T,T}}, lens::Vector{T}, 
                       use_complement::Union{Nothing, Vector{Dict{T,Bool}}}, 
                       num_motifs::T, data_matrix::Matrix{S}, pseudocount_r,
                       pseudo_count::S; smoothing=true, less_pseudocount=false) where {T<:Integer,S<:Real}
    no_comp = isnothing(use_complement);
    pfms = Vector{Matrix{S}}();    
    # preallocate temp array for MSA calculations
    max_len4 = maximum(4 .* lens);
    pos_lens = length.(positions)
    pseudo_count_arr = smoothing ? pseudocount_r .* pos_lens : fill(less_pseudocount ? S(0.1) : pseudo_count, num_motifs);
    # pseudo_count_arr = smoothing ? fill(S(20),num_motifs) : fill(pseudo_count, num_motifs);
    msa_max_cols = maximum(pos_lens);
    msa = zeros(S, (max_len4, msa_max_cols));

    @inbounds for i in 1:num_motifs
        len4 = 4*lens[i];                
        for (j,k) in enumerate(keys(positions[i]))     
            start_ = (positions[i][k]-1)*4+1;
            end_ = start_+len4-1;
            # println("i: $i, key: $k")
            if no_comp
                msa[1:len4,j] = data_matrix[start_:end_,k];   
            else
                if use_complement[i][k]                    
                    msa[1:len4,j] = submat_comlement(data_matrix,start_,end_,k,lens[i],len4);
                else
                    msa[1:len4,j] = data_matrix[start_:end_,k];   
                end
            end                  
        end
        this_msa = @view msa[1:len4,1:pos_lens[i]];
        pfm_vec = (sum(this_msa, dims=2) .+ pseudo_count_arr[i]) ./ (size(this_msa,2) .+ 4*pseudo_count_arr[i]);
        pfm = reshape(pfm_vec , (4, Int(size(this_msa,1)/4)));
        @assert all(sum(pfm, dims=1) .≈ 1.0) "pfm column need to sum to 1 "
        push!(pfms, pfm);
    end
    return pfms
end

get_pwms(pfms, bg) = [pfm2pwm(pfm, bg) for pfm in pfms];
get_thresh(pwms, pval_thresh, type_real) = [type_real(pvalue2score(pwm, pval_thresh)) for pwm in pwms];

function get_msa_init(pos_dict_singleton::Dict{T, T}, len_::Integer, 
                      data_matrix::Matrix{S}; complement=false) where {T <: Integer, S <: Real}
    len4 = 4*len_;
    msa = zeros(S, (len4, length(pos_dict_singleton)));
    @inbounds for (ind,k) in enumerate(keys(pos_dict_singleton))      
        start_ = (pos_dict_singleton[k]-1)*4+1;
        end_ = start_+len_*4-1;
        if complement             
            msa[:,ind] = submat_comlement(data_matrix,start_,end_,k,len_,len4);
        else
            msa[:,ind] = data_matrix[start_:end_,k];
        end
    end
    msa
end

# function get_msa_init(msa::Union{Nothing,Matrix{S}}, 
#                       pos_dict_singleton::Dict{T, T}, 
#                       len_::Integer, data_matrix::Matrix{S}
#                       ) where {S, T}
#     msa = isnothing(msa) ? zeros(S, (4*len_,0)) : msa; 
#     @assert size(msa, 1) == 4*len_ "size of msa matrix must agree"
#     @inbounds for k in keys(pos_dict_singleton)      
#         s = (pos_dict_singleton[k]-1)*4+1;
#         end_ = s+len_*4-1;
#         msa = hcat(msa, data_matrix[s:end_,k]);
#     end
#     return msa
# end


# function get_pfm_colwise(msa::Matrix{T}, pseudo_count=1)  where T <: Real
#     pfm_vec = (sum(msa, dims=2) .+ pseudo_count);
#     count_mat = reshape(pfm_vec , (4, Int(size(msa,1)/4)));
#     pfm = count_mat ./ sum(count_mat, dims=1);
#     @assert all(sum(pfm, dims=1) .≈ 1.0) "pfm column need to sum to 1 "
#     return pfm
# end

function get_pfm_colwise(msa::Matrix{T}, pseudo_count_r)  where T <: Real
    pfm_vec = sum(msa, dims=2);
    count_mat = reshape(pfm_vec , (4, Int(size(msa,1)/4)));
    sum_col = sum(count_mat, dims=1);
    sum_col_ps = pseudo_count_r .* sum_col;
    pfm = (count_mat .+ sum_col_ps) ./ (sum_col+4*sum_col_ps);
    @assert all(sum(pfm, dims=1) .≈ 1.0) "pfm column need to sum to 1 "
    return pfm
end

function get_ic_col_single(pfm_col::Vector{T}, bg::Vector{T}) where T <: Real
    s = T(0);
    @inbounds for i = 1:4
        s += pfm_col[i] * log2(pfm_col[i] / bg[i])
    end
    s
end

function get_ic_col_multiple(pfm::Matrix{T}, bg::Vector{T}) where T <: Real
    s = pfm .* log2.(pfm ./ bg);
    sum(s, dims=1);
end

function get_ic_col(pfm::Matrix{T}) where T <: Real
    ic_col = Float16.(dropdims(2 .+  sum(pfm .* (log2.(pfm)), dims=1), dims=1));
    return ic_col
end

function get_information_content(pfm::Matrix{T}) where T <: Real
    # assume flat background
    ic_col = get_ic_col(pfm);
    return sum(ic_col), ic_col
end

# function get_complement(ms::motifs)
#     new_pfms = [pfm for pfm in ms.pfms];    
#     new_pfms_c = [pfm_complement(pfm) for pfm in ms.pfms];
#     orig_lens = ms.orig_len;
#     append!(orig_lens, ms.orig_len);
#     append!(new_pfms, new_pfms_c);
#     motifs([pfm for pfm in new_pfms], 
#             [size(pfm,2) for pfm in new_pfms], 
#             length(new_pfms), 
#             nothing,
#             nothing,
#             Dict{findparams_tbrm_key(ms),Bool}(k=>false for k = 1:length(new_pfms)),
#             orig_lens);
# end

# function copy_motifs(ms::motifs)
#     motifs([pfm for pfm in ms.pfms], 
#         [size(pfm,2) for pfm in ms.pfms], 
#         length(ms.pfms), 
#         nothing,
#         nothing,
#         Dict{Int16,Bool}(k=>false for k = 1:length(ms.pfms)),
#         ms.orig_len); 
# end

# function rid_of_empty_motifs!(ms::motifs)    
#     len_pfms = length(ms.pfms);
#     keep = [i for i = 1:len_pfms if !isempty(ms.positions[i])]
#     return motifs([ms.pfms[i] for i in keep], 
#                    [ms.lens[i] for i in keep], 
#                    length(keep), 
#                    [ms.positions[i] for i in keep],
#                    [ms.scores[i] for i in keep],
#                    Dict{findparams_tbrm_key(ms),Bool}(k=>false for k = 1:length(keep)),
#                    [ms.orig_len[i] for i in keep]); 
# end

function scores_reevaluation(new_pfms, 
                             new_lens, 
                             new_positions, 
                             data_dat_cpu,
                             bg,
                             type_int, 
                             type_real
                             )

    new_scores = [Dict{type_int,type_real}() for _ = 1:length(new_pfms)];
    @inbounds for l = 1:length(new_pfms)
        pwm = pfm2pwm(new_pfms[l], bg);
        for (k,v) in new_positions[l]
            start_ = (v-1)*4+1; 
            len_ = 4*new_lens[l]; 
            end_ = start_+len_-1;
            view_dat = @view data_dat_cpu[start_:end_,k];
            new_scores[l][k] = sum(pwm .* reshape(view_dat, (4,new_lens[l])));
        end
    end
    return new_scores
end

function re_evaluations!(g::good_stuff{T,S,Q}, 
                        re_evaluate_pfm, 
                        re_evaluate_pwm,
                        re_evaluate_thresh,
                        re_evaluate_scores;
                        smoothing=true,
                        less_pseudocount=false
                        ) where {T,S,Q}        

    if re_evaluate_pfm    
        g.ms.pfms = posdicts2pfms(g.ms.positions , g.ms.lens, 
                                                  g.ms.use_comp, 
                                                  g.ms.num_motifs, 
                                                  g.data.data_matrix, 
                                                  g.search.pseudocounts_r,
                                                  g.search.pseudocounts;
                                                  smoothing=smoothing,
                                                  less_pseudocount=less_pseudocount);
    end                        
                     
    if re_evaluate_pwm
        # @assert re_evaluate_pfm "Must set re-evaluate-pfm to true"
        g.ms.pwms = get_pwms(g.ms.pfms, g.search.bg);
    end                                                    
    if re_evaluate_thresh
        # @assert re_evaluate_pwm "Must set re-evaluate-pwm to true"
        g.ms.thresh = get_thresh(g.ms.pwms, g.search.pval_thresh, S);
    end
    if re_evaluate_scores
        # @assert re_evaluate_pwm "Must set re-evaluate-pwm to true"
        g.ms.scores = scores_reevaluation(g.ms.pfms, g.ms.lens, g.ms.positions, 
                                    g.data.data_matrix, g.search.bg, T, S);          
    end
end