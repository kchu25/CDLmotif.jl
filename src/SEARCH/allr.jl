function allr(p, q, p_count, q_count, bg)
    allr_score = Float64[];
    for i = 1:size(p,2)
        view_p_col = Base.view(p, :, i);
        view_q_col = Base.view(q, :, i);
        nb_p = p_count .* view_p_col;
        nb_q = q_count .* view_q_col;
        a1=sum(nb_p .* pfm2pwm(view_q_col, bg)); 
        a2=sum(nb_q .* pfm2pwm(view_p_col, bg));
        push!(allr_score, (a1+a2)/(sum(nb_p)+sum(nb_q)))
    end
    # println(allr_score)
    # sort!(allr_score, rev=true);
    return sum(allr_score)
end

allr_thresh(pfm_len::Integer, fraction_of_columns, score_each_col) = floor(fraction_of_columns*pfm_len)*score_each_col;

function convolve_allr(pfm_c2, pfm,
                       counts_pfm_c2::Integer,
                       counts_pfm::Integer,                        
                       len_pfm_c2::Integer,
                       len_pfm::Integer,
                       thresh_allr::Real,
                       min_col::Integer,
                       bg                   
                       )
    #= len_pfm_c2 will always be smaller since we've select the ones
        with minimal length
    =#
    # n_cols = Int(floor(perc_col*len_pfm_c2));
    allrs = Float64[];
    # start and end indices for pwms
    # s1e2 for pfm_c2, s2e2 for pfm
    s1e1s = UnitRange{Int}[];
    s2e2s = UnitRange{Int}[];
    l_dec_1 = Int[]; l_dec_2 = Int[]; 
    r_inc_1 = Int[]; r_inc_2 = Int[];

    for i = 1:(len_pfm_c2+len_pfm-1)
        s1 = max(1, len_pfm_c2-i+1); e1 = min(len_pfm_c2, len_pfm_c2-(i-len_pfm));
        s2 = max(1, i-len_pfm_c2+1); e2 = min(i, len_pfm);
        overlap_count = e1-s1+1;
        push!(s1e1s, s1:e1); push!(s2e2s, s2:e2);
        #=  
            Note that:
            1) no need to calculate if the number of columns of the 
            pfm is less than min_col as specified
            2) no need to calculate the placements for which
            maximal value of the score is below the threshold
        =#
        if overlap_count ≥ min_col && 2*overlap_count ≥ thresh_allr 
            push!(allrs, allr(Base.view(pfm_c2,:,s1:e1), Base.view(pfm,:,s2:e2), 
                            counts_pfm_c2, counts_pfm, bg));
        else
            push!(allrs, -Inf);
        end        
        push!(l_dec_1, max(s2-1,0)); push!(l_dec_2, max(s1-1,0));
        push!(r_inc_1, max(0,len_pfm-i)); push!(r_inc_2, max(i-e2,0));
    end
    argmax_ind = argmax(allrs);
    return allrs[argmax_ind], 
           l_dec_1[argmax_ind], 
           r_inc_1[argmax_ind], 
           l_dec_2[argmax_ind], 
           r_inc_2[argmax_ind]
end

function longest_len_w_most_cover!(ms, msa_counts, grouped)
    indices = (1:ms.num_motifs )[.!grouped];
    coverings = msa_counts[.!grouped];
    sizes = ms.lens[.!grouped];
    max_size = maximum(sizes);
    mask = sizes .== max_size;
    chosen_ind = indices[mask][argmax(coverings[mask])];
    return chosen_ind
end

function most_freq_occured_len!(ms, msa_counts, grouped)
    indices = (1:ms.num_motifs)[.!grouped];
    coverings = msa_counts[.!grouped];
    sizes = ms.lens[.!grouped];
    mode_size = mode(sizes);
    mask = sizes .== mode_size;
    chosen_ind = indices[mask][argmax(coverings[mask])];
    return chosen_ind
end


function allr_merge!(g::good_stuff{T,S,Q};
                     re_evaluate_pfm=false,
                     re_evaluate_pwm=true,
                     re_evaluate_thresh=true,
                     re_evaluate_scores=false
                     ) where {T,S,Q}
    new_pfms = Matrix{S}[];
    msa_counts = length.(g.ms.positions);
    allr_threshs = allr_thresh.(g.ms.lens, 
                                g.search.allr_fraction_of_columns, 
                                g.search.allr_score_each_column);

    #= minimum columns for each PWM that's required for the allr scores to be calculated
        when it's being selected to be compare to the root (chosen one)
    =#
    min_cols = Int.(floor.(g.search.allr_fraction_of_columns .* g.ms.lens));
    # min_cols = T.(g.ms.lens .- 2);

    # root and grouped information
    root_indices = Int[];
    grouped = fill(false, g.ms.num_motifs);

    while !all(grouped)
        # choose a root
        root_ind = longest_len_w_most_cover!(g.ms, msa_counts, grouped); 
        # root_ind = most_freq_occured_len!(g.ms, msa_counts, grouped); 
        
        grouped[root_ind] = true;
       
        good_matches = Int[]; 
        ld1_matches = Int[]; ri1_matches = Int[];
        ld2_matches = Int[]; ri2_matches = Int[];
        ld1s = Int[]; indices_about_to_be_grouped = Int[];

        good_matches_c = Int[]; 
        ld1_matches_c = Int[]; ri1_matches_c = Int[];
        ld2_matches_c = Int[]; ri2_matches_c = Int[];
        ld1s_c = Int[]; indices_about_to_be_grouped_c = Int[];

        @inbounds for i = 1:g.ms.num_motifs            
            if i != root_ind && g.ms.lens[root_ind] ≥ g.ms.lens[i] && !grouped[i]
                allr_score, ld1, ri1, ld2, ri2 = convolve_allr(
                                                    g.ms.pfms[root_ind],
                                                    g.ms.pfms[i],
                                                    msa_counts[root_ind],
                                                    msa_counts[i],
                                                    g.ms.lens[root_ind],
                                                    g.ms.lens[i],
                                                    allr_threshs[i],
                                                    min_cols[i],
                                                    g.search.bg);
                allr_score_c, ld1_c, ri1_c, ld2_c, ri2_c = convolve_allr(
                                                    g.ms.pfms[root_ind],
                                                    pfm_complement(g.ms.pfms[i]),
                                                    msa_counts[root_ind],
                                                    msa_counts[i],
                                                    g.ms.lens[root_ind],
                                                    g.ms.lens[i],
                                                    allr_threshs[i],
                                                    min_cols[i],
                                                    g.search.bg); 
                                                                                                                       
                if allr_score ≥ allr_score_c && allr_score > allr_threshs[i]                     
                    push!(indices_about_to_be_grouped, i);              
                    push!(good_matches, i);
                    push!(ld1_matches, ld1);
                    push!(ld2_matches, ld2);
                    push!(ld1s, ld1+ld2);
                    push!(ri1_matches, ri1);
                    push!(ri2_matches, ri2);       
                elseif allr_score_c ≥ allr_score && allr_score_c > allr_threshs[i]   
                    push!(indices_about_to_be_grouped_c, i);              
                    push!(good_matches_c, i);
                    push!(ld1_matches_c, ld1_c);
                    push!(ld2_matches_c, ld2_c);
                    push!(ld1s_c, ld1_c+ld2_c);
                    push!(ri1_matches_c, ri1_c);
                    push!(ri2_matches_c, ri2_c); 
                end                
            end
        end

        @inbounds for i in indices_about_to_be_grouped grouped[i] = true; end
        @inbounds for i in indices_about_to_be_grouped_c grouped[i] = true; end
        
        push!(root_indices, root_ind);     

        @inbounds if length(good_matches) == 0
            push!(new_pfms, g.ms.pfms[root_ind]);            
        else                   
            merged_len = 4*g.ms.lens[root_ind];
            msa = zeros(S, (merged_len,0));

            for (ind,good_match) in enumerate(good_matches)
                left = 4*(ld2_matches[ind])+1;
                right = 4*(g.ms.lens[root_ind]-ri2_matches[ind]-1)+4;

                pos_dict_match = Dict{T,T}(
                    k=>(g.ms.positions[good_match][k]+ld1_matches[ind])
                    for k in keys(g.ms.positions[good_match]));
           
                msa_match = get_msa_init(pos_dict_match, 
                            g.ms.lens[good_match]-ri1_matches[ind]-ld1_matches[ind], 
                            g.data.data_matrix);         

                mas_match_expand = zeros(S, (merged_len, length(g.ms.positions[good_match])));
                mas_match_expand[left:right,:] = msa_match;
                msa = hcat(msa, mas_match_expand);
            end

            for (ind,good_match) in enumerate(good_matches_c)
                left = 4*(ld2_matches_c[ind])+1;
                right = 4*(g.ms.lens[root_ind]-ri2_matches_c[ind]-1)+4;

                pos_dict_match = Dict{T,T}(
                    k=>(g.ms.positions[good_match][k]+ri1_matches_c[ind])
                    for k in keys(g.ms.positions[good_match]));
           
                msa_match = get_msa_init(pos_dict_match, 
                            g.ms.lens[good_match]-ld1_matches_c[ind]-ri1_matches_c[ind], 
                            g.data.data_matrix; complement=true);         
                                      
                mas_match_expand = zeros(S, (merged_len, length(g.ms.positions[good_match])));
                mas_match_expand[left:right,:] = msa_match;
                msa = hcat(msa, mas_match_expand);
            end

            msa_root = get_msa_init(g.ms.positions[root_ind], 
                                    g.ms.lens[root_ind], 
                                    g.data.data_matrix);             
            msa = hcat(msa, msa_root);
            
            # ps=S(g.sea*size(msa_root,2))
            # ps=g.search.pseudocounts;

            merged_pfm = get_pfm_colwise(msa, g.search.pseudocounts_r);
            # merged_pfm = get_pfm_colwise(msa, g.search.pseudocounts);
            # merged_pfm = get_pfm_colwise(msa);

            if !isempty(indices_about_to_be_grouped)
                min_lds = minimum(ld2_matches);
                min_ris = minimum(ri2_matches);
                l = size(merged_pfm,2);
                ### print ####                
                # q=(l-min_ris)-(min_lds+1)+1;
                # println("truncated! o:$(ms.lens[root_ind]) t:$q")
                ##############
                push!(new_pfms, merged_pfm[:,(min_lds+1):(l-min_ris)]);
            else
                push!(new_pfms, merged_pfm);
            end
        end
    end    

    g.ms.pfms = new_pfms;
    g.ms.lens = size.(new_pfms,2);
    g.ms.num_motifs = length(new_pfms);
    # g.ms.orig_len = g.ms.orig_len[root_indices];
    g.ms.positions = nothing;
    g.ms.scores = nothing;
    re_evaluations!(g, re_evaluate_pfm, re_evaluate_pwm, re_evaluate_thresh, re_evaluate_scores);          
end
