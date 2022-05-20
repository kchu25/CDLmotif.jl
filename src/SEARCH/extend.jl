function explore_sides!(offset, 
                        L,
                        type_real, 
                        positions, 
                        included, 
                        num_included, 
                        data_matrix,
                        bg, 
                        mlen,
                        pseudocount;
                        right=false)
    this_col = zeros(type_real, 4); 
    inc_right = mlen-1+offset;
    to_exclude = Int[];
    for (ind,k) in enumerate(keys(positions))
        if included[ind]
            pos = right ? positions[k] + inc_right : positions[k]-offset;
            # pos > L && println("Here! key $k")
            (pos < 1 || pos > L) && (push!(to_exclude,ind); num_included-=1; continue;)
            s_ = (pos-1)*4+1; e_ = s_+4-1;
            dat = @view data_matrix[s_:e_,k];
            this_col = this_col + dat;                
        end
    end   
    this_col = (this_col .+ pseudocount) ./ (sum(this_col)+ 4*pseudocount);
    ic_col = get_ic_col_single(this_col, bg);

    return ic_col, to_exclude
end

function extension!(L, positions, mlen,
                    data_matrix, bg, type_int, type_real, pseudocount; 
                    perc_no_less_than=0.95, 
                    ic_extension_thresh=0.075
    )
    len_pos = length(positions);
    dec = 0; inc = 0; included = fill(true, len_pos); 
    num_included = len_pos; usage_lowerbdd = Int(floor(perc_no_less_than*len_pos));
    
    # controls
    explore_left = true; explore_right = true;
    ic_col_left = nothing; ic_col_right = nothing;
    
    # results
    new_positions = positions; len_inc = 0;

    while true 
        # left
        if explore_left
            ic_col_left, to_exclude = explore_sides!(dec+1, L, type_real, positions, 
                                      included, num_included, 
                                      data_matrix, bg, mlen, 
                                      pseudocount; right=false);
            if ic_col_left ≤ ic_extension_thresh 
                explore_left = false 
            else 
                to_exclude!(included, to_exclude); len_inc += 1; dec+=1;                
            end
        end
        # right
        if explore_right
            ic_col_right, to_exclude = explore_sides!(inc+1, L, type_real, positions, 
                                      included, num_included, 
                                      data_matrix, bg, mlen, 
                                      pseudocount; right=true);
        
            if ic_col_right ≤ ic_extension_thresh
                explore_right = false 
            else 
                to_exclude!(included, to_exclude); len_inc += 1; inc+=1;                
            end
        end
        (!explore_left && !explore_right) && break;
        (num_included ≤ usage_lowerbdd) && break;
    end    
    
    (dec > 0 || inc > 0) && (new_positions = Dict{type_int,type_int}(k=>positions[k]-dec 
        for (ind,k) in enumerate(keys(positions)) if included[ind]);)
            
    return new_positions, mlen+len_inc
end

function explore_sides!(offset, 
                        L,
                        type_real, 
                        positions, 
                        comp,
                        included, 
                        num_included, 
                        data_matrix,
                        bg, 
                        mlen,
                        pseudocount;
                        right=false)
    this_col = zeros(type_real, 4); 
    inc_right = mlen-1+offset;
    to_exclude = Int[];
    for (ind,k) in enumerate(keys(positions))
        if included[ind]
            pos = nothing;            
            if (right && !comp[k]) || (!right && comp[k])
                pos = positions[k]+inc_right; 
            else # (right && comp[k]) || (!right && !comp[k])
                pos = positions[k]-offset;
            end            
            (pos < 1 || pos > L) && (push!(to_exclude,ind); num_included-=1; continue;)
            s_ = (pos-1)*4+1; e_ = s_+4-1;
            dat = data_matrix[s_:e_,k];
            dat = comp[k] ? dummy_comp(dat) : dat;
            this_col = this_col + dat;                
        end
    end   
    this_col = (this_col .+ pseudocount) ./ (sum(this_col)+ 4*pseudocount);
    ic_col = get_ic_col_single(this_col, bg);
    return ic_col, to_exclude
end

function to_exclude!(included::Vector{Bool}, to_exclude)
    for i in to_exclude included[i] = false; end
end

function extension!(L, positions, mlen, comp,
                    data_matrix, bg, type_int, type_real, pseudocount; 
                    perc_no_less_than=0.95, 
                    ic_extension_thresh=0.075
    )
    len_pos = length(positions);
    dec = 0; inc = 0; included = fill(true, len_pos); 
    num_included = len_pos; usage_lowerbdd = Int(floor(perc_no_less_than*len_pos));
    
    # controls
    explore_left = true; explore_right = true;
    ic_col_left = nothing; ic_col_right = nothing;
    
    # results
    new_positions = positions; len_inc = 0;

    while true 
        # left
        if explore_left
            ic_col_left, to_exclude = explore_sides!(dec+1, L, type_real, positions, comp,
                                      included, num_included, 
                                      data_matrix, bg, mlen, 
                                      pseudocount; right=false);
            if ic_col_left ≤ ic_extension_thresh 
                explore_left = false 
            else 
                to_exclude!(included, to_exclude); len_inc += 1; dec+=1;                
            end
        end
        # right
        if explore_right
            ic_col_right, to_exclude = explore_sides!(inc+1, L, type_real, positions, comp,
                                      included, num_included, 
                                      data_matrix, bg, mlen, 
                                      pseudocount; right=true);
        
            if ic_col_right ≤ ic_extension_thresh
                explore_right = false 
            else 
                to_exclude!(included, to_exclude); len_inc += 1; inc+=1;                  
            end
        end
        (!explore_left && !explore_right) && break;
        (num_included ≤ usage_lowerbdd) && break;
    end    
    
    (dec > 0 || inc > 0) && (new_positions = Dict{type_int,type_int}(
        k=>(comp[k] ? positions[k]-inc : positions[k]-dec)
        for (ind,k) in enumerate(keys(positions)) if included[ind]);)
            
    return new_positions, mlen+len_inc
end

function extend!(g::good_stuff{T,S,Q}; re_evaluate_pfm=true,
                    re_evaluate_pwm=false,
                    re_evaluate_thresh=false,
                    re_evaluate_scores=false, 
                    smoothing=true) where {T,S,Q}
    @inbounds for l = 1:g.ms.num_motifs
        ps = g.search.pseudocounts;
        if !isnothing(g.ms.use_comp)
            g.ms.positions[l], g.ms.lens[l] = extension!(g.data.L, 
                                g.ms.positions[l], g.ms.lens[l], g.ms.use_comp[l], 
                                g.data.data_matrix, g.search.bg, T, S, 
                                ps;
                                perc_no_less_than=g.search.perc_no_less_than, 
                                ic_extension_thresh=g.search.ic_extension_thresh);
        else
            g.ms.positions[l], g.ms.lens[l] = extension!(g.data.L, 
                            g.ms.positions[l], g.ms.lens[l], g.data.data_matrix,
                            g.search.bg, T, S, ps;
                            perc_no_less_than=g.search.perc_no_less_than, 
                            ic_extension_thresh=g.search.ic_extension_thresh);
        end
    end

    indicator = (length.(g.ms.positions) .> 0) .& (g.ms.lens .> g.cdl.filter_size+2);
    g.ms.positions=g.ms.positions[indicator];
    g.ms.pfms=g.ms.pfms[indicator];
    g.ms.lens=g.ms.lens[indicator];
    g.ms.num_motifs=length(g.ms.pfms);
    !isnothing(g.ms.use_comp) && (g.ms.use_comp = g.ms.use_comp[indicator])

    re_evaluations!(g, 
                    re_evaluate_pfm, 
                    re_evaluate_pwm, 
                    re_evaluate_thresh, 
                    re_evaluate_scores;
                    smoothing=smoothing);                   
end

