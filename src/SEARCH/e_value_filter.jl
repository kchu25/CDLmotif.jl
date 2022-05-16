function filter_using_evalue!(g::good_stuff; cpu=false, non_overlap=false, get_evalues_only=false)
    if cpu 
        if non_overlap
            bg_pos, bg_score = overlapping_scan_bg!(g);
            bg_pos, _ = non_overlapping_scan_bg!(g, bg_pos, bg_score);
        else
            bg_pos, _ = overlapping_scan_bg!(g);
        end
    else
        bg_dat_gpu = cu(g.data.data_matrix_bg);        
        bg_pos = scan_w_gpu!(g, bg_dat_gpu; scan_bg=true);
        bg_dat_gpu = nothing; 
    end

    len_positions = length.(g.ms.positions);
    len_positions_bg = length.(bg_pos);    
    evalues = fill(0.0, g.ms.num_motifs);    

    @inbounds for i = 1:g.ms.num_motifs
        # multiplier = g.data.L-g.ms.lens[i]+1;
        a = len_positions[i]; b = len_positions_bg[i];
        q = FisherExactTest(promote_i(a, g.data.N, b, g.data.N)...);     
        evalues[i] = HypothesisTests.pvalue(q)*g.ms.num_motifs;
    end    
    
    if get_evalues_only
        return evalues
    else
        evalue_filtering!(g, evalues, cpu);    
    end
end

function evalue_filtering!(g::good_stuff, evalues::Vector{T}, cpu::Bool) where T <: Real
    @assert g.ms.num_motifs == length(evalues)
    pass_indicator = cpu ?  evalues .< g.search.e_value_thresh_2 : evalues .< g.search.e_value_thresh_1;
    pass_indicator = pass_indicator .& (g.ms.lens .> 5); # width of PWM must be bigger than 5
    g.ms.pfms = g.ms.pfms[pass_indicator];
    g.ms.pwms = g.ms.pwms[pass_indicator];
    g.ms.thresh = g.ms.thresh[pass_indicator];
    g.ms.lens = g.ms.lens[pass_indicator];
    g.ms.num_motifs = length(g.ms.pfms);
    g.ms.positions = g.ms.positions[pass_indicator];
    g.ms.scores = nothing;
    # g.ms.orig_len = g.ms.orig_len[pass_indicator];
end

function rid_of_len_less_than_six!(g::good_stuff)
    pass_indicator = (g.ms.lens .> 5); # width of PWM must be bigger than 5
    g.ms.pfms = g.ms.pfms[pass_indicator];
    g.ms.pwms = g.ms.pwms[pass_indicator];
    g.ms.thresh = g.ms.thresh[pass_indicator];
    g.ms.lens = g.ms.lens[pass_indicator];
    g.ms.num_motifs = length(g.ms.pfms);
end

function get_evalues_each(g::good_stuff)
    bg_pos = overlapping_scan_bg!(g);
    len_positions = length.(g.ms.positions);
    len_positions_bg = length.(bg_pos);    
    evalues = fill(0.0, g.ms.num_motifs);  
    @inbounds for i = 1:g.ms.num_motifs
        a = len_positions[i]; b = len_positions_bg[i];
        q = FisherExactTest(promote_i(a, g.data.N, b, g.data.N)...);     
        evalues[i] = HypothesisTests.pvalue(q)*g.ms.num_motifs;
    end    
    return evalues
end