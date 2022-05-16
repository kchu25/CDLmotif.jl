function cdl_optimize(g::good_stuff{T,S,Q}) where {T,S,Q}
    cover_starts = ones(S, g.data.N); # TODO: maybe eliminate this line later
    ar = activation_record{T, S}(g.data.N, g.data.L, g.cdl.M);
    for train_num = 1:g.cdl.num_train
        _cdl_ = init_cdl(g.cdl, g.data);        
        for iter = 1:g.cdl.iter_num 
            code_update!(_cdl_, g.cdl.csc_iter_num, Q);
            update_best!(_cdl_, g.data.data_matrix_gpu);
            dict_update!(_cdl_, g.cdl.cdl_iter_num, Q);
            iter % 5 == 0 && println("obj: ", round(objective(_cdl_, g.data.data_matrix_gpu),digits=2));            
        end
        code_thresh = g.cdl.std_threshold_multiplier*std(_cdl_.X_best[_cdl_.X_best .> 0]);
        record_training_result!(_cdl_, 
                                code_thresh, 
                                [g.cdl.filter_size for _ = 1:g.cdl.M],
                                cover_starts,
                                ar, train_num
                                );
        _cdl_ = nothing; 

        avail_bytes, total_bytes = Mem.info();
        if avail_bytes/total_bytes < 0.5 
            # println("clean up gpu mem");        
            # CUDA.memory_status()
            GC.gc(); GC.gc(true); GC.gc(false); CUDA.reclaim(); 
            GC.gc(); GC.gc(true); GC.gc(false); CUDA.reclaim();
        end

    end
    @time at = get_spatial_info!(ar); 
    return ar, at
end

function get_ssc(g::good_stuff{T,S,Q}) where {T,S,Q}
    ssc = sum([sum([g.ms.scores[k][n] 
                                for n in keys(g.ms.scores[k])])
                                for k = 1:g.ms.num_motifs]) 
    # println("SSC: $(round(ssc,digits=1)), # pfms: $(length(t.ms.pfms))");
    return ssc
end

# save_found_results_sim(g)

function find_motif(g::good_stuff)

    ar, at = cdl_optimize(g);

    avail_bytes, total_bytes = Mem.info();
    if avail_bytes/total_bytes < 0.5 
        GC.gc(); GC.gc(true); GC.gc(false); CUDA.reclaim(); 
        GC.gc(); GC.gc(true); GC.gc(false); CUDA.reclaim();
    end
    
    @time get_pair_pfms!(at, ar, g);    
    
    @time for _ = 1:g.search.default_max_iter
        scan_w_gpu!(g, g.data.data_matrix_gpu; re_evaluate_pfm=false,
                                               re_evaluate_pwm=false,
                                               re_evaluate_thresh=false);
        extend!(g; re_evaluate_pfm=true,
                   re_evaluate_pwm=false,
                   re_evaluate_thresh=false,
                   smoothing=true);
        masked_trim!(g; re_evaluate_pfm=false,
                        re_evaluate_pwm=true,
                        re_evaluate_thresh=true);
    end

    @time begin 
        filter_using_evalue!(g);
        allr_merge!(g; re_evaluate_pwm=false,
                       re_evaluate_thresh=false);
        masked_trim!(g; re_evaluate_pfm=false,
                        re_evaluate_pwm=true,
                        re_evaluate_thresh=true,
                        update_positions=false);
    end

    # set up a higher ic threshold as most extensions is done beforehand
    g.search.ic_extension_thresh=0.6;

    @time begin 
        last_ssc = Inf; cur_ssc = nothing;
        for j = 1:g.search.max_score_iter
            # println("$j-th iter")
            overlapping_scan!(g)
            non_overlap_scan!(g; re_evaluate_pfm=false, 
                                re_evaluate_pwm=false, 
                                re_evaluate_thresh=false)
            extend!(g; re_evaluate_pfm=true,
                    re_evaluate_pwm=false,
                    re_evaluate_thresh=false,
                    smoothing=false)
            masked_trim!(g; re_evaluate_pwm=true,
                            re_evaluate_thresh=true)
            cur_ssc = get_ssc(g);
            abs(last_ssc - cur_ssc) < g.search.entropy_score_tol && break;
            last_ssc = cur_ssc;        
        end
    end

    @time begin
        last_ssc = Inf; cur_ssc = nothing;
        filter_using_evalue!(g; cpu=true, non_overlap=true);
        allr_merge!(g)
        rid_of_len_less_than_six!(g)
        for _ = 1:g.search.max_score_iter
            overlapping_scan!(g)
            non_overlap_scan!(g;re_evaluate_pfm=true, 
                                re_evaluate_pwm=true, 
                                re_evaluate_thresh=true,
                                smoothing=false,
                                less_pseudocount=true)
            cur_ssc = get_ssc(g);
            abs(last_ssc - cur_ssc) < g.search.entropy_score_tol && break;
            last_ssc = cur_ssc;
        end
    end
end

function try_to_find_motif(data; 
                           number_trials=10,
                           pval_thresh_inc=5e-4,
                           eval_thresh_inc_1=2.5e-1,
                           eval_thresh_inc_2=1e-2)                              
    g=nothing; motif_found = false; 
    pval_thresh = dat_t(2.7e-4);     
    eval_thresh1 = dat_t(5e-1);                                       
    eval_thresh2 = dat_t(1e-3);
    num_trial_left = number_trials;
    while(num_trial_left > 0)
        try
            g = good_stuff{int_t,dat_t,cdl_c}(data, 
                                              CDL_setup{int_t, dat_t, cdl_c}(), 
                                              SEARCH_setup{int_t, dat_t}());
            g.search.pval_thresh = pval_thresh;    
            g.search.e_value_thresh_1 = eval_thresh1;
            g.search.e_value_thresh_2 = eval_thresh2;                                        
            @time find_motif(g)
            motif_found = true;
            break
        catch e
            if isa(e, ArgumentError)
                pval_thresh += pval_thresh_inc;
                eval_thresh1 += eval_thresh_inc_1;
                eval_thresh2 += eval_thresh_inc_2;
                println("Did not find motifs")
                println("Relax search p-value to $(round(pval_thresh,digits=4)) for the score threshold");
                println("Relax search e-value cutoff-1 to $(round(eval_thresh1,digits=4)) ");
                println("Relax search e-value cutoff-2 to $(round(eval_thresh2,digits=4)) ");
                num_trial_left -= 1;
            else
                break
            end        
        end            
    end
    if motif_found
        return g
    else
        return nothing
    end
end


# j=5;
# while(j>0)
#     try
#         g = good_stuff{int_t,dat_t,cdl_c}(data, cdl_setup, search_setup);
#         @time find_motif(g)
#         break
#     catch e
#         if isa(e, ArgumentError)
#             println("haha")
#             j-=1
#         else
#             break
#         end        
#     end
# end