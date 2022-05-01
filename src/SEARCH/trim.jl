function trim_this_pfm_by_masking!(pfm, bg, mask_len, mask_thresh)
    ic_col = get_ic_col_multiple(pfm, bg); len_ = size(ic_col,2);                               
    trim_info = fill(false, len_); 
    @inbounds for (mlen,mthresh) in zip(mask_len,mask_thresh)
        for i =1:len_-mlen+1
            if trim_info[i] 
                continue
            elseif median(ic_col[i:i+mlen-1]) < mthresh
                trim_info[i] = true;
            else
                break
            end
        end

        for j = len_:-1:mlen
            if trim_info[j] 
                continue
            elseif median(ic_col[j:-1:j-mlen+1]) < mthresh
                trim_info[j] = true;
            else
                break
            end
        end
    end
    front_trim = 0; back_trim = 0;
    for j = 1:length(trim_info) 
        if trim_info[j] == 1 front_trim += 1; else break; end
    end
    for j = length(trim_info):-1:1
        if trim_info[j] == 1 back_trim += 1; else break; end
    end
    if front_trim > 0 || back_trim > 0
        pfm = pfm[:, (front_trim+1):(end-back_trim)];
    end
    return pfm, front_trim, back_trim
end

function masked_trim!(g; re_evaluate_pfm=false,
                         re_evaluate_pwm=true,
                         re_evaluate_thresh=true,
                         re_evaluate_scores=false,
                         update_positions=true,
                         after_merged=false
                         )

    @inbounds for i = 1:length(g.ms.pfms)
        g.ms.pfms[i], front_trim, back_trim = trim_this_pfm_by_masking!(g.ms.pfms[i], g.search.bg, g.search.mask_len, g.search.mask_thresh);
        if !isnothing(g.ms.positions)
            for k in keys(g.ms.positions[i]) 
                if after_merged
                    if g.ms.use_comp[i][k]
                        g.ms.positions[i][k] += back_trim
                    else
                        g.ms.positions[i][k] += front_trim; 
                    end
                else
                    g.ms.positions[i][k] += front_trim; 
                end
            end
        end
    end
    
    # rid of motifs with len less than 7
    indicator = size.(g.ms.pfms,2) .> g.cdl.filter_size+1; 
    # rid of motifs with no activations
    if !isnothing(g.ms.positions)
        indicator_pos = length.(g.ms.positions) .> 0;
        indicator = indicator .& indicator_pos;
    end

    g.ms.pfms = g.ms.pfms[indicator];
    g.ms.lens = size.(g.ms.pfms,2);
    g.ms.num_motifs = length(g.ms.pfms);
    update_positions ? g.ms.positions=g.ms.positions[indicator] : g.ms.positions=nothing;    
    
    re_evaluations!(g, re_evaluate_pfm, re_evaluate_pwm, re_evaluate_thresh, re_evaluate_scores);          
end
