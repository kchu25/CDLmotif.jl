function save_found_motifs(target_folder::String, ms)
    for (ind,pfm) in enumerate(ms.pfms) 
        CSV.write(target_folder*"/d_pfm$ind.csv",  Tables.table(pfm), writeheader=false) 
    end
end

function save_found_results_sim(g::good_stuff)    
    # save the discovered motifs for PWM plots
    save_found_motifs("expr", g.ms);
    # make the csv on ground truth motifs to show how they bind
    make_motif_type_file(g.data.motif, "expr");
    # get binding info, score contribution, and perf-coeff-d; return the max lr score
    max_lr_score = gt_cover_by(g, "expr");
    # get the performance coefficient
    perf_coeff = get_perform_coeff(g.ms, g.data.motif, g.data);
    println("perf coeff: $perf_coeff")      
    ###### save the pickle ###########################
    # data_matrix = g.data.data_matrix; 
    gms = g.data.motif; dms = g.ms;
    @save "expr/gt_motif.jld2" gms
    @save "expr/d_motif.jld2" dms
    
    script_loc_logo = "scripts_python/get_logos.py";
    script_loc_render = "scripts_python/render.py";
    script_loc_template = "scripts_python/jinja_templates/";

    run(`python3 $script_loc_logo expr`);
    run(`python3 $script_loc_render expr $perf_coeff $max_lr_score $(g.data.prob_per_seq) $script_loc_template`);
end


function save_result_fasta(g::good_stuff,target_folder::String)
    
    ###### html and d3 stuff #########################    
    !isdir(target_folder) && mkdir(target_folder);
    # save_found_motifs(target_folder, g.ms);
    max_lr_score = get_cover_by_fasta(g, target_folder);
    num_seq = g.data.N;
    ###### save the pickle ###########################
    # @save target_folder*"/"*"good_stuff.jld2" g

    script_loc_logo = "scripts_python/get_logos.py";
    script_loc_render = "scripts_python/render_fasta.py";
    script_loc_template = "scripts_python/jinja_templates/";

    run(`python3 $script_loc_logo $target_folder`)
    run(`python3 $script_loc_render $target_folder $max_lr_score $script_loc_template $num_seq`)
end