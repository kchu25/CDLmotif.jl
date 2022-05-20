function save_found_motifs(target_folder::String, ms)
    for (ind,pfm) in enumerate(ms.pfms) 
        CSV.write(target_folder*"/d_pfm$ind.csv",  Tables.table(pfm), writeheader=false) 
    end
end

const pkg_dir = pkgdir(CDLmotif);
const script_loc_template = pkg_dir*"/scripts_python/jinja_templates/";

function save_found_results_sim(output_folder::String, g::Union{good_stuff, Nothing}, data)        
    if !isnothing(g)
        # save the discovered motifs for PWM plots
        save_found_motifs(output_folder, g.ms);
        # make the csv on ground truth motifs to show how they bind
        make_motif_type_file(g.data.motif, output_folder);
        # get binding info, score contribution, and perf-coeff-d; return the max lr score
        max_lr_score = gt_cover_by(g, output_folder);
        # get the performance coefficient
        perf_coeff = get_perform_coeff(g.ms, g.data.motif, g.data);
        println("perf coeff: $perf_coeff")      
        ###### save the pickle ###########################
        # data_matrix = g.data.data_matrix; 
        gms = g.data.motif; dms = g.ms;
        @save output_folder*"/gt_motif.jld2" gms
        @save output_folder*"/d_motif.jld2" dms
        # @save output_folder*"/gt_data_matrix.jld2" data_matrix             
        script_loc_logo = pkg_dir*"/scripts_python/get_logos.py";
        script_loc_render = pkg_dir*"/scripts_python/render.py";  
        run(`python3 $script_loc_logo $output_folder`);
        run(`python3 $script_loc_render $output_folder $perf_coeff $max_lr_score $(g.data.prob_per_seq) $script_loc_template`);
    else
        # script_loc_render = pkg_dir*"/scripts_python/render_not_found.py";  
        # run(`python3 $script_loc_render $output_folder $(data.prob_per_seq) $script_loc_template`);
    end
end
################ visualization for co-occuring parts ######################
function distance_plot(target_folder::String,
                       g::good_stuff
    )
    num_discover = length(g.ms.pfms);
    for i = 1:num_discover
        d_i_info = [];
        for j = 1:num_discover
            if i != j
                vs = Int[];
                for k in keys(g.ms.positions[i])
                    if haskey(g.ms.positions[j],k)
                        diff = g.ms.positions[j][k] - g.ms.positions[i][k];
                        dist = nothing;
                        if diff > 0
                            dist = g.ms.positions[j][k]-1-(g.ms.positions[i][k]+g.ms.lens[i]-1);
                        elseif diff < 0
                            dist = -(g.ms.positions[i][k]-1-(g.ms.positions[j][k]+g.ms.lens[j]-1));
                        end
                        push!(d_i_info, (discovered="D$j", distance=dist));
                    end    
                end
            end        
        end
        df = DataFrame(d_i_info); rename!(df, :distance => :num_nucleotides_apart);
        gadfly_plot = Gadfly.plot(df, x="num_nucleotides_apart", 
                        y="discovered", Geom.histogram2d,
                        Scale.color_continuous(colormap=Scale.lab_gradient("white", "black")));
        draw(PNG(target_folder*"/distance_D$i.png", 4inch, 3inch), gadfly_plot);
    end
end


function get_occ_pairs!(target_folder::String, 
                       nuc_apart::Integer, 
                       g::good_stuff,
                       max_occ_num::Integer,
                       max_rows::Integer;
                       omit_num::Integer=30
    )                   
    record_d = Dict{String, Integer}();
    for i = 1:g.data.N
        this_i = Int[];
        this_i_motif_ind = Integer[];
        for (ind,p) in enumerate(g.ms.positions)
            if haskey(p, i)
                push!(this_i, p[i]);
                push!(this_i_motif_ind, ind);
            end
        end
        perm_ind = sortperm(this_i);
        occured_pattern = this_i_motif_ind[perm_ind];
        occured_locations = this_i[perm_ind];
        len_occ = length(occured_pattern);
        for j = 1:len_occ
            for k = j+1:len_occ                
                occured_locations_j_end = occured_locations[j]+g.ms.lens[occured_pattern[j]]-1;
                nucleotides_apart = occured_locations[k]-occured_locations_j_end;
                @assert nucleotides_apart ≥ 0 "must be non-overlapping: $i, $j, $k"
                _key_ = "(D$(occured_pattern[j]), D$(occured_pattern[k]))";
                if nucleotides_apart < nuc_apart
                    haskey(record_d, _key_) ? record_d[_key_] += 1 : record_d[_key_] = 1;                
                end                
            end
        end
    end
    df = DataFrame([(pair=k,occur=v) for (k,v) in record_d]);
    !isempty(df) && filter!(:occur => x -> x > omit_num, df)
    if !isempty(df)         
        sort!(df,[:occur], rev=true);
        max_occ_num = df.occur[1] > max_occ_num ? df.occur[1] : max_occ_num;
        max_rows = size(df,1) > max_rows ? size(df,1) : max_rows;
        CSV.write(target_folder*"/pair_occ_within_$(nuc_apart).csv", df)
    end
    return max_occ_num, max_rows
end


function save_jaspar(target_folder::String, g::good_stuff)
    if !isnothing(g)           
        ###### html and d3 stuff #########################    
        !isdir(target_folder) && mkdir(target_folder);
        save_pfm_msa_fasta(target_folder, g);
        # save_found_motifs(target_folder, g.ms);
        # max_lr_score = get_cover_by_fasta(g, target_folder);
        _ = get_cover_by_fasta(g, target_folder);
        N = g.data.N;
        ###### save    ###########################
        script_loc_render = pkg_dir*"/scripts_python/render_jaspar.py";
        script_loc_logo = pkg_dir*"/scripts_python/get_logos.py";
        ###### plot the distance information and save ####
        length(g.ms.pfms) > 1 && distance_plot(target_folder, g);
        max_occ_num = 0; max_rows = 0;
        for nuc_apart = 5:30
            max_occ_num, max_rows = get_occ_pairs!(target_folder, nuc_apart, g, 
                                        max_occ_num, max_rows);
        end        
        u = [length(i) for i in g.ms.positions];

        run(`python3 $script_loc_logo $target_folder`)
        # run(`python3 $script_loc_render $target_folder $max_lr_score $script_loc_template $num_seq`)
        run(`python3 $script_loc_render $target_folder $script_loc_template $N $max_occ_num $max_rows $u`)
    else
        script_loc_render = pkg_dir*"/scripts_python/render_fasta_not_found.py";
        run(`python3 $script_loc_render $target_folder $script_loc_template`)
    end
end

###########################################################################


function rev_comp(str::String)
    join([DNA_complement[a] for a in reverse(str)])
end

function save_pfm_msa_fasta(target_folder::String, g::good_stuff)
    if !is_simulated_data(g.data)
        for i = 1:g.ms.num_motifs
            strings = Vector{String}();
            for k in keys(g.ms.positions[i])
                start_ = g.ms.positions[i][k];
                end_ = start_+g.ms.lens[i]-1;
                str = g.data.raw_data[k][start_:end_];
                if g.ms.use_comp[i][k]
                    push!(strings, rev_comp(str))
                else
                    push!(strings, str)
                end            
            end
            rev_comp_strings = rev_comp.(strings);
            open(target_folder*"/d$i.fa","w") do file
                for (ind,s) in enumerate(strings)
                    write(file, string(">sequence_", string(ind),"\n"));
                    write(file, string(s,"\n"))
                end
            end
            open(target_folder*"/d$(i)_c.fa","w") do file
                for (ind,s) in enumerate(rev_comp_strings)
                    write(file, string(">sequence_", string(ind),"\n"));
                    write(file, string(s,"\n"))
                end
            end
            run(`weblogo -f $(target_folder)/d$i.fa --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 17 --color-scheme classic -Y NO -o $(target_folder)/d$i.png`);
            run(`weblogo -f $(target_folder)/d$(i)_c.fa --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 17 --color-scheme classic -Y NO -o $(target_folder)/d$(i)_c.png`);
        end
    end
end

function save_pfms_as_transfac(logo_folder::String, g::good_stuff, sort_perm::Vector{Int})
    for (i,ind) in enumerate(sort_perm)
        io = open(logo_folder*"/d$i.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = g.ms.pfms[ind] .* 1000; # make it a count matrix
        for j = 1:size(g.ms.pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        run(`weblogo -D transfac -f $(logo_folder)/d$(i).transfac --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 17 --color-scheme classic -o $(logo_folder)/d$(i).png`);

        # do it for the reverse complement as well
        io = open(logo_folder*"/d$(i)_c.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = pfm_complement(g.ms.pfms[ind]) .* 1000;
        for j = 1:size(g.ms.pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        run(`weblogo -D transfac -f $(logo_folder)/d$(i)_c.transfac --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 17 --color-scheme classic -o $(logo_folder)/d$(i)_c.png`);
    end
end

function get_configurations(g::Union{good_stuff, Nothing})
    configurations = Dict{Vector{Tuple{Int, Bool}}, Int}();
    gap_configurations = Dict{Vector{Tuple{Int, Bool}}, Vector{Vector{Int}}}();
    for n = 1:g.data.N
        occured_n = Tuple{Int, Bool}[];
        position_n = Int[];
        for k = 1:g.ms.num_motifs
            if haskey(g.ms.positions[k],n)
                push!(occured_n, (k, g.ms.use_comp[k][n]));
                # push!(occured_n, (sort_perm_indices[k], g.ms.use_comp[k][n])); # since we've permuted by using e-values
                push!(position_n, g.ms.positions[k][n])            
            end
        end

        if !isempty(occured_n)
            sort_perm = sortperm(position_n);
            occured_n = occured_n[sort_perm];
            if haskey(configurations, occured_n) 
                configurations[occured_n] += 1;
            else configurations[occured_n] = 1;
            end
            if length(occured_n) > 1
                position_n = position_n[sort_perm];
                gap_vec = Int[];
                num_gaps = length(occured_n)-1;
                for i = 1:num_gaps
                    # which_one = inv_sort_perm_map[occured_n[i][1]];
                    gap_len = position_n[i+1]-(position_n[i]+g.ms.lens[occured_n[i][1]]-1);
                    push!(gap_vec, gap_len)
                end
                @assert all(gap_vec .≥ 0) "gap between PWM positions should not be negative"
                if haskey(gap_configurations, occured_n)
                    push!(gap_configurations[occured_n], gap_vec);
                else gap_configurations[occured_n] = [gap_vec];
                end
            end
        end
    end
    return configurations, gap_configurations
end

tufte_bar = Theme(
    default_color=colorant"gray",
    background_color=colorant"white", 
    highlight_width=0.3mm,
    bar_spacing=15pt,
    major_label_font_size=24pt,
    point_label_font_size=15pt,
    minor_label_font_size=15pt,
    # key_label_color="gray",
    # minor_label_color="gray",
    # major_label_color="gray",
    # point_label_color="gray",
    grid_line_width=1pt, 
    minor_label_font="Times",
    major_label_font="Times", 
)

kde_theme = Theme(
    line_width=1.25mm,
    key_label_color="black",
    key_label_font_size=25pt,
    minor_label_font_size=35pt,
    major_label_font_size=55pt,
    minor_label_font="Times",
    major_label_font="Times", 
)

function map_key_update(selected_keys, sort_perm_map)
    selected_keys_update = typeof(selected_keys)();
    for i = 1:length(selected_keys)
        push!(selected_keys_update, [(sort_perm_map[j[1]], j[2]) for j in selected_keys[i]]);
    end
    selected_keys_update
end

function make_configurations_plots_2(configurations, 
                                   gap_configurations,
                                   pics_folder,
                                   sort_perm_map; 
                                   threshold_weight=0.001,
                                   kde_bandwidth=1.25)
    config_counts = [configurations[k] for k in keys(configurations)];
    config_weights = config_counts ./ sum(config_counts);
    # take only configurations that have "significant" weights
    # might do a more rigorious way to evaluate what is significant later
    # for now, simply discard all the configuration that has weight less than threshold_weight
    selected_keys = [k for k in keys(configurations)][config_weights .> threshold_weight];
    selected_configurations = Dict(k=>configurations[k] for k in selected_keys);
    selected_gap_configurations = Dict(k=>gap_configurations[k] 
                                for k in selected_keys if length(k) > 1); 
    # re-estimate the weights                                                                      
    selected_config_counts = [selected_configurations[k] for k in selected_keys];
    selected_config_weights = selected_config_counts ./ sum(selected_config_counts);
    # sort the patterns according to their weights
    sp_ind = sortperm(selected_config_weights, rev=true);
    selected_config_weights_sorted = selected_config_weights[sp_ind];
    selected_config_counts_sorted = selected_config_counts[sp_ind];
    selected_keys_sorted = selected_keys[sp_ind];
    selected_keys_sorted_update = map_key_update(selected_keys_sorted, sort_perm_map)
    
    # plot the gap kde plots
    # from the largest weight to the smallest weight
    for ind = 1:length(sp_ind)
        key = selected_keys_sorted[ind];
        for gap_bt_ind = 1:(length(key)-1)     
            # println("ind:", ind)
            # println("length key: $(length(key))")
            # remember gap_config is a vector of vector of integers...
            # [[gap_len_1, gap_len_2]₁, [gap_len_1, gap_len_2]₂, ...]
            q=[i[gap_bt_ind] for i in selected_gap_configurations[key]];
            p=Gadfly.plot(DataFrame(q=q), x=:q, 
                          Geom.density(bandwidth=kde_bandwidth),
                          Guide.xlabel("Number of nucleotides in between"), 
                          Guide.xticks(ticks=0:10:maximum(q)),
                          kde_theme);
            # save plot
            draw(PNG(pics_folder*"/gap_$(gap_bt_ind)_mode_$ind.png", 12inch, 7inch), p)
        end        
    end
    
    # plot the weights
    # patterns = [join([j[2] ? "$(j[1])" : "rc($(j[1]))" for j in i],", ") for i in keys(selected_configurations)]

    patterns = [join([!j[2] ? "D$(j[1])" : "rc(D$(j[1]))" for j in i],", ") 
                    for i in selected_keys_sorted_update]
    df=DataFrame(pattern=patterns, weights=selected_config_weights_sorted); df = sort(df, :weights);
    p = Gadfly.plot(df, y=:pattern, x=:weights,
        Geom.bar(orientation=:horizontal),
        Guide.xlabel("Weights", orientation=:horizontal),
        Guide.ylabel("    "),
        # Guide.ylabel("Binding Patterns"),
        Guide.title("Enriched Patterns and their weights"),
        tufte_bar);
    draw(PNG(pics_folder*"/mb.png", 12inch, 12inch), p)
    return selected_keys_sorted_update, 
           selected_config_weights_sorted, 
           selected_config_counts_sorted,
           maximum(length.(selected_keys_sorted_update))-1
end    


function save_result_fasta(g::Union{good_stuff, Nothing}, target_folder::String)
    logo_folder_name = "logos"
    pics_folder_name = "other_pics"
    logo_folder = target_folder*"/"*logo_folder_name
    pics_folder = target_folder*"/"*pics_folder_name
    if !isnothing(g)
        # make folders
        !isdir(target_folder) && (mkdir(target_folder);)
        !isdir(logo_folder) && (mkdir(logo_folder);)
        !isdir(pics_folder) && (mkdir(pics_folder);)
        ###### save msa file for each pwm ################    
        evalues = filter_using_evalue!(g; cpu=true, non_overlap=true, get_evalues_only=true);            
        sort_perm = sortperm(evalues); # sort it according to e-values (small to big)
        sort_perm_map = Dict(item=>index for (index,item) in enumerate(sort_perm));
        evalues = round.(evalues[sort_perm], sigdigits=3);
        scores = [round(sum(values(g.ms.scores[i])),digits=1)
                    for i = 1:g.ms.num_motifs][sort_perm]; # scores of each discovered motifs
        msa_counts = [length(g.ms.positions[j]) for j = 1:g.ms.num_motifs][sort_perm];

        configurations, gap_configurations = get_configurations(g);
        selected_keys_sorted_update, weights_sorted, count_sorted, max_gap_num = make_configurations_plots_2(configurations, gap_configurations, pics_folder, sort_perm_map);
        weights_sorted = round.(weights_sorted, digits=3);

        labels = ["D$j" for j = 1:g.ms.num_motifs];
        logos = ["d$(j)" for j = 1:g.ms.num_motifs];
            
        save_pfms_as_transfac(logo_folder, g, sort_perm);

        # render the main result page
        df = DataFrame(label=labels, count=msa_counts, slrs=scores, eval=evalues, logo=logos);
        out = Mustache.render(html_template, target_folder=target_folder, logo_folder=logo_folder_name, num_seq=g.data.N, DF=df);
        ###### html stuff #########################            
        io = open(target_folder*"/summary.html", "w")
        print(io, out);
        close(io)

        # render the binding pattern page
        tpl=create_template_bindind_patterns(max_gap_num)
        out = Mustache.render(tpl, DF=get_df_for_binding_pattern(selected_keys_sorted_update, 
                                                                 logo_folder_name, 
                                                                 pics_folder_name,
                                                                 weights_sorted,
                                                                 count_sorted))
        io = open(target_folder*"/bp.html", "w")
        print(io, out);
        close(io)
    end
end

