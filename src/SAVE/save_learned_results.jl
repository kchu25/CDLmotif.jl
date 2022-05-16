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

# distance_plot(target_folder, g);


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
    for i in sort_perm   
        io = open(logo_folder*"/d$i.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = g.ms.pfms[i] .* 1000; # make it a count matrix
        for j = 1:size(g.ms.pfms[i],2)
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
        q = pfm_complement(g.ms.pfms[i]) .* 1000;
        for j = 1:size(g.ms.pfms[i],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        run(`weblogo -D transfac -f $(logo_folder)/d$(i)_c.transfac --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 17 --color-scheme classic -o $(logo_folder)/d$(i)_c.png`);
    end
end

function save_result_fasta(g::Union{good_stuff, Nothing}, target_folder::String)
    logo_folder_name = "logos"
    logo_folder = target_folder*"/"*logo_folder_name
    if !isnothing(g)
        !isdir(target_folder) && (mkdir(target_folder);)
        !isdir(logo_folder) && (mkdir(logo_folder);)
        ###### save msa file for each pwm ################    
        evalues = filter_using_evalue!(g; cpu=true, non_overlap=true, get_evalues_only=true);    

        sort_perm = sortperm(evalues); # sort it according to e-values (small to big)        
        evalues = round.(evalues[sort_perm], sigdigits=3);
        scores = [round(sum(values(g.ms.scores[i])),digits=1)
                    for i = 1:g.ms.num_motifs][sort_perm]; # scores of each discovered motifs
        msa_counts = [length(g.ms.positions[j]) for j = 1:g.ms.num_motifs][sort_perm];
        labels = ["D$j" for j = 1:g.ms.num_motifs];
        logos = ["d$(j)" for j = 1:g.ms.num_motifs];
            
        save_pfms_as_transfac(logo_folder, g, sort_perm);

        df = DataFrame(label=labels, count=msa_counts, slrs=scores, eval=evalues, logo=logos)
        out = render(html_template, target_folder=target_folder, logo_folder=logo_folder_name, num_seq=g.data.N, DF=df);
        ###### html stuff #########################            
        io = open(target_folder*"/summary.html", "w")
        print(io, out);
        close(io)
    end
end