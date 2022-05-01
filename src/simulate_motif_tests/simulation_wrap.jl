
######### gapped motifs tests ##########
const mg1 = ([8,8],[4]);
const mg2 = ([9,9],[4]);
const mg3 = ([9,9],[5]);
const mg4 = ([8,8],[5]);
const mg5 = ([10,10],[3]);
const mg6 = ([13,13],[3]);
const mg7 = ([15,15],[4]);
const mg8 = ([17,17],[4]);
const mg9 = ([10,10],[5]);
const mg10 = ([13,13],[5]);
const mg11 = ([15,15],[5]);
const mg12 = ([17,17],[6]);
const mg13 = ([9,9,9],[4,4]);
const mg14 = ([12,12],[3]);
const mg15 = ([12,12],[6]);
const mg16 = ([11,10,9],[5,5]);
const mg17 = ([11,10,9],[4,3]);
const mg18 = ([11,10,10],[3,6]);
const mg19 = ([10,10,10],[3,3]);
const mg20 = ([10,10,10],[5,5]);

const gap_tests0 =[mg1, mg2];
const gap_tests = [mg1, mg2, mg3, mg4, mg5, mg6, mg7, mg8, mg9, mg10, 
        mg11, mg12, mg13, mg14, mg15, mg16, mg17, mg18, mg19, mg20];
# const gap_tests2 = [mg6]; # ,mg7,mg11,mg10,mg14,mg15
const gap_tests2 = [mg6,mg7,mg11,mg10,mg14,mg15]; # ,
const gap_tests3 = [mg14]; # ,
########################################

######### mixture no gap ###############
const mk1 = [1:8, 5:12];           # 8 long, overlap 4
const mk2 = [1:9, 6:14];           # 9 long, overlap 4
const mk3 = [1:10, 6:15];          # 10 long, overlap 5
const mk4 = [1:12, 7:18];          # 12 long, overlap 6
const mk5 = [1:12, 6:17];          # 12 long, overlap 7 
const mk6 = [1:15, 8:22];          # 15 long, overlap 8
const mk7 = [1:17, 10:26];         # 17 long, overlap 8
const mk8 = [1:17, 8:24];          # 17 long, overlap 10
const mk9 = [1:10, 6:15, 11:20];   # 10 long, overlap 5
const mk10 = [1:12, 6:18, 12:24];  # 12 long, overlap 6 
const mk11 = [1:13, 8:19, 14:26];  # 13 long, overlap 6
const mk12 = [1:15, 6:20, 10:25];  # 15 long, overlap 10
const mk13 = [1:16, 9:26];         # 17 long, overlap 8
const mk14 = [1:16, 7:22];        
const mk15 = [1:13, 8:20];
const mixture_test = [mk1, mk2, mk3, mk4, mk5, mk6, mk7, mk8, mk9, 
            mk10, mk11, mk12, mk13, mk14, mk15];
const mixture_test_triples = [mk9, mk10, mk11, mk12];
const mixture_test_s = [mk14];
const mixture_test_s2 = [mk15];
########################################

######### mixture with gap #############

const kg1 = ([5,5,5,5],[0,4,0],[1:2,1:4,3:4])
const kg2 = ([6,6,6,6],[0,4,0],[1:2,1:4,3:4])
const kg3 = ([7,7,7,7],[0,4,0],[1:2,1:4,3:4])
const kg4 = ([8,8,8,8],[0,4,0],[1:2,1:4,3:4])
const kg5 = ([9,9,9,9],[0,4,0],[1:2,1:4,3:4])
const kg6 = ([5,5,5,5],[0,6,0],[1:2,1:4,3:4])
const kg7 = ([6,6,6,6],[0,6,0],[1:2,1:4,3:4])
const kg8 = ([7,7,7,7],[0,6,0],[1:2,1:4,3:4])
const kg9 = ([8,8,8,8],[0,6,0],[1:2,1:4,3:4])
const kg10 = ([9,9,9,9],[0,6,0],[1:2,1:4,3:4])
const kg11 = ([8,8,8,8],[4,0,4],[1:3,2:3,2:4]);
const kg12 = ([9,9,9,9],[4,0,4],[1:3,2:3,2:4]);
const kg13 = ([10,10,10,10],[4,0,4],[1:3,2:3,2:4]);
const kg14 = ([8,8,8,8],[6,0,6],[1:3,2:3,2:4]);
const kg15 = ([9,9,9,9],[6,0,6],[1:3,2:3,2:4]);
const kg16 = ([10,10,10,10],[6,0,6],[1:3,2:3,2:4]);
const kg17 = ([10,10], [5], [1:1,1:2,2:2])
const kg18 = ([15,15], [7], [1:1,1:2,2:2])
const kg19 = ([9,13,9],[4,4], [1:2, 2:2, 2:3, 1:3])
const kg20 = ([13,9,13],[4,4], [1:2, 2:2, 2:3, 1:3])

const kg_tests = [kg1, kg2, kg3, kg4, kg5, 
            kg6, kg7, kg8, kg9, kg10, 
            kg11, kg12, kg13, kg14, 
            kg15, kg16, kg17, kg18, 
            kg19, kg20];
 
########################################

# function use_default(data::Simulated_DNA_Dataset,target_folder_expr::String)
#     cdl_setup = CDL_setup{int_t, dat_t, cdl_c}();
#     search_setup = SEARCH_setup{int_t, dat_t}();
#     g = good_stuff{int_t,dat_t,cdl_c}(data, cdl_setup, search_setup);
#     @time find_motif(g)
#     save_found_results_sim(target_folder_expr, g)
# end

function use_default(data::Simulated_DNA_Dataset,target_folder_expr::String)
    g = try_to_find_motif(data);
    save_found_results_sim(target_folder_expr, g, data)
end

function test_mixture_triple(output_folder::String, 
                             num_seq=num_seq, 
                             num_expr_each=num_expr_each; 
                             data_pt_len::Integer=100, 
                             bern_prob::Real=1.0)
    for _range_ in mixture_test_triples, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        m_label = join(split(split(string(_range_),"}")[2], ", "),"_");
        expr_name = "mixture_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/mixture/"*expr_name;
        _m_ = mixture_k_parts_motifs(_range_);                
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end

function test_mixture_s(output_folder::String, 
                      num_seq=num_seq, 
                      num_expr_each=num_expr_each; 
                      data_pt_len::Integer=100, 
                      bern_prob::Real=1.0)
    for _range_ in mixture_test_s, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        m_label = join(split(split(string(_range_),"}")[2], ", "),"_");
        expr_name = "mixture_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/mixture/"*expr_name;
        _m_ = mixture_k_parts_motifs(_range_);                
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end

function test_mixture_s2(output_folder::String, 
                      num_seq=num_seq, 
                      num_expr_each=num_expr_each; 
                      data_pt_len::Integer=100, 
                      bern_prob::Real=1.0)
    for _range_ in mixture_test_s2, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        m_label = join(split(split(string(_range_),"}")[2], ", "),"_");
        expr_name = "mixture_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/mixture/"*expr_name;
        _m_ = mixture_k_parts_motifs(_range_);                
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end


function test_mixture(output_folder::String, 
                      num_seq=num_seq, 
                      num_expr_each=num_expr_each; 
                      data_pt_len::Integer=100, 
                      bern_prob::Real=1.0)
    for _range_ in mixture_test, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        m_label = join(split(split(string(_range_),"}")[2], ", "),"_");
        expr_name = "mixture_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/mixture/"*expr_name;
        _m_ = mixture_k_parts_motifs(_range_);                
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end

function test_single(output_folder::String, 
                     num_seq::Integer=num_seq, 
                     num_expr_each::Integer=num_expr_each; 
                     data_pt_len::Integer=100, 
                     bern_prob::Real=1.0)
    for l in 8:22, expr_k = 1:num_expr_each     
        println("======================================"); println("working on: $l")
        # save the data first
        expr_name = "single_$(l)_$(expr_k)"
        target_folder_expr = output_folder*"/single/"*expr_name;
        # data preparation
        _m_ = single_part_motif(l);                        
        data = Sim_DNA{int_t, dat_t}(_m_,num_seq,target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end
end

function test_single_short(output_folder::String, 
                     num_seq::Integer=num_seq, 
                     num_expr_each::Integer=num_expr_each; 
                     data_pt_len::Integer=100, 
                     bern_prob::Real=1.0)
    for l in 7:9, expr_k = 1:num_expr_each     
        println("======================================"); println("working on: $l")
        # save the data first
        expr_name = "single_$(l)_$(expr_k)"
        target_folder_expr = output_folder*"/single/"*expr_name;
        # data preparation
        _m_ = single_part_motif(l);                        
        data = Sim_DNA{int_t, dat_t}(_m_,num_seq,target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end
end

function test_single_long(output_folder::String, 
                     num_seq::Integer=num_seq, 
                     num_expr_each::Integer=num_expr_each; 
                     data_pt_len::Integer=100, 
                     bern_prob::Real=1.0)
    for l in 19:22, expr_k = 1:num_expr_each     
        println("======================================"); println("working on: $l")
        # save the data first
        expr_name = "single_$(l)_$(expr_k)"
        target_folder_expr = output_folder*"/single/"*expr_name;
        # data preparation
        _m_ = single_part_motif(l);                        
        data = Sim_DNA{int_t, dat_t}(_m_,num_seq,target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end
end

function test_gapped(output_folder::String, 
                     num_seq=num_seq, 
                     num_expr_each=num_expr_each; 
                     data_pt_len::Integer=100, 
                     bern_prob::Real=1.0)
    for _range_ in gap_tests, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        strs = split(string(_range_),")");
        m_label = join(split(strs[1][2:end], ", "),"_");
        expr_name = "gaps_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/gap/"*expr_name;
        # data preparation
        _m_ = gapped_k_parts_motif(_range_...);                
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end

function test_gapped0(output_folder::String, 
                     num_seq=num_seq, 
                     num_expr_each=num_expr_each; 
                     data_pt_len::Integer=100, 
                     bern_prob::Real=1.0)
    for _range_ in gap_tests0, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        strs = split(string(_range_),")");
        m_label = join(split(strs[1][2:end], ", "),"_");
        expr_name = "gaps_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/gap/"*expr_name;
        # data preparation
        _m_ = gapped_k_parts_motif(_range_...);                
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end


function test_gapped2(output_folder::String, 
                     num_seq=num_seq, 
                     num_expr_each=num_expr_each; 
                     data_pt_len::Integer=100, 
                     bern_prob::Real=1.0)
    for _range_ in gap_tests2, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        strs = split(string(_range_),")");
        m_label = join(split(strs[1][2:end], ", "),"_");
        expr_name = "gaps_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/gap/"*expr_name;
        # data preparation
        _m_ = gapped_k_parts_motif(_range_...);                
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end

function test_gapped3(output_folder::String, 
                     num_seq=num_seq, 
                     num_expr_each=num_expr_each; 
                     data_pt_len::Integer=100, 
                     bern_prob::Real=1.0)
    for _range_ in gap_tests3, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        strs = split(string(_range_),")");
        m_label = join(split(strs[1][2:end], ", "),"_");
        expr_name = "gaps_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/gap/"*expr_name;
        # data preparation
        _m_ = gapped_k_parts_motif(_range_...);                
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end

function test_kg(output_folder::String, 
                 num_seq=num_seq, 
                 num_expr_each=num_expr_each; 
                 data_pt_len::Integer=100, 
                 bern_prob::Real=1.0)
    for _range_ in kg_tests, expr_k = 1:num_expr_each
        println("======================================"); println("working on: $_range_")
        # save the data first
        strs = split(string(_range_),", UnitRange{Int64}");
        str1 = join(split(strs[1][2:end], ", "),"_");
        str2 = join(split(strs[2][1:end-1],", "), "_");
        m_label = string(str1, "_", str2);

        expr_name = "kparts_gaps_$(m_label)_$(expr_k)"
        target_folder_expr = output_folder*"/kg/"*expr_name;
        # data preparation
        _m_ = mixture_gapped_k_parts_motifs(_range_...);      
        data = Sim_DNA{int_t, dat_t}(_m_, num_seq, target_folder_expr,data_pt_len,bern_prob);
        use_default(data, target_folder_expr);
    end        
end


