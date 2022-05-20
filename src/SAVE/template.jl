html_template=mt"""<!DOCTYPE html>
<html>
<meta charset="utf-8">
<head>
	<title></title>
	<style>
	table, td {
		border-collapse: collapse;
		margin: 15px 15px;
		padding: 5px 5px;
        table-layout: fixed;
        min-width: 85px;
	}
	.top_row {
	  font-weight: bold;
	  color: #808080;
	}

	thead,tfoot {
		font-weight: bold;
	    background-color: #333;
	    color:white;
	}

	.info {
    	background-color: #E2E2E2;
    	margin:5px;
    	padding:5px;
    }
	</style>
	 <script id="MathJax-script" async
          src="https://cdn.jsdelivr.net/npm/mathjax@3.0.0/es5/tex-mml-chtml.js">
  </script>
</head>


<body>
	<div style="display:flex;">
		<div class="info" style="border:1px solid black; max-width: 500px">
			Explanation for the PWM table:	
			<ul>
				<li><b>Label:</b> An assigned label for each motif.</li><br>
				<li><b>Count:</b> The number of subsequences (substrings) used in the multiple sequence alignment (MSA) to estimate the PWM. Each subsequence in MSA comes from distinct input sequences in the dataset.</li><br>
				<li><b>E-value</b>: The p-value from the Fisher exact test multiplied by the number of PWMs found. The Fisher exact test test on how much a PWM is <i>activated</i> in the training dataset vs. how much a PWM is activated the control dataset. The "<i>PWM activation</i>" means that the PWM can score above a fixed score-threshold, which we calculate using methods from <a href="https://almob.biomedcentral.com/articles/10.1186/1748-7188-2-15">Touzet et al.</a>. The control dataset is constructed by <a href="https://link.springer.com/article/10.1186/1471-2105-9-192">shuffling the nucleotides in the training dataset</a>; this shuffling preserves the frequency of the 3-mers in the training dataset. </li><br>
				<ul>
					<li>The less the e-value, the more statistically significant the discovered PWM: a statistically significant PWM is much more enriched in the training dataset in comparison to the control dataset. </li>
				</ul>
			</ul>
			Explanation for the bar plot <i>Enriched Patterns and their weights</i>:
			<ul>
				<li>The labels on the y-axis are the <i><b>patterns</b></i>. The patterns show how combinations of PWMs are activated in the training dataset.</li><br>
				<ul>
					<li>Each pattern shows the order of activations of each PWM. They are expressed in the order of scanning a sequence from the left to the right. For example, a label (D5, D2) means that D5 is activated on the left of the activation of D2.</li><br>
					<li>rc(\( \cdot \)) means that the PWM was activated in its reverse complement orientation.</li><br>					
				</ul>
				<li>The weights shows how each patterns are enriched relative to each other in the training dataset. The weights sum to one.</li><br>
				<li>We discard all the patterns that have negligible weights (weight less than 0.01). The rest of the weights are renormalized and presented in the plot. </li>
				<br>
			</ul>
            Other Notes:
                <ul>
                    <li>All the positions in the training dataset that used to estimate the PWMs are mutually disjoint. 
                    </li><br>
                </ul>
				<br>
				<br>	
<!-- 			<center><h4>Likelihood ratio scores:</h4></center>
			Let \(P_j\) be the position frequency matrix estimated from the \(j\)th learned motif and \(N_j\) be the number of sequences used in the estimation of \(P_j\).
			The <i>likelihood ratio score</i> of the \(j\)th learned motif of length \(L\) is 
			$$ \sum_{n=1}^{N_j}\sum_{\ell=1}^L \sum_{\alpha} \unicode{x1D7D9}\left[s_n[\ell]=\alpha\right] \, P_j[\alpha,\ell]\, \ln \frac{P_j[\alpha,\ell]}{B[\alpha]}$$
			where \(B[\alpha]\) is the background frequency of nucleotide \(\alpha\), \(s_n\) the \(n\)th substring used in estimating \(P_j\),
			and \(\unicode{x1D7D9}[\cdot]\) is the indicator function. In this experiment, \(B[\alpha]=1/4,\,\forall \alpha\).	 -->	
		</div>
		<div style="float:left; margin:25px; border:1px solid black; max-width:500px; padding:10px;" > 			
			<table>
				<thead>
					<tr>
						<th colspan="100%">
							Discovered motifs
						</th>
					</tr>
				</thead>
			<tbody>	
					Number of input sequences: {{:num_seq}}
					<tr class="top_row">
						<td>Label</td><td>Count</td><td>E-value</td><td>Logo</td>
					</tr>		

                    {{#:DF}}
                    <tr>
                        <td>{{:label}}</td>
                        <td>{{:count}}</td>
                        <td>{{:eval}}</td>
                        <td><img id="d_logo_{{:label}}" width="165" src="{{:logo_folder}}/{{:logo}}.png"><br>
                            <div id="d_orientation_{{:label}}">Learned PWM</div><br>
							<button type="button" onclick="discovered_{{:label}}_changeToRC()">Reverse complement</button>
                        </td>
                        <script type="text/javascript">					
									function discovered_{{:label}}_changeToRC() {
										var image = document.getElementById("d_logo_{{:label}}");
										if (image.src.match("_c")) {
							                image.src = "{{:logo_folder}}/{{:logo}}.png";
							            } else {
							            	image.src = "{{:logo_folder}}/{{:logo}}_c.png";
							            }
							            var orientation = document.getElementById("d_orientation_{{:label}}");
							            if (orientation.innerHTML === "Learned PWM"){
							            	orientation.innerHTML = "Learned PWM's reverse-complement";
							            } else {
							            	orientation.innerHTML = "Learned PWM";
							            }
									}						 
                        </script>	
                    </tr>
                    {{/:DF}}		
				</tbody>
			</table>
			<br><br>	
		</div>	
        <div style="margin=25px; border:1px solid black;"> 		
            <img src="other_pics/mb.png" style="width:90%; height:auto"><br><br><br>
            &nbsp&nbsp Check <a href="bp.html">here</a> to see more details on the enriched patterns.
        </div>
    </div>
</body>
</html>
"""


html_head2="<!DOCTYPE html>
        <html>
        <meta charset=\"utf-8\">
        <head>
            <style>
                th, td {
                 padding: 15px;
                 text-align: center;
                }
                .top_row {
                font-weight: bold;
                color: #808080;
                }

                table {
                border-collapse: collapse;
                }
                tr { 
                border: solid;
                border-width: 1px 0;
                }
                tr:first-child {
                  border-top: none;
                }
                tr:last-child {
                  border-bottom: none;
                }
                .info {
                    background-color: #E2E2E2;
                    margin:5px;
                    padding:5px;
                }
               </style>
        </head>
        <div class=\"info\" style=\"border:1px solid black; max-width: 500px\">
                            test test test
        </div>
        <body>
            <table>
                    <thead></thead>
            ";

html_end="</tbody>
        </table>
        </body>"

function print_table_body(num_gap::Int)
    top_row = "<td>Weights</td><td># Occurence</td><td>Part 1</td>";
    for i = 1:num_gap
        top_row = top_row*"<td>Gap $i</td><td>Part $(i+1)</td>"
    end
    top_row = "<tr class=\"top_row\">"*top_row*"</tr>"

    cells = "<td>{{{:weights}}}</td><td>{{{:counts}}}</td><td>{{{:part_1}}}</td>";
    for i = 1:num_gap
        cells = cells*"<td>{{{:gap_$i}}}</td><td>{{{:part_$(i+1)}}}</td>"
    end
    top_row*"{{#:DF}}"*"<tr>"*cells*"</tr>{{/:DF}}"
end

function create_template_bindind_patterns(num_gap::Int)
    html_head2*print_table_body(num_gap)*html_end
end

img_str_logo(x::String) = "<img width=\"195\" src=\""*x*"\">"
img_str_gap(x::String) = "<img width=\"265\" src=\""*x*"\">"

function get_logo_str(x::Int, comp::Bool, logo_folder::String) 
    if comp
        return logo_folder*"/d$(x)_c.png"
    else
        return logo_folder*"/d$x.png"
    end
end

function get_gap_str(k::Int, i::Int, pics_folder::String)
    pics_folder*"/gap_$(k)_mode_$i.png"
end

function return_binding_pattern_name(gap_num)
    names = String[]; push!(names, "weights"); push!(names, "counts")
    k = 1;
    ks = Int[];
    for i = 1:gap_num*2+1
        push!(ks, k)
        if isodd(i)             
            push!(names, "part_$k");
        else
            push!(names, "gap_$k");
            k+=1;
        end
    end
    return names, ks
end

function get_df_for_binding_pattern(selected_keys_, logo_folder_name, pics_folder_name, weights, counts)
    gap_num = maximum(length.(selected_keys_))-1;
    column_num = gap_num*2+1;
    names, ks = return_binding_pattern_name(gap_num);

    rows = Vector{NamedTuple{Tuple(Symbol.(names))}}();
    for i = 1:length(selected_keys_)
        row = Vector{String}(); push!(row, string(weights[i])); push!(row, string(counts[i]));
        part_ind = 1; 
        for ind_j = 1:length(selected_keys_[i])*2-1
            if isodd(ind_j)
                label = selected_keys_[i][part_ind][1];
                html_string_1 = img_str_logo(
                    get_logo_str(label, 
                    selected_keys_[i][part_ind][2], 
                                 logo_folder_name));
                html_string_2 = selected_keys_[i][part_ind][2] ? "rc(D$label)" : "D$label";
                push!(row, html_string_1*"<br>"*html_string_2
                        );
                part_ind+=1
            else
                push!(row, img_str_gap(get_gap_str(ks[ind_j],i,pics_folder_name)))
            end
        end
        row_len = length(row);
        if row_len-2 < column_num
            for _ = 1:column_num-(row_len-2) push!(row, "") end
        end
        row_tuple = NamedTuple{Tuple(Symbol.(names))}(r for r in row);
        push!(rows, row_tuple);
    end
    return DataFrame(rows)
end

