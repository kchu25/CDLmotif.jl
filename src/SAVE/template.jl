const html_template=mt"""<!DOCTYPE html>
<html>
<meta charset="utf-8">
<head>
	<title></title>
	<style>
	table, td {
		border-collapse: collapse;
		margin: 15px 15px;
		padding: 5px 5px;
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
			Explanation for each column:	
			<ul>
				<li>Label: An assigned label for each motif.</li><br>
				<li>count: The number of subsequences used in the MSA to estimate the PWM. This is \(N_j\) as defined by Likelihood ratio scores. Each subsequence in MSA comes from distinct input sequences in the dataset.</li>
				<br>
				<li>sum-LRS: Sum of the likelihood ratio scores. </li>
				<br>				
			</ul>
			<center><h4>Likelihood ratio scores:</h4></center>
			Let \(P_j\) be the position frequency matrix estimated from the \(j\)th learned motif and \(N_j\) be the number of sequences used in the estimation of \(P_j\).
			The <i>likelihood ratio score</i> of the \(j\)th learned motif of length \(L\) is 
			$$ \sum_{n=1}^{N_j}\sum_{\ell=1}^L \sum_{\alpha} \unicode{x1D7D9}\left[s_n[\ell]=\alpha\right] \, P_j[\alpha,\ell]\, \ln \frac{P_j[\alpha,\ell]}{B[\alpha]}$$
			where \(B[\alpha]\) is the background frequency of nucleotide \(\alpha\), \(s_n\) the \(n\)th substring used in estimating \(P_j\),
			and \(\unicode{x1D7D9}[\cdot]\) is the indicator function. In this experiment, \(B[\alpha]=1/4,\,\forall \alpha\).		
		</div>
		<div style="float:left; margin:10px"> 
			
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
						<td>Label</td><td>count</td><td>sum-LRS</td><td>Logo</td>
					</tr>		

                    {{#:DF}}
                    <tr>
                        <td>{{:label}}</td>
                        <td>{{:count}}</td>
                        <td>{{:slrs}}</td>
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
</body>
</html>
"""