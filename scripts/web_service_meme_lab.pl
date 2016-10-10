  ### MEME web service script ###

#input argument order:
# - promoter sequence file
# - promoter length
# - cluster membership file
# - number of motifs to generate per cluster
# - minimum motif width
# - maximum motif width


BEGIN {
  push(@INC, '/scripts'); #This should be the location of Meme_Lab.pm
}

#needs the following perl packages : MooseX::Declare (cpanned), Bio::Seq (apt-getted), GD (trainwreck shell install)
#plus libgd (libgd-dev) and maybe libpng (libpng-dev)

use Meme_Lab;

my $meme_binary = '/root/meme/bin/meme';


#These are the required input parameters
my $output_tmp_dir ="output/";							#results go here. Needs trailing slash
my $core_promoter_max_length = $ARGV[1];					# (50-1000 bp) Search for motifs within the region of the promoter from the transcription start site to this number of bases upstream in the 5-prime direction. Overlapping sequences are merged to avoid redundancy. As a result some sequences in the final set may be longer than the specified maximum.
my $num_motifs = $ARGV[3];									# (1-10) Num motifs to generate for each cluster
my $min_motif_width =$ARGV[4];								# (6-20 bp) 
my $max_motif_width = $ARGV[5];								# (6-20 bp) Range for sizes allowed for motifs
my $input_file = $ARGV[2];						#cluster file consisting of 2 tab separated cols, col1 = cluster num, col2 = gene ID.
my $sequence_file = $ARGV[0]; 				#fasta file with sequnces for the gene ids in the above file
my $background_file = ""; 								#Optional file. If present must match specification defined here http://meme-suite.org/doc/bfile-format.html

my $core_promoter_min_length = 50; 				#In webtool user can't change this


if ($max_motif_width < $min_motif_width) {
	print "Maximum provided motif width smaller than minimum motif width. Reversing.\n";
	$holder = $max_motif_width;
	$min_motif_width = $max_motif_width;
	$max_motif_width = $holder;
}

if ($background_file) {
	$MEME_parameters = MEME_Parameters->new(-minw => $min_motif_width, -maxw => $max_motif_width, -nmotifs =>$num_motifs, -bfile => $background_file);
}else{
	$MEME_parameters = MEME_Parameters->new(-minw => $min_motif_width, -maxw => $max_motif_width, -nmotifs =>$num_motifs);
}


if (!-e $output_tmp_dir) {
    mkdir ($output_tmp_dir) or die 'Error creating output directories!';
}

my $logo_directory = 'logos/';
if (!-e $output_tmp_dir.$logo_directory) {
    mkdir ($output_tmp_dir.$logo_directory) or die 'Error creating output directories!';
}
my $meme_directory = $output_tmp_dir.'meme/';
if (!-e $meme_directory) {
    mkdir ($meme_directory) or die 'Error creating output directories!';
}
my $gene_list_path = 'which_gene_ids/';
if(!-e $output_tmp_dir.$gene_list_path){
    mkdir($output_tmp_dir.$gene_list_path) or die("Error creating output directories");
}



#On submission the gene IDs will be checked to see that corresponding sequences can be retrieved from the fasta file
#Any that are invalid will be removed. Should more than 50% be invalid, or should any individual cluster contain no valid IDs,
#then job submission will fail.


#Returns zero on success or non-zero on failure

my $output_file = $input_file; # NB. will overwite original input file


#First validate input file
print "Validating input file...\n";

if (lc(substr($sequence_file, -6)) ne ".fasta"){
	print "Your fasta file of corresponding promoter sequences must have the .fasta extension, or it will not be read correctly. Please rename this file and try again.\n"; 
	exit 1; # return status: FAIL		
}	

my ($result, $cluster_ids_ref , $clusters_hash_ref, $sequences_from_db_ref) = valid_clusters($input_file, $sequence_file);

if ($result == 0){
	exit 1; #invalid inputs
}


my @clusters = @$cluster_ids_ref; #list of valid cluster names in input order
my %clusters_hash = %$clusters_hash_ref; #hash of gene lists with cluster names as keys
my $num_clusters = scalar (@clusters);

my %sequences_from_db = %$sequences_from_db_ref; #array of hashes with keys NAME, SEQUENCE

 my @meme_arguments = keys %{$MEME_parameters};
 my $meme_inputs = '';
      
 foreach my $arg (@meme_arguments) {
	 	 # If argument is -revcomp or -pal or -dna then because they have no required parameters do not print the value associated with that key
	  	print "ARG= $arg\n" ;
	  	if( ($arg eq '-revcomp') or  ($arg eq '-pal') or ($arg eq '-dna') ) {
	      	if ($$MEME_parameters{$arg}) {
		  		$meme_inputs .= " $arg";	
	      	}
	  	}
	  	else{  
	      	$meme_inputs .= " $arg $$MEME_parameters{$arg}";	 
	  	}			
}
print $meme_inputs;

my %all_meme_results = {};

foreach my $cluster (@clusters) {
	
        print "\n\nPROCESSING CLUSTER: ".$cluster."\n";

        my $full_logo_directory = $output_tmp_dir.$logo_directory.$cluster.'/';
        if (!-e $full_logo_directory) {
        	mkdir ($full_logo_directory) or die 'Error creating output folders!';
        }
        
        my $full_meme_directory = $meme_directory.$cluster.'/';
        if (!-e $full_meme_directory) {
        	mkdir ($full_meme_directory) or die 'Error creating output folders!';
        }
  
        my @geneIDs = @{ $clusters_hash{$cluster} }; # get gene ids for that cluster from %clusters
        my $number_of_IDs = 1+$#geneIDs;
        print "Length of gene-ID list is: ".$number_of_IDs."\n";
	  	
	  	print "Writing meme input file\n" ;
	  	my $meme_input_fasta_file = $full_meme_directory."meme_input.fasta";
	  	make_gi_set(\%sequences_from_db, \@geneIDs, $core_promoter_max_length, $meme_input_fasta_file);
			
		print "Running MEME now...\n\n";
      
	  	my $cmd = $meme_binary." ".$meme_input_fasta_file." ".$meme_inputs." -oc ".$full_meme_directory;
	  	system(`$cmd`);
	  	
	  	
    	#elements are array refs of hashes, one for each motif found by meme
						
		$all_meme_results{$cluster} =  parse_MEME_text_output($full_meme_directory."meme.txt");	
			
}

#MEME done, so output results

### write a file "summary.out" in output_tmp_dir and fill with the output below:

my $summary = $output_tmp_dir."summary.out"; # full path to file
open (SUMMARY, ">$summary") or die "cannot write a summary file at $summary\n";

  
print "Generating results files...\n";


#PEB - Variables created to hold values that will be written to javascript in results page
my @jsClusterID;
my @jsNumMotifs;
my @jsMotifID;
my @jsNumSites;
my @jsPCSites;
my @jsWidth;
my @jsInfoContent;
my @jsEValue;
my @jsPBias;
my @jsSBias;
my @jsWhichGeneIDs = ();
#These should be 2/3 D arrays
##########################

my $counter = 0;  


foreach my $cluster (@clusters) {

	
    print SUMMARY "CLUSTER ".$cluster."\n";
    $jsClusterID[$counter] = $cluster;       #cluster name
       
    my $motif_num = 0; 
    
    my @all_motifs = @{$all_meme_results{$cluster}};  #This is now an array of hashes
    
    foreach my $motif (@all_motifs) {
    	
				my $wm = $motif->{WM};
                # add simple pseudo-counts "by hand", attempt to get around problem with weight matrix distance function DELETE LATER
                my $wm_length = $wm->get_pssm_length();
                foreach (my $i=0; $i<$wm_length;$i++) {
                	foreach (my $nucleotide=0; $nucleotide<4; $nucleotide++) {
                                my $previous_freq = ${${$wm->wm_freq_dist}[$i]}[$nucleotide];
                                my $new_freq = ($previous_freq+0.01)/1.04;
                                ${${$wm->wm_freq_dist}[$i]}[$nucleotide] = $new_freq;
                	}
                }
                my $previous_ID = $wm->wm_identifier;       

                $wm->wm_identifier($clusters[$counter].'-'.$previous_ID);
                my $full_logo_directory = $output_tmp_dir.$logo_directory.$clusters[$counter].'/';
                my $filepath = $wm->generate_logo($full_logo_directory);
                print "\nMotif ".$previous_ID.":\n";
                print SUMMARY "\nMotif ".$previous_ID.":\n";
                print "--------\n";
                print SUMMARY "--------\n";
                print "Has ".$motif->{num_sites}." sites (";
                print SUMMARY "Has ".$motif->{num_sites}." sites (";
                my $percentage = 100*$motif->{occurrence_ratio};
                print $percentage." percent)\n";
                print SUMMARY $percentage." percent)\n";
                my @found_in_genes = @{$motif->{which_genes}};
                if(@found_in_genes){            	

                	for (my $g=0; $g <= $#found_in_genes; $g++){
                        if(!($g % 2)){
                        	print SUMMARY "\n";
                        }
                        print SUMMARY "Gene id ".$found_in_genes[$g]{'id'}.", strand ".$found_in_genes[$g]{'strand'}.", position ".$found_in_genes[$g]{'position'}."/".$found_in_genes[$g]{'promoter_length'}."\t";
                	}
                	print SUMMARY "\n\n";
                }

                print "Width: ".$motif->{motif_width}."\n";
                print SUMMARY "Width: ".$motif->{motif_width}."\n";
                my $info_content = $wm->get_information_content();

                print "Information Content: ".$info_content."\n";
                print SUMMARY "Information Content: ".$info_content."\n";
                print "E-value: ".$motif->{e_value}."\n";
                print SUMMARY "E-value: ".$motif->{e_value}."\n";

            	print "Positional bias: ".$motif->{positional_bias_pvalue}."\n";
            	print SUMMARY "Positional bias: ".$motif->{positional_bias_pvalue}."\n";

                print "Strand bias: ".$motif->{strand_bias_pvalue}."\n\n";
                print SUMMARY "Strand bias: ".$motif->{strand_bias_pvalue}."\n\n";
   
                #PEB added this to save the results in arrays that can be written to the javascript in the output html file

                $jsMotifID[$counter][$motif_num] = $previous_ID;
                $jsNumSites[$counter][$motif_num] = $motif->{num_sites};
                $jsPCSites[$counter][$motif_num] = $percentage;
                $jsWidth[$counter][$motif_num] = $motif->{motif_width};
                $jsInfoContent[$counter][$motif_num] = $info_content;
                $jsEValue[$counter][$motif_num] = $motif->{e_value};
                $jsPBias[$counter][$motif_num] = $motif->{positional_bias_pvalue};
                $jsSBias[$counter][$motif_num] = $motif->{strand_bias_pvalue};               
                $jsWhichGeneIDs[$counter][$motif_num] = [@found_in_genes];


                $motif_num++;

    }#end of foreach motif
    $jsNumMotifs[$counter] = $motif_num;
    $counter++;  
    

}#end of foreach $cluster

print "\n\n\n";
print SUMMARY "\n";
close (SUMMARY);


#PEB added this bit to generate an html output file

print "Generating html...\n";

#put input params into display format
$input_file = substr($input_file, rindex($input_file, "/")+1, length($input_file)-rindex($input_file, "/"));
my $sequence_location = substr($sequence_file, rindex($sequence_file, "/")+1, length($sequence_file)-rindex($sequence_file, "/"));
if ($background_file){
	$background_file = substr($background_file, rindex($background_file, "/")+1, length($background_file)-rindex($background_file, "/"));
}

$num_clusters = $counter;

#process gene id lists

my $htmlfile = $output_tmp_dir."results.html"; # full path to file
open (HTML, ">$htmlfile") or die "cannot write a results file at $htmlfile\n";

print HTML << "RESULTS" ;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>

<head>

<title>
        MEME-LaB Results
</title>

<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.5/jquery.min.js"></script>

<script type = 'text/javascript'>

        var numClusters = $num_clusters;
        var numMotifs = $num_motifs;

        \$(document).ready(function(){
                window.onresize = function() {sizeControls()}
                sizeControls();
                clearfilter();
        });

        function filter()
        {
                var numVisible = new Array(numClusters);

                var filterNumSites = document.getElementById('numSites').value;
                var filterPercentSites = document.getElementById('percentSites').value;
                var filterInfoContent = document.getElementById('infoContent').value;
                var filterPBias = document.getElementById('pBias').value;
                var filterSBias = document.getElementById('sBias').value;


                //remove any white space at start or end
                filterNumSites = filterNumSites.replace(/^\\s+|\\s+\$/g, '');
                if(filterNumSites.length == 0){
                        //nothing entered so set to low value, eliminating it from filter
                        filterNumSites = 0;
                }else{
                        filterNumSites = parseInt(filterNumSites);
                        if (isNaN(filterNumSites)){
                                //something, but not a number
                                alert('Please enter a valid number of sites.');
                                return;
                        }
                }
                filterPercentSites = filterPercentSites.replace(/^\\s+|\\s+\$/g, '');
                if(filterPercentSites.length == 0){
                        //nothing entered so set to low value, eliminating it from filter
                        filterPercentSites = 0;
                }else{
                        filterPercentSites = parseFloat(filterPercentSites);
                        if (isNaN(filterPercentSites) || (filterPercentSites > 100)){
                                //something, but not a number
                                alert('Please enter a valid percentage of sites.');
                                return;
                        }
                }
                filterInfoContent = filterInfoContent.replace(/^\\s+|\\s+\$/g, '');
                if(filterInfoContent.length == 0){
                        //nothing entered so set to low value, eliminating it from filter
                        filterInfoContent = 0;
                }else{
                        filterInfoContent = parseFloat(filterInfoContent);
                        if (isNaN(filterInfoContent)){
                                //something, but not a number
                                alert('Please enter a valid information content.');
                                return;
                        }
                }


                filterPBias = filterPBias.replace(/^\\s+|\\s+\$/g, '');
                if(filterPBias.length == 0){
                        //nothing entered so set to one, eliminating it from filter
                        filterPBias = 1;
                }else{
                        filterPBias = parseFloat(filterPBias);
                        if (isNaN(filterPBias) || (filterPBias > 1) || (filterPBias < 0)){
                                //something, but not a number
                                alert('Please enter a valid p value for positional bias (0 <= p <= 1).');
                                return;
                        }
                }


                filterSBias = filterSBias.replace(/^\\s+|\\s+\$/g, '');
                if(filterSBias.length == 0){
                        //nothing entered so set to one, eliminating it from filter
                        filterSBias = 1;
                }else{
                        filterSBias = parseFloat(filterSBias);
                        if (isNaN(filterSBias) || (filterSBias > 1) || (filterSBias < 0)){
                                //something, but not a number
                                alert('Please enter a valid p value for strand bias (0 <= p <= 1).');
                                return;
                        }
                }



                \$("[name='cluster']").each(function(i){
                        \$(this).data('motifs', 0);
                });

                var filterOn = false;

                \$("[name='motif']").each(function(i){
                        var showMotif = true;
                        var sites = \$(this).find("[name='sites']").html();
                        var result = /(\\d+)\\s*\\(([0-9\\.]+)\\s*%\\)/.exec(sites);
                        try{
                                var numSites = parseInt(result[1]);
                                var percentSites = parseFloat(result[2]);
                                showMotif = showMotif && (numSites >= filterNumSites) && (percentSites >= filterPercentSites);
                        }catch(err){
                                //ignore these params
                        }
                        if(showMotif == true){
                                //ok so far
                                var infoContent = parseFloat(\$(this).find("[name='infco']").html());
                                showMotif = showMotif && (infoContent >= filterInfoContent);
                                if(showMotif == true){
                                        var strandBias = parseFloat(\$(this).find("[name='strb']").html());
                                        showMotif = showMotif && (strandBias <= filterSBias);
                                        if(showMotif == true){
                                                var posBias = parseFloat(\$(this).find("[name='posb']").html());
                                                showMotif = showMotif && (posBias <= filterPBias);
                                        }

                                }
                        }

                        if(showMotif == true){
                                \$(this).css('display', 'inline');
                                var num = \$(this).parents("[name='cluster']").data('motifs');
                                \$(this).parents("[name='cluster']").data('motifs', num+1);
                        }else{
                                \$(this).css('display', 'none');
                                filterOn = true;
                        }

                });

                //no motifs, hide cluster
                var total = 0;
                \$("[name='cluster']").each(function(i){

                        var motifs = \$(this).data('motifs');
                        if(motifs == 0){
                                \$(this).css('display', 'none');
                        }else{
                                \$(this).css('display', 'inline');
                                total += motifs;
                        }

                });


                //update select to reflect number visible
                var ctrl = document.getElementById('select');
                while(ctrl.length){
                        ctrl.remove(0);
                }


                \$("[name='cluster']").each(function(i){

                        var motifs = \$(this).data('motifs');
                        if (motifs > 0){
                                var opt = document.createElement("option");
                                var name = \$(this).attr('id');
                                name = name.substr(0, 7)+" "+name.substr(7);
                                opt.text = name+" ("+motifs+")";

                                try{
                                        //for IE7
                                        ctrl.add(opt, ctrl.options[null]);
                                }catch (e){
                                        ctrl.add(opt, null);
                                }
                        } 

                });


                if (filterOn == true){
                        document.getElementById('title').innerHTML = "Showing "+total+" motifs (filtered)";
                }else{
                        document.getElementById('title').innerHTML = "Showing all "+total+" motifs";
                }
                if(total == 0){
                        document.getElementById('warning').style.display = 'block';
                }else{
                        document.getElementById('warning').style.display = 'none';
                }
        }

        function clearfilter(){

                document.getElementById('numSites').value = "";
                document.getElementById('percentSites').value = "";
                document.getElementById('infoContent').value = "";
                document.getElementById('pBias').value = "";
                document.getElementById('sBias').value = "";
          
                filter();
        }

        function jumpto(val){

                var cluster = /cluster\\s*(\\S+)\\s+\\(\\d+\\)/.exec(val);
                cluster = document.getElementById('cluster'+cluster[1]);
                window.scrollTo(0, cluster.offsetTop);
        }

        function sizeControls(){

                var leftcol = document.getElementById('leftcol');
                var height = (document.documentElement.clientHeight - leftcol.offsetTop) + "px";
                leftcol.style.height = height;
                var centrecol = document.getElementById('centrecol');
                var width = (document.documentElement.clientWidth - centrecol.offsetLeft - 20) + "px";
                centrecol.style.width = width;
        }

        function showInfo(n){

                var msg;

                if (n == 1){
                        msg = 'The information content measures the overall specificity of the motif. It is larger when nucleotides in the motif are more determined, or when the motif is longer.';
                }else if (n == 2){
                        msg = 'The Kolmogorov-Smirnov test is used to compare the observed distribution of motif start sites with an expected linear distribution along the length of the promoter. A strong p-value indicates a non-uniform distribution, providing additional evidence that the motif is linked to transcription.';
                }else if (n == 3){
                        msg = 'This statistic measures whether a motif occurs more frequently on one strand that the other. A strong p-value indicates that the motif has a strand preference, providing additional evidence that the motif is not a randomly occurring feature.';
                }

                alert(msg);

        }

        function showWhichGeneIDs(e, motif, title){

                if (!e) {
                                var e = window.event;
                }

                var top = e.screenY;
                var left = e.screenX-800;

                try{
                        window.open("$gene_list_path" + motif + ".html", title, "top="+top+",left="+left+",menubar=no,status=no,toolbar=no,width=800,scrollbars=yes");
                }catch(err){
                        window.open("$gene_list_path" + motif + ".html", title, "menubar=no,status=no,toolbar=no,width=800,scrollbars=yes");
                }
        }


</script>

<style type = "text/css">

        body {font-family:sans-serif;
                        background-color:white
                }

        .title {font-size:24px;
                        font-weight:bold;
                }

        .paramTable {margin-left:auto;
                                margin-right:auto;
                                padding-left:10px;
                                padding-right:10px;
                                font-size:16px;
                        }

        .paramTable p {font-weight:bold;}

        .paramTable table {white-space:nowrap;
                                                border:0px none white}

        .paramTable td {padding:2px 20px 2px 2px;}

        .findBox {position:fixed;
                                left:0px;
                                top:225px;
                                width:225px;
                                padding-left:10px;
                                padding-right:10px;
                                font-size:14px;
                                }

        .leftCol  {position:fixed;
                                left:0px;
                                top:295px;
                                width:225px;
                                padding-left:10px;
                                padding-right:10px;
                                overflow-y:auto;
                                font-size:14px;
                                }

        .leftCol div {position:absolute;
                                        text-align:right;}

        .leftCol input {text-align:right;
                                        width:40px;
                                        }

        .popupLink {cursor:pointer;
                                color:blue;
                                text-decoration:underline;
                                }

        .warning {font-size:18px;
                                margin-left:5px;
                                margin-right:5px;
                                color:red;
                                text-align:left;
                        }

        .centreCol {position:absolute;
                                left:245px;
                                top:250px;
                                padding:10px;
                        }

        .clusterBox {background-color:#CCD8DD;
                                border:1px solid grey;
                                padding:10px;
                                }

        .clusterTitle {
                                        font-size:16px;
                                        font-weight:bold;
                                }

        .motifBox {
                                        background-color:white;
                                        border:1px solid grey;
                                        padding:5px;
                                }

        .motifTable {
                                        border:0px none white;
                                }

        .motifTable td {
                                                vertical-align:middle;
                                                border:0px none white;
                                        }

        .motifDataTable {
                                                border:1px solid black;
                                                margin:10px;
                                                border-collapse:collapse
                                        }

        .motifDataTable th {
                                                background-color:#CCD8DD;
                                                font-weight:bold;
                                                font-size:14px;
                                                padding:5px;
                                                border:1px solid black;
                                                }

        .motifDataTable td {
                                                        padding:2px;
                                                        border-style:none;

                                                }

</style>

</head>

<body>

<p align='center' class = 'title'>MEME-LaB Results</p>

<div class = 'paramTable' align='center'><!-- showns input params -->

        <p>This analysis was run with the following input parameters</p>

        <table><colgroup span='4' width='0*'></colgroup>
                <tr>
                		<td style='text-align:right'>Gene list</td>
                        <td style='text-align:left'>$input_file</td>
                        <td style='text-align:right'>Sequences from</td>
                        <td style='text-align:left'>$sequence_location</td>
                </tr>
                <tr>    
                        <td style='text-align:right'>Promoter maximum length</td>
                        <td style='text-align:left'>$core_promoter_max_length</td>	
						<td style='text-align:right'>Number of motifs per cluster</td>
                        <td style='text-align:left'>$num_motifs</td>
                </tr>
                <tr>
                        <td style='text-align:right'>Minimum motif width</td>
                        <td style='text-align:left'>$min_motif_width</td>
                        <td style='text-align:right'>Maximum motif width</td>
                        <td style='text-align:left'>$max_motif_width</td>
                </tr>

	
RESULTS
	if ($background_file){
		print HTML << "RESULTS" ;

				<tr>
					<td style='text-align:right'>Background model</td>
					<td style='text-align:left'>$background_file</td>
				</tr>
		
RESULTS
	
	}
	
print HTML << "RESULTS" ;

        </table>
</div>

<div class = 'findBox'>
        <p>
        <span id = 'title' style='font-style:italic'></span>
        </p>
        <p>
                <b>Find</b>
                <span style='float:right'>
                        <select id = 'select' onchange='jumpto(this.value)'>

                        </select>
                </span>
        </p>

</div>

<div id = 'leftcol' class = 'leftCol'>

        <p><b>Filter motifs</b></p>
        <p><i>You can filter the motifs displayed by the following criteria</i></p>
        <ul>
                <li>Minimum number of promoters that the motif appears in</li>
                <li>Minimum percentage of promoters that the motif appears in</li>
                <li>Minimum information content <span class = 'popupLink' onclick = 'showInfo(1)'>more...</span></li>
                <li>Maximum p-value for positional bias <span class = 'popupLink'  onclick = 'showInfo(2)'>more...</span></li>
                <li>Maximum p-value for strand bias <span class = 'popupLink' onclick = 'showInfo(3)'>more...</span></li>
                
        </ul>

        <div style = 'top:290px;left:0px;width:140px;'>
                Promoters &gt=
        </div>
        <div style = 'top:285px;left:145px'>
                <input type='text' id='numSites'/>
        </div>

        <div style = 'top:320px;left:0px;width:140px;'>
                % Promoters &gt=
        </div>
        <div style = 'top:315px;left:145px'>
                <input type='text' id='percentSites' />
        </div>

        <div style = 'top:350px;left:0px;width:140px;'>
                Info content &gt=
        </div>
        <div style = 'top:345px;left:145px'>
                <input type='text' id='infoContent'/>
        </div>

        <div style = 'top:380px;left:0px;width:140px;'>
                Positional bias &lt=
        </div>
        <div style = 'top:375px;left:145px'>
                <input type='text' id='pBias' />
        </div>

        <div style = 'top:410px;left:0px;width:140px;'>
                Strand bias &lt=
        </div>
        <div style = 'top:405px;left:145px'>
                <input type='text' id='sBias'/>
        </div>


        <div style = 'top:455px;left:0px;width:140px;'>
                <button type='button' onclick='clearfilter()'>Clear</button>
        </div>
        <div style = 'top:455px;left:140px;margin-left:5px;'>
                <button type='button' onclick='filter()'>Apply</button>
        </div>

        <br /><br />
</div>

<div class = 'centreCol' id = 'centrecol'><!-- clusters -->

    <span class = 'warning' id = 'warning' style='text-align:center;display:none;'>No motifs meet the specified criteria</span>

RESULTS


for($c = 0; $c <= $#jsClusterID; $c++){

        my $cluster_id = "cluster".$jsClusterID[$c];
        print HTML "\t<span id = '",$cluster_id,"' name = 'cluster'>\n";
        print HTML "\t\t<div class = 'clusterBox'>\n";
        my $cnum = $c+1;

        for ($m = 0; $m < $jsNumMotifs[$c]; $m++){

                my $motif_id = $cluster_id."_motif".$jsMotifID[$c][$m];
                my $motif_title = "Cluster  ".$jsClusterID[$c]." Motif ".$jsMotifID[$c][$m];
                print HTML "\t\t\t<span id = '",$motif_id, "' name = 'motif'>\n";
                print HTML "\t\t\t\t<p class = 'clusterTitle'>", $motif_title, "</p>\n";
                print HTML "\t\t\t\t<div class = 'motifBox'>\n";
                print HTML "\t\t\t\t\t<table class = 'motifTable'><tr>\n";

                #relative path
                my $logo_file = $logo_directory.$jsClusterID[$c].'/'.$jsClusterID[$c].'-'.$jsMotifID[$c][$m].'.gif' ;  #Was .png

                print HTML "\t\t\t\t\t\t<td><img src='$logo_file' id='", $motif_id, "_image'/></td>\n";
                print HTML "\t\t\t\t\t\t<td>\n";
                print HTML "\t\t\t\t\t\t\t<table class='motifDataTable'>\n";
                print HTML "\t\t\t\t\t\t\t\t<tr><th colspan='2'>Motif Features";

                my @gene_list = @{$jsWhichGeneIDs[$c][$m]};
                if(@gene_list){
                        #create html files listing genes associated with each motif
                        write_which_gene_ids_file($output_tmp_dir.$gene_list_path, $motif_id, $motif_title, \@gene_list, $jsWidth[$c][$m],  sprintf("%.2g", $jsPBias[$c][$m]),  sprintf("%.2g", $jsSBias[$c][$m]));
                        print HTML "<span style='margin-left:10px' class = 'popupLink' onclick = 'showWhichGeneIDs(event, \"", $motif_id, "\", \"", $motif_title, "\")'>more...</span>";
                }

                print HTML "</th></tr>\n";

                print HTML "\t\t\t\t\t\t\t\t<tr><td>Sites</td><td id = '", $motif_id, "_sites' name='sites'>", $jsNumSites[$c][$m], " (", $jsPCSites[$c][$m], "%)</td></tr>\n";
                print HTML "\t\t\t\t\t\t\t\t<tr><td>Width</td><td id = '", $motif_id, "_width' name = 'width'>", $jsWidth[$c][$m], "</td></tr>\n";
                print HTML "\t\t\t\t\t\t\t\t<tr><td>Information content</td><td id = '", $motif_id, "_infco' name = 'infco'>"; printf HTML "%.1f", $jsInfoContent[$c][$m]; print HTML "</td></tr>\n";
                print HTML "\t\t\t\t\t\t\t\t<tr><td>E-value</td><td id = '", $motif_id, "_eval'>"; printf HTML "%.2g", $jsEValue[$c][$m]; print HTML "</td></tr>\n";
                print HTML "\t\t\t\t\t\t\t\t<tr><td>Positional bias</td><td id = '", $motif_id, "_posb' name = 'posb'>"; printf HTML "%.2g", $jsPBias[$c][$m]; print HTML "</td></tr>\n";
                print HTML "\t\t\t\t\t\t\t\t<tr><td>Strand bias</td><td id = '", $motif_id, "_strb' name = 'strb'>"; printf HTML "%.2g", $jsSBias[$c][$m]; print HTML "</td></tr>\n";
                print HTML "\t\t\t\t\t\t\t</table>\n";
                print HTML "\t\t\t\t\t\t</td>\n";
                print HTML "\t\t\t\t\t</tr></table>\n";
                print HTML "\t\t\t\t</div>\n";
                print HTML "\t\t\t</span>\n";
                
              
        }
        print HTML "\t\t</div>\n";
        print HTML "\t<br/>\n";
        print HTML "\t</span>\n";

}

print HTML << "RESULTS" ;

</div>

</body>

</html>

RESULTS

close(HTML);

print "Done\n\n";
    
#finished


sub write_which_gene_ids_file {

        #creates an html file that lists all the gene ids associated with a specific motif

        my $file_path = $_[0];          #directory to save file in
        my $id = $_[1];                         #name for file (no spaces) 'clusterX_motifY'
        my $title = $_[2];                      #title for file 'Cluster X Motif Y'
        my @gene_list = @{$_[3]};       #list of hashes with fields 'id', 'strand' and 'position'
        my $width = $_[4];                      #motif width, needed to calculate last possible start position
        my $pos_bias = $_[5];           #positional bias p value
        my $strand_bias = $_[6];        #strand bias p value

        if (substr($file_path, -1, 1) ne "/"){
                $file_path = $file_path."/";
        }

        my $htmlfile = $file_path.$id.".html";
        open (FH, ">$htmlfile") or die "cannot write a results file at $htmlfile\n";
        print FH << "RESULTS" ;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">

<html>
<head>
        <title>$title</title>

        <script type="text/javascript">
                function sizeWindow(){

						var totalWidth = 0;
                        for(var col = 1; col<=3; col++){
                                var colWidth = document.getElementById('c'+col).clientWidth;
                                var titleWidth = document.getElementById('t'+col).clientWidth;
                                var width = Math.max(colWidth, titleWidth);
                                document.getElementById('t'+col).style.width = width+"px";
                                document.getElementById('c'+col).style.width = width+"px";
                                totalWidth += width;
                        }
                        
                        document.getElementById('note').style.width = (totalWidth+30)+"px";


                        var tbl = document.getElementById('tbl');
                        var width = tbl.offsetWidth+20;
                        document.getElementById('scrollBox').style.width = width+"px";
                        width+=25;
                        document.getElementById('distCol').style.left = width+"px";

                        window.resizeTo(width + 560, document.getElementById('distCol').clientHeight+document.getElementById('distCol').offsetTop+60);
                        document.getElementById('leftCol').style.visibility='visible';
                        document.getElementById('distCol').style.visibility='visible';
                }

                window.onresize = function(){
                        var height = (document.documentElement)? document.documentElement.clientHeight : document.clientHeight;
                        document.getElementById('scrollBox').style.height = (height - document.getElementById('leftCol').offsetTop-100)+"px";
                }
                
                function highlightMarker(id, show){
                	
					if (show == true){
						document.getElementById(id).style.backgroundColor="red";
						document.getElementById(id).style.height = "30px";
						document.getElementById(id).style.top = "-10px";
					}else{
						document.getElementById(id).style.backgroundColor="black";
						document.getElementById(id).style.height = "10px";
						document.getElementById(id).style.top = "0px";

					}                	
                	
                }


        </script>

        <style type = "text/css">

                body    {
                                font-family:sans-serif;
                                background-color:white;
                                padding:10px;
                        }

                table   {
                                border-width:0px;
                                border-style:none;
                                font-size:14px;
                                border-collapse:collapse;
                                table-layout:fixed;
                                white-space:nowrap;
                        }

                th      {
                                padding:5px 0px 5px 5px;
                                vertical-align:top;
                }

                td      {
                                padding:5px 0px 5px 5px;
                                vertical-align:top;
                }

                .strand {    
                                position:relative;
                                margin-top:15px;
                                background-color:grey;
                                height:10px;
                                width:520px;
                }

                .strandLabel {
                                position:absolute;
                                background-color:white;
                                width:20px;
                                height:10px;
                                left:0px;
                                top:0px;
                                margin:0px;
                }

                .strandLabel span {
                        font-size:16px;
                        position:relative;
                        top:-4px;
                }

                .marker {
                                position:absolute;
                                background-color:black;
                                height:10px;
                                top:0px;
                                margin:0px;
                                width:5px;
                }
                .leftCol {

                        left:0px;
                        top:75px;
                        position:fixed;
                        padding-left:20px;
                        visibility:hidden;

                }

                .distCol {
                        top:75px;
                        width:550px;
                        position:fixed;
                        visibility:hidden;
                }

                .scrollBox {
                        overflow-y:auto;
                        overflow-x:hidden;
                }


        </style>

</head>

<body onload = "setTimeout(function(){sizeWindow()}, 100)">
        <div style='font-size:18px;background-color:#CCD8DD;padding:10px;'>
                $title is found in the following promoters
        </div>
        <div id='leftCol' class='leftCol'>
                <table id = "tbl">
                        <tr>
                                <td style='padding:0px'>
                                        <table>
                                                <tr>
                                                        <th id='t1'>Gene ID</th><th id='t2'>Strand</th><th id='t3'>Position*</th>
                                                </tr>
                                        </table>
                                </td>
                        </tr>
                        <tr>
                                <td style='padding:0px'>
                                        <div class='scrollBox' id='scrollBox'>
                                                <table>

RESULTS
        my $num_genes = scalar (@gene_list);

        for($g=0; $g<= $#gene_list; $g++){

                print FH "\t\t\t\t\t\t\t<tr onmouseover='highlightMarker(\"", $gene_list[$g]{'id'}, "\", true)' onmouseout='highlightMarker(\"", $gene_list[$g]{'id'}, "\", false)'>\n";
                if(!$g){
                        print FH "\t\t\t\t\t\t\t\t<td id='c1'>", $gene_list[$g]{'id'},"</td><td align='center' id='c2'>", $gene_list[$g]{'strand'},"</td><td align='right' id='c3'>", $gene_list[$g]{'position'}, "/", $gene_list[$g]{'promoter_length'}, "</td>\n";
                }else{
                        print FH "\t\t\t\t\t\t\t\t<td>", $gene_list[$g]{'id'},"</td><td align='center'>", $gene_list[$g]{'strand'},"</td><td align='right'>", $gene_list[$g]{'position'}, "/", $gene_list[$g]{'promoter_length'}, "</td>\n";
                }
                print FH "\t\t\t\t\t\t\t</tr>\n";
        }

        print FH << "RESULTS" ;
                                                </table>
                                        </div>
                                </td>
                        </tr>
                        <tr>
                        	<td>
                        		<div id = 'note' style='font-style:italic;font-size:12px;padding-top:10px;white-space:normal'>
                        			<b>*</b> Overlapping sequences are merged to avoid redundancy. As a result some sequences in the final set may be longer than the specified maximum.
                        		</div>
                        	</td>
                        </tr>
                </table>
        </div>
        <div id='distCol' class='distCol'>
                <table>
                        <tr>
                                <th>Distribution</th>
                        </tr>
                        <tr>
                                <td>
                                        <div class='strand'>
                                                <div class = 'strandLabel'><span>+</span></div>
RESULTS
        my $pos_count = 0;
        for($i=0; $i< $num_genes; $i++){
                if($gene_list[$i]{'strand'} eq '+'){
                        my $pos = $gene_list[$i]{'position'}/($gene_list[$i]{'promoter_length'}-$width+1);
                        my $left = (20+int($pos*(500-$width+1)))."px";
                        print FH "\t\t\t\t\t\t<div id = '", $gene_list[$i]{'id'}, "' class='marker' style='left:$left;'></div>\n";
                        $pos_count++;
                }
        }
        print FH << "RESULTS" ;
                                        </div>
                                </td>
                        </tr>
                        <tr>
                                <td>
                                        <div class='strand'>
                                                <div class = 'strandLabel'><span>-</span></div>
RESULTS
        my $neg_count = 0;
        for($i=0; $i< $num_genes; $i++){
                if($gene_list[$i]{'strand'} eq '-'){
                        my $pos = $gene_list[$i]{'position'}/($gene_list[$i]{'promoter_length'}-$width+1);
                        my $left = (20+int($pos*(500-$width+1)))."px";
                        print FH "\t\t\t\t\t\t<div  id = '", $gene_list[$i]{'id'}, "' class='marker' style='left:$left;'></div>\n";
                        $neg_count++;
                }
        }
        print FH << "RESULTS" ;
                                        </div>
                                </td>
                        </tr>
                        <tr>
                                <td>
                                        <div style = 'font-style:italic;margin-top:20px;margin-left:20px;white-space:normal'>
                                                <p>Position values represent number of bases from the 5 prime end of the sequence provided to MEME and the length of that sequence. It is these figures expressed as proportions that are used for statistical analysis.</p>
                                                <p>Probability value of comparison of these values to a linear distribution, p = $pos_bias</p>
                                                <p>The motif is found $pos_count times on the + stand and $neg_count times on the - strand.</p>
                                                <p>Probability value that it is evenly distributed between the two strands, p = $strand_bias</p>
                                        		<p>Move the mouse pointer over the gene ids to highlight the corresponding marker</p>
                                        </div>
                                </td>
                        </tr>
                </table>
        </div>
</body>
</html>

RESULTS


        close(FH);

}

sub valid_clusters {
	
	use Bio::SeqIO;

	### parsing input: read in file of gene IDs
	### check validity of listed gene IDs
	### remove any invalid gene IDs
	### if more than half of gene IDs overall are invalid, report an error "check file format" and don't run
	### if all IDs in one or more cluster(s) are invalid, report an error and don't ru

	my $input = $_[0]; # must be tab separated with no trailing white space
	my $output = $input;  # overwrites input file (with invalid IDs removed)
	my $fasta_file = $_[1];

	my @verified_ids_array;
	my %all_sequences; #This will be returned

	#read all gene ids in fasta file
	my $seqfile = Bio::SeqIO->new(-file => $fasta_file, -format => "fasta");	
	my $count = 0;
	while($seq_obj = $seqfile->next_seq){
		my $gene_id = $seq_obj->display_id;
		if($seq_obj->seq ne ""){
			push(@verified_ids_array, $gene_id);	
	  	   $all_sequences{$gene_id}=$seq_obj->seq;
				
		}
		if (++$count%500 == 0){
			print ".";
		}
	}
	
	my $ns = scalar (@verified_ids_array);
	print "\nFound ".$ns." sequences in FASTA file\n";	


	my %verified_ids; #hash keys will be ids in fasta file
	foreach my $id (@verified_ids_array) {
		$verified_ids{$id} = 13; #placeholder
	}


	# open and parse input file
	
	open(INPUTFILE, '<'.$input) or die "Cannot open '$input': $!";
	my %clusters = (); # a hash of clusters with key being cluster number and values being an array of gene ids in the cluster
	
	my %seen;
	my @clusters_ids_in_input_order = (); #A list of unique names
	
	while (<INPUTFILE>) {
 		chomp;
 		my ($cluster_id, $gene_to_add) = split("\t");
 		if (defined $gene_to_add && defined $cluster_id)  {
     		push @ { $clusters{$cluster_id} }, $gene_to_add;
     		push (@cluster_ids_in_input_order, $cluster_id) unless ($seen{$cluster_id}++);
 		}
	}
	close INPUTFILE;

	# check IDs in each cluster
	my $number_of_clusters = scalar @cluster_ids_in_input_order;
	print "number of clusters: ".$number_of_clusters."\n";

	my $total_validated_gene_ids = 0;
	my $total_good_ids = 0;
	my $there_was_a_rejected_gene_id_without_non_whitespace_character = FALSE;
	my @voided_clusters_empty;
	my @rejected_gene_ids;

	# First remove any invalid ids
	my %clusters_valid_ids_only = ();
	my @cluster_ids_in_input_order_valid_only = ();

	foreach my $clustername (@cluster_ids_in_input_order) {
  		print "checking IDs of cluster ". $clustername."\n";
  		my @ids_to_check = @{ $clusters{$clustername} };
  		my $number_of_ids_in_cluster = scalar @ids_to_check;
  		print "IDs to check: ".$number_of_ids_in_cluster."\n";
		my @tmp_cluster = ();
  
  		foreach my $element(@ids_to_check) {
			# for each id in cluster
  			if (defined $verified_ids{$element}) {
      			push (@tmp_cluster, $element); #id is in fasta file
    		} else {
				if ($element =~ m/\S/) {  # only rejected IDs with non-whitespace characters are shown to the user
	    			push (@rejected_gene_ids, $element);
				} else {
	    			$there_was_a_rejected_gene_id_without_non_whitespace_character = TRUE;
				}
    		}
  		}
  		
  		my $valid_cluster_size = scalar @tmp_cluster;
 		# if a cluster contains no valid IDs or just 1, push it into the voided cluster array which will be reported to the user
  		if ($valid_cluster_size < 2) {
    		push (@voided_clusters_empty, $clustername);
    		print "Invalid cluster\n";
  		}else{
  			#valid cluster so save all valid ids
  			$clusters_valid_ids_only{$clustername} = [@tmp_cluster];
  			push(@cluster_ids_in_input_order_valid_only, $clustername);
  			print "Valid cluster\n";
  		}
  	
  		print $valid_cluster_size." elements exist\n";
  		$total_validated_gene_ids = $total_validated_gene_ids + $number_of_ids_in_cluster;
  		print "running total of gene IDs that have been checked is: ". $total_validated_gene_ids."\n";
  		$total_good_ids = $total_good_ids + $valid_cluster_size;
  		print "running total of correct gene IDs is: ".$total_good_ids."\n\n";
	}

	# Now know all we need to know to decide how to proceed
	# if any cluster has no valid IDs, report an error and do not run

    if (scalar @voided_clusters_empty > 0) {
		print "ERROR: input file contains the following clusters with fewer than 2 valid gene-IDs: \n"; 
		foreach my $name (@voided_clusters_empty) {
	    	print $name."\n";
		}
		print "Please check the input file.\n";
		return (0, \@cluster_ids_in_input_order_valid_only, \%clusters_valid_ids_only, \%all_sequences); # return status: FAIL
    }
    
	# Check great enough proportion of all ids valid

	my $ratio_for_output = 0;
	if ($total_validated_gene_ids>0) {
    	$ratio_for_output = $total_good_ids/$total_validated_gene_ids;
	}
	print "good IDs: ".$total_good_ids." / all IDs: ".$total_validated_gene_ids." = ". $ratio_for_output * 100 ."%\n";

	if ($ratio_for_output < 0.5) {
  		# if overall % valid IDs is less than 50%, report an error and do not run
  		print "ERROR: more than 50% of the gene IDs are invalid. Please check your input file.\n";
  		return (0, \@cluster_ids_in_input_order_valid_only, \%clusters_valid_ids_only, \%all_sequences); # return status: FAIL
	}


	# Now decide if warnings are required

	if ($ratio_for_output < 1){
		# Warn about invalid ids

		print "WARNING: some gene-IDs were not found in the fasta file and were removed.\n";
		my $formatted_ratio = int(10000*$ratio_for_output)/100;
		print $formatted_ratio . "% of gene-IDs remain.\n";
		my $number_of_rejected_ids = scalar @rejected_gene_ids;
		
		if ($number_of_rejected_ids>0) {
			my $number_to_print = 10;
			my $printing_all = FALSE;
			if ($number_of_rejected_ids<=$number_to_print) {
				$number_to_print = $number_of_rejected_ids;
				$printing_all = TRUE;
			}
			if ($printing_all) {
				print("These gene-IDs were not found:\n");
			} else {
				print("Some examples of gene-IDs that were not found:\n");
			}

			for (my $x=0;$x<$number_to_print;$x++) {
				print $rejected_gene_ids[$x]."\n";
			}
			if ($there_was_a_rejected_gene_id_without_non_whitespace_character) {
				print "There was also at least one gene-ID that was made of whitespace characters only.\n";
			}
		} else {
			print "All invalid gene-IDs were made of whitespace characters only.\n";
		}

	}else{
		print "All gene IDs are valid.\n";
	}
	
	# create output file form valid ids and clusters

	# open output filehandle (for the cleaned-up input file). 
	# (THIS WILL OVERWRITE THE ORIGINAL INPUT FILE).


	open FILE, ">", $output or die $!;

	foreach my $clustername (@cluster_ids_in_input_order_valid_only) {
 
  		my @gene_ids = @{ $clusters_valid_ids_only{$clustername} };
  
  		foreach my $element(@gene_ids) {
	  		print FILE $clustername."\t".$element."\n";
      	}
  	}

	close FILE;
	
	#return flag to indicate success and an array of unique cluster names in same order as input file, and 
	#a hash of gene lists with cluster names as keys. Only vlaid clusters with  > 1 vlaid gene id are included
	return (1, \@cluster_ids_in_input_order_valid_only, \%clusters_valid_ids_only, \%all_sequences);
	
}


sub make_gi_set {
	
	
 my %sequences = %{$_[0]};
 my @ids = @{$_[1]};
 my $length = $_[2];
 my $fasta_file = $_[3];
 

my $n = 0;
open FH, ">".$fasta_file;	
	 
# take specified length of sequence, unless this is longer than the database sequence

foreach my $gene_id(@ids){
		
	$sequence = $sequences{$gene_id};
	my $length_to_take = $length;
    my $db_seq_length = length($sequence);
    if ($db_seq_length < $length) {
       	$length_to_take = $db_seq_length;
    }
	
    my $subseq = substr $sequence, -$length_to_take;
		
	print "Retrieved gene id ".$gene_id." from the FASTA database\n";
    print $subseq."\n";
      	
   	print FH ">".$gene_id."\n";
	print FH $subseq."\n";
	$n = $n+1;
  }
  
 close FH;	
 print "File created at ".$fasta_file." with ".$n." ids\n";

} # make_gi_set #

 sub parse_MEME_text_output {	
 	
 	my $inputFile = $_[0];
	
	open (FILEHANDLE, "< " . $inputFile) or die "oops: Cannot open $inputFile $!";
    my $lineNumber = 0;
    my $nmotifs = 0;
    my $numSeqs;
    my %promotor_lengths = ();
    my $longest_promoter = 0;
    my @motif_data;
    #my $cn = substr($d, -1,1);
    #$GU->user_info( 3, "\nCLUSTER $cn\n" ); # where is $d defined?
    while (<FILEHANDLE>) {
    	my $lineOfInput = $_;
	
      	$lineNumber++;
      	
      	if ($lineOfInput =~ m/^TRAINING SET/) {
      		
      		#get promotor lengths
  
      		 while (<FILEHANDLE>){
      		 	
      		 	my $lineOfInput = $_;
      		 	
      		 	if ($lineOfInput =~ m/\S+\s+[\d\.]+\s+\d+\s+/) {
      		 		
      		 		while($lineOfInput =~ m/(\S+)\s+[\d\.]+\s+(\d+)\s+/g) {
      		 		
      		 			$promotor_lengths{$1} = $2;
      		 			if($2 > $longest_promoter){
      		 				$longest_promoter = $2;
      		 			}
      		 			
      		 		}
      		 		
      		 	}
      		 	
      		 	if ($lineOfInput =~ m/^COMMAND/) {
      		 		last;
      		 	}
      		 	
      		 }
      	}
	
      	if ($lineOfInput =~ m/^model:/) {
		
			my @split = split(/\s+/,$_);
			$nmotifs = $split[4];
      	}
	
      	if ($lineOfInput =~ m/^data/) {
			
			my @split = split(/\s+/,$_);
			$numSeqs = $split[4];
			
			for my $m (1 .. $nmotifs) {
				
	 	 		my $rightCount = 0;
	  			my $leftCount = 0;
	  			my $plusStrandCount = 0;
	  			my $minusStrandCount = 0;
	  			#my $png;
				
	  			# Initialise variables for individual motif data
	  			my $motifNum;
	  			my $width;
	  			my $sites;
	  			my $evalue;
	  			my $ratio;
	  			my @matrix;
	  			my $pos_bias;
	  			my $strand_bias_pvalue;
	  			my @which_gene_ids = ();
	  
	  			while (<FILEHANDLE>) {
	    
	  				my $lineOfInput = $_;
	    
	    			if ($lineOfInput =~ m/^MOTIF/) {
	      
	      				my @split = split(/\s+/,$_);
	      
	      				$motifNum = $split[1];
	      				
	      				
	      				#webtool uses meme 4.0.0, line looks like this
	      				#MOTIF  1	width =    8   sites =   6   llr = 58   E-value = 1.1e+004
	      				
	      				#meme 4.11.2 has an extra word
	      				#MOTIF  1 MEME	width =   8  sites =  32  llr = 342  E-value = 8.7e-006
	      				
	      				#can get round this by counting backwards
	      
	      				#$width = $split[4];
	      				#$sites = $split[7];
	      				#$evalue = $split[13];
	      				
	      				$max_elem = $#split;
	      				
	      				$width = $split[$max_elem-9];
	      				$sites = $split[$max_elem-6];
	      				$evalue = $split[$max_elem];
	      				
	      				
	      				$ratio = $sites / $numSeqs;
	      				$ratio = substr($ratio, 0, 5);
	      				print "MOTIF $motifNum\tWidth: $width\tSites: $sites\tRatio: $ratio\tE-value:$evalue\n";
	      
	   				 }
	     
	   				 if ($lineOfInput =~ m/\s{2,}(\+|-)\s{2,}/) {  # Sascha: changed this as Richard's previous pattern ($lineOfInput =~ m/\s{10,}(\+|-)/))
                                                          #         did not work on fungal data where gene-IDs are too long to leave 10 white-space
                                                          #         characters
						print "****$_****" ;
	      				my @split = split(/\s+/,$_);
	      				my $strand = $split[1];
	     				my $position = $split[2];
	     				my $gene_id = $split[0];
	     				print $position."\n" ;
	     				my $scaled_position = 0;
	      				# check if numeric
	      				unless($position =~ m/[^0-9.]/) {
	      					#scale position as each promotor a different length
							$scaled_position = ($position/($promotor_lengths{$gene_id}-$width+1)) * ($longest_promoter-$width+1);
				#			if ($position >= ($promotor_lengths{$gene_id}/2)) {
		  		#				$rightCount++;
				#			} else {
		  		#				$leftCount++;
				#			}	    
						}
	      				# update strand counts
	      				if ($strand eq '+') {
		  					$plusStrandCount = $plusStrandCount + 1;
	      				} else {
		  					if ($strand eq '-') {
		      					$minusStrandCount = $minusStrandCount + 1;	  
		  					} else {
		      					print "Strand was neither plus nor minus! (it was: ".$strand.")\n";
		  					}
	      				}
	      				push(@which_gene_ids, {'id', $gene_id, 'strand', $strand, 'position', $position, 'scaled_position', $scaled_position, 'promoter_length', $promotor_lengths{$gene_id} });
	    			}
					
	    			if ($lineOfInput =~ m/Time/) {
	      				# Quick! Work out the binomial probabilty of getting positional distribution
	      				#my $p;
						my $statistics = Statistics->new();
						
						my @positions = ();
						for(my $g=0;$g<=$#which_gene_ids;$g++){
							unless($which_gene_ids[$g]{'position'} =~ m/[^0-9.]/){
								push(@positions, $which_gene_ids[$g]{'scaled_position'});
							}
						}
						#Use K-S test against a linear (random) distribution
						my @cdf = (1 .. ($longest_promoter-$width+1));	#theoretical linear distribution range
						$pos_bias = $statistics->one_sided_ks_test(\@positions, \@cdf, 'linear');
						
				#		$GU->user_info( 3, "Right count: ".$rightCount."\n" );
				#		$GU->user_info( 3, "Left count: ".$leftCount."\n" );
				#		if ($rightCount >= $leftCount) {
		    	#			$p = $statistics->get_binomial_probability($rightCount+$leftCount,$rightCount,0.5);
				#		} else {
		    	#			$p = $statistics->get_binomial_probability($rightCount+$leftCount,$leftCount,0.5);
				#		}
		
						# Compute strand bias p-value
						print "Positive strand count: ".$plusStrandCount."\n";
						print "Negative strand count: ".$minusStrandCount."\n";
						if ($plusStrandCount >= $minusStrandCount) {
		    				$strand_bias_pvalue = $statistics->get_binomial_probability($plusStrandCount+$minusStrandCount,$plusStrandCount,0.5);
						} else {
		    				$strand_bias_pvalue = $statistics->get_binomial_probability($plusStrandCount+$minusStrandCount,$minusStrandCount,0.5);
						}
	     	      						
						print "Prob. of distribution: ".$pos_bias."\n" ;
						print "Prob. of strand bias: ".$strand_bias_pvalue."\n" ;
						
	     	 			last;
						
	    			}
					
	    			if ($lineOfInput =~ m/letter-probability matrix/) {
							
	      				print "Hit Matrix Data\n" ;
	      				while (<FILEHANDLE>) {
							chomp;
							
							if ($_ =~ m/^-/) {
			  					last;	
							}
								
							my @bases = split(/\s+/, $_);
							shift(@bases);
							foreach (@bases) {
		  						print "...$_...\t" ;	
							}
							print "\n" ;
								
							push(@matrix, \@bases);
								
	      				}	
	    			}
				} #end while
				
				#Sort gene ids
				my @sorted_gene_ids = sort {$a->{'id'} cmp $b->{'id'}} @which_gene_ids;
				

	  			# Store all motif info in anonomous hash
	  			my $wm = Generic_Weight_Matrix->new(wm_identifier=>$motifNum, wm_freq_dist=>\@matrix);
	  			
	  			my $motif_record = {
			      	motif_width  => $width,
			      	num_sites  => $sites,
			      	occurrence_ratio  => $ratio,
			      	e_value => $evalue,
			      	positional_bias_pvalue => $pos_bias,
			      	WM => $wm,
			      	strand_bias_pvalue => $strand_bias_pvalue,
			      	which_genes => \@sorted_gene_ids,
				};
	  			print "$motifNum\n" ;
	  			push(@motif_data, $motif_record);
			}	#end for 
		} # end if	
	} # end while
	return \@motif_data;
} # parse_MEME_text_output #
	
	


