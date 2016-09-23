use MooseX::Declare;


class MEME_Parameters {
	
	use constant {FALSE => 0,
		  TRUE	=> 1};	
   
    has '-oc' => (is => 'rw', isa => 'Str'); # name of directory for output files will replace existing directory
    has '-dna' => (is => 'rw', isa => 'Int', default => TRUE); #sequences use DNA alphabet
    has '-mod' => (is => 'rw', isa => 'Str', default => 'zoops');		# [-mod oops|zoops|anr]		distribution of motifs
    has '-nmotifs' => (is => 'rw', isa => 'Int', default => '3');	# [-nmotifs <nmotifs>]		maximum number of motifs to find
    has '-evt' => (is => 'rw', isa => 'Int');		# [-evt <ev>]				stop if motif e-value greater than <evt>
    has '-nsites' => (is => 'rw', isa => 'Int');	# [-nsites <sites>]			number of sites for each motif
    has '-minsites' => (is => 'rw', isa => 'Int');	# [-minsites <minsites>]	minimum number of sites for each motif
    has '-maxsites' => (is => 'rw', isa => 'Int');	# [-maxsites <maxsites>]	maximum number of sites for each motif
    has '-w' => (is => 'rw', isa => 'Int');			# [-w <w>]					motif width
    has '-minw' => (is => 'rw', isa => 'Int');		# [-minw <minw>]			minumum motif width
    has '-maxw' => (is => 'rw', isa => 'Int');		# [-maxw <maxw>]			maximum motif width
    has '-revcomp' => (is => 'rw', isa => 'Int', default => TRUE); 	# [-revcomp]				allow sites on + or - dna strands
    has '-pal' => (is => 'rw', isa => 'Int', default => FALSE);		# [-pal]					force palindromes (requires -dna)
    has '-maxiter' => (is => 'rw', isa => 'Int', default => 5000);	# [-maxiter <maxiter>]		maximum em iterations to run
    has '-distance' => (is => 'rw', isa => 'Int');	# [-distance <distance>]	em convergence criterion
    has '-cons' => (is => 'rw', isa => 'Str'); 		# [-cons <cons>]			consensus sequence to start EM from
    has '-maxsize' => (is => 'rw', isa => 'Int', default => 600000); # must be above size of sequence file for MEME to accept input
	has '-bfile' => (is => 'rw', isa => 'Str');		# [-bfile <file> ] bckground model file
} # MEME_Parameters #


class Statistics {
	
	
	method get_binomial_probability(Int $n, Int $k, Num $p) {# Changed $p type form Int to Num
	    if ( $p == 0 ) {	# If the p-value is 0...
		if ( $k == 0 ) {	# AND the next position is 0...
		    return 1;		# Return 1.
		}
		else {				# But if the next site exists
		    return 0;		# Return 0
		}
	    }
	    if ( $p == 1 ) {	# If the p-value is 1
		return 1;			# Return 1
	    }		
	    my @firstsum;		# If the p-value NE 0 or 1, do stuff.
	    my @secondsum;
	    $#firstsum = $n;
	    $#secondsum = $n;
	    for ( my $i = 0 ; $i <= $n ; $i++ ) {
		if ( $i == 0 ) {
		    $firstsum[$i] = 0;
		    $secondsum[$i] = 0;
		}
		else {
		    $firstsum[$i] = $firstsum[ $i - 1 ] + log( ( $n - $i ) + 1 );
		    $secondsum[$i] = $secondsum[ $i - 1] + log($i);
		}
	    }
	    my $result = 0;
	    for ( my $i = $n ; $i >= $k ; $i-- ) {
		my $addme = $i * log($p);
		#$GU->user_info(3,"1. $addme\n");
		$addme += ( $n - $i ) * log( 1 - $p );
		#$GU->user_info(3,"2. $addme\n");
		$addme -= $secondsum[$i];
		#$GU->user_info(3,"3. $addme\n");
		$addme += $firstsum[$i];
		#$GU->user_info(3,"4. $addme\n");
		$addme = exp($addme);
		#$GU->user_info(3,"5. $addme\n");
		$result += $addme;
		#$GU->user_info(3,"R. $result\n\n");
		
	    }
	    #$GU->user_info(3,"RESULT = ".$result."\n");
	    return $result;
	} # get_binomial_probability #
	
	
	method one_sided_ks_test(ArrayRef $data, ArrayRef $range, Str $type){
		
		#Added by PEB, June 2012. Tests whether a list of numbers is distributed between two specified values
		#according to the specified theoretical distribution, using one sided Kolmogorov-Smirnov test
		#Used to test if the distribution of motif start positions in a promotor region is random, ie linear
		
		#data - measured distribution to be tested
		#range - ascending list of possible values of theoretical distribution
		#type - expected shape of  theoretical distribution. 'linear' means data randomly distributed over range
		my @range = @$range;
		my $lower = $range[0];
		my $upper = $range[$#range];
		my $distribution_type = lc($type);
		
		if($lower >= $upper){
			die("ERROR: statistics::one_sided_ks_test: Invalid cdf range specified\n");	
		}
		
		my @sorted_data = sort {$a <=> $b} @$data;
		my $num_values = scalar (@sorted_data);
			
		if($sorted_data[0] < $lower){
			print "ERROR: statistics::one_sided_ks_test: ".$sorted_data[0]." less than ".$lower;
		}
		if($sorted_data[$num_values-1] > $upper){
			print "ERROR: statistics::one_sided_ks_test: ".$sorted_data[$num_values-1]." greater than ". $upper;
		}
			

		print 'Distribution over range '.$lower.' to '.$upper.' expected to be '.$distribution_type.'\n';

		#Compare distributions to find largest difference
		my $D_value = 0;
		my $idx = 0;
		
		for(my $i = 0; $i <= $#range; $i++){
			my $cumulative_expected;
			if($distribution_type eq "linear"){
				#The expected cumulative distributon if linear
				$cumulative_expected = ($i+1)/($#range+1);
			}else{
				die("ERROR: statistics::one_sided_ks_test: Invalid distribution type ".$distribution_type." specified\n");		
			}
			
			#This is the proportion of values in data expected to be <= $range[$i]
			#find actual proportion
			my $val = $range[$i];
			while($idx < $num_values){
				if($sorted_data[$idx] <= $val){
					$idx++;
				}else{
					last;
				}
			}
			my $cumulative_observed = $idx/$num_values;
			
			if(abs($cumulative_observed-$cumulative_expected) > $D_value){
				#record largest difference between the two
				$D_value = 	abs($cumulative_observed-$cumulative_expected);
			}	
		}
		print 'Test statistic: '.$D_value.'\n';
		
		#Now calculate the corresponding probability of getting this deviation from linear
		#according to algorithm from Marsaglia et al, Journal of Statisitcal Software, Vol 8:18
		
		my $probability_value;
		my $s = $D_value**2 * $num_values;

		if(($s>7.24) || (($s>3.76 && ($num_values>99)))){
			#shortcut for big D and high values of n, ie highly siginificant p value
			$probability_value = 2 * exp(-(2.000071+0.331/sqrt($num_values)+1.409/$num_values) * $s)
		}else{
			my $k = int($num_values * $D_value) + 1;
			my $m = int(2 * $k - 1);
			my $h = $k - $num_values * $D_value;
		
			my @H = ();
			my @Q = ();
			my $eQ;
	
			for(my $i = 0; $i < $m; $i++){
				for(my $j = 0; $j < $m; $j++){
					if($i-$j+1 < 0){
						$H[$i*$m+$j] = 0;
					}else{
						$H[$i*$m+$j] = 1;
					}
				}
			}
	
			for(my $i = 0; $i < $m; $i++){
				$H[$i*$m] -= $h**($i+1);
				$H[($m-1)*$m+$i] -= $h**($m-$i); 	
			}
			$H[($m-1)*$m] += ((2*$h-1 > 0) ? (2*$h-1)**$m : 0);
	
			for(my $i = 0; $i < $m; $i++){
				for(my $j = 0; $j < $m; $j++){
					if($i-$j+1 > 0){
						for(my $g=1; $g<=($i-$j+1); $g++ ){
							$H[$i*$m+$j] /= $g;
						}
					}
				}
			}
	
			my $eH = 0;
			$self->mPower(\@H, $eH, \@Q, \$eQ, $m, $num_values);
			$s = $Q[($k-1)*$m+$k-1];
			for(my $i=1; $i<=$num_values;$i++){
				$s = $s*$i/$num_values;
				if($s<1e-140){
					$s*=1e140;
					$eQ-=140;
				}
			}
	
			$s*=10**$eQ;
			$probability_value = 1-$s;
			
			print 'P value: '.$probability_value.'\n';
			
			return $probability_value;
		}
			
	} #one_sided_ks_test


	method mMultiply(ArrayRef $A, ArrayRef $B, ArrayRef $Result, Int $m) {
		
		#Helper function for K-S test
	
		for(my $i=0;$i<$m;$i++){
			for(my $j=0;$j<$m;$j++){
				my $s = 0.0;
				for(my $k=0;$k<$m;$k++){
					$s+=$$A[$i*$m+$k]*$$B[$k*$m+$j];
				}
				$$Result[$i*$m+$j] = $s;	
			}
		}	
	} #mMultiply

	method mPower(ArrayRef $A, Int $eA, ArrayRef $Result, ScalarRef $eV, Int $m, Int $n) {
	
		#Helper function for K-S test
	
		my $eB;

		if($n==1){
			for(my $i=0;$i<($m*$m);$i++){
				$$Result[$i] = $$A[$i];	
			}	
			$$eV = int($eA);
			return;
		}
	
		#calls itself recursively
		$self->mPower($A, $eA, $Result, $eV, $m, int($n/2));
	
		my @B = ();
		$self->mMultiply($Result, $Result, \@B, $m);
		$eB = int(2 * $$eV);
	
		if($n%2 == 0){
			for(my $i=0; $i<($m*$m);$i++){
				$$Result[$i] = $B[$i];	
			}
			$$eV = int($eB);
		}else{
			$self->mMultiply($A, \@B, $Result, $m);
			$$eV = int($eB+$eA);
		}
		if($$Result[($m/2)*$m + ($m/2)] > 1e140){
			for(my $i=0;$i<($m*$m);$i++){
				$$Result[$i] = $$Result[$i] * 1e-140;	
			}
			$$eV += 140;	
		}			
	} # mPower
	
} # Statistics #

class Generic_Weight_Matrix {
	
	use constant {FALSE => 0,
		  TRUE	=> 1};	
	
	has 'wm_identifier' => (is => 'rw', isa => 'Str', required => 1); # Trigger (trigger => \&private_clear_matrices) removed, as it prevents the adding of counts to a custom matrix. 
	has 'wm_freq_dist' => (is => 'rw', isa => 'ArrayRef[ArrayRef]', clearer => 'private_clear_wm_freq_dist');
	
	
	method generate_logo(Str $output_dir){

	    use GD;
	    
	    if (!defined $self->wm_freq_dist) {
			die 'cannot produce a sequence logo if the frequency matrix is not defined.';
	    }
	
	    my @icm_matrix = $self->private_compute_icm();
	    my $pssm_id = $self->wm_identifier();
		
		my @transposed_icm_matrix;
		
		for (my $i=0; $i < 4; $i++){
			
			for (my $j=0; $j < scalar(@icm_matrix); $j++){
				
				$transposed_icm_matrix[$i][$j] = $icm_matrix[$j][$i];
			}
		}
		
		my $logo_filepath = $output_dir."/$pssm_id.gif";  #changed from .png as couldn't get libpng to work with perl on Mac
		
		draw_logo(\@transposed_icm_matrix, -file=> $logo_filepath, 
						-full_scale =>2.25,
						-xsize=>500,
						-ysize =>250,  
						-x_title=>"position", 
						-y_title=>"bits");
		
		return $logo_filepath;
		
		sub draw_logo {
			
			#no strict;
			my $matrixref = shift;
			my %args = (-xsize      => 600,
			-full_scale => 2.25,
			-graph_title=> "",
			-x_title    => "",
			-y_title    => "",
			@_);
			# Other parameters that can be specified:
			#       -ysize -line_width -margin
			# do not have a fixed default value 
			#   - they are calculated from xsize if not specified
			
			# draw postscript logo if asked for
			#if ($args{'-ps'} || $args{'-pdf'}){
			# return _draw_ps_logo($self, @_);  
			#}
			
			#require GD;
			
			
			my ($xsize,$FULL_SCALE, $x_title, $y_title)  = @args{qw(-xsize -full_scale -x_title y_title)} ;
			
			my $PER_PIXEL_LINE = 300;
			
			# calculate other parameters if not specified
			
			my $line_width = ($args{-line_width} or int ($xsize/$PER_PIXEL_LINE) or 1);
			my $ysize      = ($args{-ysize} or $xsize/1.6); 
			# remark (the line above): 1.6 is a standard screen x:y ratio
			my $margin     = ($args{-margin} or $ysize*0.15);
			
			my $image = GD::Image->new($xsize, $ysize);
			my $white = $image->colorAllocate(255,255,255);
			my $black = $image->colorAllocate(0,0,0);
			#my $motif_size = $self->pdl_matrix->getdim(0);
			my $motif_size = scalar(@{$$matrixref[0]});
			my $font = ((&GD::gdTinyFont(), &GD::gdSmallFont(), &GD::gdMediumBoldFont(), 
			&GD::gdLargeFont(), &GD::gdGiantFont())[int(($ysize-50)/100)]
			or &GD::gdGiantFont());
			my $title_font = ((&GD::gdSmallFont(), &GD::gdMediumBoldFont(), 
			&GD::gdLargeFont(), &GD::gdGiantFont())[int(($ysize-50)/100)]
			or &GD::gdGiantFont());
			
			
			# WRITE LABELS AND TITLE
			
			# graph title   #&GD::Font::MediumBold
			$image->string($title_font,
			$xsize/2-length($args{-graph_title})* $title_font->width() /2,
			$margin/2 - $title_font->height()/2,
			$args{-graph_title}, $black);
			
			# x_title
			$image->string($font,
			$xsize/2-length($args{-x_title})*$font->width()/2,
			$ysize-( $margin - $font->height()*0  - 5*$line_width)/2 
			- $font->height()/2*0,
			$args{-x_title}, 
			$black);
			# y_title
			$image->stringUp($font,
			($margin -$font->width()- 5*$line_width)/2 
			- $font->height()/2 ,
			$ysize/2+length($args{'-y_title'})*$font->width()/2,
			$args{'-y_title'}, $black);
			
			
			# DRAW AXES
			
			# vertical: (top left to bottom right)
			$image->filledRectangle($margin-$line_width, $margin-$line_width, 
			$margin-1, $ysize-$margin+$line_width, 
			$black);
			# horizontal: (ditto)
			$image->filledRectangle($margin-$line_width, $ysize-$margin+1, 
			$xsize-$margin+$line_width,$ysize-$margin+$line_width,
			$black);
			
			# DRAW VERTICAL TICKS AND LABELS
			
			# vertical axis (IC 1 and 2) 
			my $ic_1 = ($ysize - 2* $margin) / $FULL_SCALE;
			foreach my $i (1..$FULL_SCALE)  {
				$image->filledRectangle($margin-3*$line_width, 
				$ysize-$margin - $i*$ic_1, 
				$margin-1, 
				$ysize-$margin+$line_width - $i*$ic_1, 
				$black);
				$image->string($font, 
				$margin-5*$line_width - $font->width,
				$ysize - $margin - $i*$ic_1 - $font->height()/2,
				$i,
				$black);
			}
			
			# DRAW HORIZONTAL TICKS AND LABELS, AND THE LOGO ITSELF 
			
			# define function refs as hash elements
			my %draw_letter = ( A => \&draw_A,
			C => \&draw_C,
			G => \&draw_G,
			T => \&draw_T );
			
			my $horiz_step = ($xsize -2*$margin) / $motif_size;		

			foreach my $x (0..$motif_size)  {
				
			    if ($x>70) {
				
				print "cut motif length down to avoid crash (".$args{-file}.")\n";
				   # Kashi, Sascha: trying to print a motif-logo for R02192 (TRANSFAC) caused segmentation
				   #                faults. The error occurred in the draw_G-function when it tries to fill
				   #                the arc of the G at one of the last positions in this consensus sequence.
                                   #                Resolution was not obvious, therefore restricted this code. No need for
                                   #                logos this long at this point in time.
				last;
			    }

				$image->filledRectangle($margin + $x*$horiz_step, 
				$ysize-$margin+1, 
				$margin + $x*$horiz_step+ $line_width, 
				$ysize-$margin+3*$line_width, 
				$black);
				last if $x==$motif_size;								
				
				# get the $i-th column of matrix
				my %ic; 
				#($ic{A}, $ic{C}, $ic{G}, $ic{T}) = list $self->pdl_matrix->slice($x);
				$ic{A} = $$matrixref[0][$x];
				$ic{C} = $$matrixref[1][$x];
				$ic{G} = $$matrixref[2][$x];
				$ic{T} = $$matrixref[3][$x];
				
				# sort nucleotides by increasing information content
				my @draw_order = sort {$ic{$a}<=>$ic{$b}} qw(A C G T);
				#$GU->user_info(3,Dumper (%ic));
				# draw logo column
				my $xlettersize = $horiz_step /1.1;
				my $ybottom = $ysize - $margin;		
			
				foreach my $base (@draw_order)  {
					#$GU->user_info(3,$base."\n");
					my $ylettersize = int($ic{$base}*$ic_1 +0.5);
					#$GU->user_info(3,$ylettersize."\n");
					next if $ylettersize ==0;
					
	    			# draw letter				
					$draw_letter{$base}->($image,
							      $margin + $x*$horiz_step,
							      $ybottom - $ylettersize,
							      $xlettersize, $ylettersize, $white);
					$ybottom = $ybottom - $ylettersize-1;
						    			
				}			
				
				if ($args{'-error_bars'} and ref($args{'-error_bars'}) eq "ARRAY")  {
					my $sd_pix   = int($args{'-error_bars'}->[$x]*$ic_1);
					my $yt     = $ybottom - $sd_pix+1;
					my $yb  = $ybottom + $sd_pix-1;
					my $xpos     = $margin + ($x+0.45)*$horiz_step;
					my $half_width;
					
					if ($yb > $ysize-$margin+$line_width)  {
						$yb = $ysize-$margin+$line_width
					}
					else {
						$image->line($xpos - $xlettersize/8, $yb, 
						$xpos + $xlettersize/8, $yb, 
						$black);
					}
					
					$image->line($xpos, $yt, $xpos, $yb, $black);
					$image->line($xpos - 1 , $ybottom, $xpos+1, $ybottom, $black);
					$image->line($xpos - $xlettersize/8, $yt, 
					$xpos + $xlettersize/8, $yt, 
					$black);
					
					
				}
				
				# print position number on x axis
				$image->string($font,
				$margin + ($x+0.5)*$horiz_step - $font->width()/2,
				$ysize - $margin +5*$line_width,
				$x+1,
				$black);
			}
			
			# print $args{-file};
			if  ($args{-file}) {  
				open (GIFFILE, ">".$args{-file})
				or die;
				print GIFFILE $image->gif; #changed from ->png as couldn't get libpng to work with perl on Mac
				close GIFFILE;
			}
			return $image;
		}
		
		
		sub draw_C  {
			my ($im, $x, $y, $xsize, $ysize, $white) = @_;
			my $blue = $im->colorAllocate(0,0,204);
			$im->arc($x+$xsize*0.54, $y+$ysize/2,1.08*$xsize,$ysize,0,360,$blue);
			$im->fill($x+$xsize/2, $y+$ysize/2, $blue);
			if ($ysize>12) {
				$im->arc($x+$xsize*0.53, $y+$ysize/2, 
				0.55*$xsize, (0.625-0.625/$ysize)*$ysize,
				0,360,$white);
				$im->fill($x+$xsize/2, $y+$ysize/2, $white);
				$im->filledRectangle($x+$xsize/2, $y+$ysize/2.8+1, 
				$x+$xsize*1.1, $y+(2.8*$ysize/4)-1,
				$white);
			}
			elsif ($ysize>3)  {
				$im->arc($x+$xsize*0.53, $y+$ysize/2, 
				(0.75-0.75/$ysize)*$xsize, (0.725-0.725/$ysize)*$ysize,
				0,360,$white);
				$im->fill($x+$xsize/2, $y+$ysize/2, $white);
				$im->filledRectangle($x+$xsize*0.25, $y+$ysize/2, 
				$x+$xsize*1.1, $y+$ysize/2,
				$white);
				
			}
			return 1;
		}
		
		sub draw_G  {
			my ($im, $x, $y, $xsize, $ysize, $white) = @_;
			my $yellow = $im->colorAllocate(255,179,0);
			$im->arc($x+$xsize*0.54, $y+$ysize/2,1.08*$xsize,$ysize,0,360,$yellow);
			$im->fill($x+$xsize/2, $y+$ysize/2, $yellow);
			if ($ysize>20) {
				$im->arc($x+$xsize*0.53, $y+$ysize/2, 
				0.55*$xsize, (0.625-0.625/$ysize)*$ysize,
				0,360,$white);
				$im->fill($x+$xsize/2, $y+$ysize/2, $white);
				$im->filledRectangle($x+$xsize/2, $y+$ysize/2.8+1, 
				$x+$xsize*1.1, $y+$ysize/2-1,
				$white);
			}
			elsif($ysize>3)  {
				$im->arc($x+$xsize*0.53, $y+$ysize/2, 
				(0.75-0.75/$ysize)*$xsize, (0.725-0.725/$ysize)*$ysize,
				0,360,$white);
				$im->fill($x+$xsize/2, $y+$ysize/2, $white);
				$im->filledRectangle($x+$xsize*0.25, $y+$ysize/2, 
				$x+$xsize*1.1, $y+$ysize/2,
				$white);
				
			}
			$im->filledRectangle($x+0.85*$xsize, $y+$ysize/2,
			$x+$xsize,$y+(3*$ysize/4)-1,
			$yellow);
			$im->filledRectangle($x+0.6*$xsize, $y+$ysize/2,
			$x+$xsize,$y+(5*$ysize/8)-1,
			$yellow);
			return 1;
		}
		
		sub draw_A {
			
			my ($im, $x, $y, $xsize, $ysize, $white) = @_;
			my $green = $im->colorAllocate(0,204,0);
			my $outPoly = GD::Polygon->new();
			$outPoly->addPt($x, $y+$ysize);
			$outPoly->addPt($x+$xsize*.32, $y);
			$outPoly->addPt($x+$xsize*.68, $y);
			$outPoly->addPt($x+$xsize, $y+$ysize);
			$outPoly->addPt($x+0.75*$xsize, $y+$ysize);
			$outPoly->addPt($x+0.635*$xsize, $y+0.75*$ysize);
			$outPoly->addPt($x+0.355*$xsize, $y+0.75*$ysize);
			$outPoly->addPt($x+0.25*$xsize, $y+$ysize);
			$im->filledPolygon($outPoly, $green);
			if ($ysize>8)  {
				my $inPoly = GD::Polygon->new();
				$inPoly->addPt($x+$xsize*.5, $y+0.2*$ysize);
				$inPoly->addPt($x+$xsize*.40, $y+0.6*$ysize-1);
				$inPoly->addPt($x+$xsize*.60, $y+0.6*$ysize-1);
				$im->filledPolygon($inPoly, $white);
			}
			return 1;
		}
		
		sub draw_T {
			
			my ($im, $x, $y, $xsize, $ysize, $white) = @_;
			my $red = $im->colorAllocate(204,0,0);
			$im->filledRectangle($x, $y, $x+$xsize, $y+0.16*$ysize, $red);
			$im->filledRectangle($x+0.35*$xsize, $y, $x+0.65*$xsize, $y+$ysize, $red);
			return 1;
		}
		
	} # generate_logo #
	
	
	method get_information_content() {
	    # returns information content for weight matrices that have a frequency matrix
	    
	    if (!defined $self->wm_freq_dist) {
		die 'require frequency matrix to compute information content';
	    }
	    my $matrix = $self->wm_freq_dist;
	    my $result = 0;
	    foreach my $position (@{$matrix}) {
		my $tic = $self->private_compute_total_information_content_for_position($position);
		$result = $result + $tic;
	    }
	    return $result;
	} # get_information_content #

	method private_compute_icm(){
		
		my $matrix = $self->wm_freq_dist;
		my @icm;
		
		foreach my $position (@{$matrix}){
			
			my $tic = $self->private_compute_total_information_content_for_position($position);
			my @ic4p;
			
			foreach my $base (@{$position}){
				
				
				my $ic = $base * $tic;
				#$GU->user_info(3,$base." * ".$tic." = ". $ic." ic= ".$ic."\n");
				#$GU->user_info(3,$ic."\t");
				push(@ic4p, $ic);
			}
			push(@icm, \@ic4p);
			#$GU->user_info(3,"\n");
			
		}
		
		return @icm;
		
	} # compute_icm #
	
	method private_compute_total_information_content_for_position(ArrayRef $position) {
	    
	    my $addme = 0;
	    
	    foreach my $base (@{$position}){
		if ($base != 0) {
		    $addme += ( $base * ( log($base) / log(2) )  );
		}
	    }
	    
	    $addme = $addme + 2;
	    
	    return $addme;
	    
	} # private_compute_total_information_content_for_position #
	
	
	method get_pssm_length () {
	    my $wm_length;
	    
	    if (!defined $self->wm_freq_dist) {		
		    die "cannot establish weight matrix length.";
	    } else {
			$wm_length = @{$self->wm_freq_dist};
	    }
	    return $wm_length;
	} # get_pssm_length #


	sub private_clear_matrices() {
	    my ($self, $arg) = @_; 

	    $self->private_clear_wm_freq_dist(); # clearer methods are provided by Moose if a name is provided (as done above)
	} # private_clear_matrices #

} # Generic_Weight_Matrix #


