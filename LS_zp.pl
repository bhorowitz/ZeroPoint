#!/usr/bin/perl -w

#
# LaSilla-Quest Calibration Program 
# Last Revised: July 12, 2011, Ben Horowitz (benjamin.a.horowitz@yale.edu, horowitz.ben@gmail.com)
#

# This program is designed to take a list (RRablist.dat) with columns: [ID] [RA] [Dec] [UT] [JD] [MAG] and output a file which lists the zero point
# and error for each ID. Although it worked well as of July 12, 2011, changes in dependancies (QUEST pipeline output and exptime.list) should be examined 
# if problems occur when trying to run.


# Note, as of Nov 28th, 2013, the SDSS webserver doesn't process requests from DR 7. An updated version of this code has been folded into the pipeline
# programs managed by David Rabinowitz.

# loading packages

use strict;
use warnings;
use FileHandle;
use File::Basename;
use Math::Trig;
use POSIX qw(log10);

#global variables

my $idlist;
if ($#ARGV == -1) {
    $idlist = new FileHandle "< calib190_195.dat" or die "Couldn't open survey_10good.txt:$!\n";
} elsif ($#ARGV == 0) {
    $idlist = new FileHandle "< $ARGV[0]" or die "Couldn't open $ARGV[0]:$!\n";
} else {
    print STDERR "Usage: LS_zp.pl <targetlist>\n";
    exit(1);
}

#__________________________________________________________
# subroutines
#
sub get_coords {
  #looks at the ID list to find RA and DEC to be used to find the field.
	my ($id) = @_;
	my $candlist = new FileHandle "< calib190_195.dat"  or die "Couldn't open RRLcands:$!\n";
	foreach (<$candlist>) {
		my @line = split ' ', $_;
		my $id_cand = $line[0];
		my $ra = $line[1];
		my $dec = $line[2];
#		print "$ra\t$dec\n";
		if ($id eq $id_cand) {
			return($ra, $dec)
		}
	}
}

	
sub get_catalog {
  #Uses the JD and UT to find the entire QUEST catalog for that time/date.
	my ($ut, $jd) = @_;
	print "$ut\n";
	my @split_ut = split (/(.{8})/, $ut, 2);
	my $date = $split_ut[1];
	my $image = $split_ut[2];
#	my @split_jd = split (/\./, $jd);
	my @split_jd = split (/(.{5})/, $jd, 2);
	my $decimal = $split_jd[2];
	print "This is the decimal $decimal\t and jd $jd\t and $ut\n";
	print "$date\t$image\n";
	if ($decimal < .99) {
		
		my $catalog = "/group/astro/pipeout/2010/$date/images$image*/*cat.gz";
		return ($catalog);
		}
	if ($decimal > .990) {
		my $date1 = $date +1;
		my $catalog = "/group/astro/pipeout/2010/$date1/images$image*/*cat.gz";
		return ($catalog);
		}
}

sub get_quest_fs_zps {
        # This section of code is to select the potential quest field stars and compare them against SDSS.
	my ($id, $sloan, $ut, $ccd) = @_;
	print "catalog\n";
	my $zpsfile = new FileHandle "> $id/$id.zps" or die "Couldn't create $id.quest:$!\n";
	
#	my $testfile = new FileHandle "> $id/$id.tst" or die; 
#	printf ( $testfile " ra/dec, mag_im, corr_mag, mag_sloan, zp \n"); 
	my $im = new FileHandle "< tmpcat"  or die "Couldn't open tmpcat:$!\n";
	print "1\n";
	while (<$im>) {
		
		#print "CCD: $ccd\n";
		my @line_im = split ' ', $_;
		my $ccd_im = $line_im[0];
		my $ra_im = $line_im[1];
		my $dec_im = $line_im[2];
		my $mag_im = $line_im[10];
		my $mag_error = $line_im[20];
		
		foreach (<$sloan>){
		        # Here Sloan is being looked at to find the matching sloan stars to the quest stars. 
			my @line_sloan = split ' ', $_;
			my $ra_sloan = $line_sloan[0];
			my $dec_sloan = $line_sloan[1];
			my $mag_sloan = $line_sloan[2];
			my $sloan_error = $line_sloan[3];
			my $g_r = $line_sloan[6];
			
			my $dra = ($ra_sloan - $ra_im)*3600;
			my $ddec = ($dec_sloan - $dec_im)*3600;
			#print $logfile2 "Coordinate: $dra \t $ddec \t  CCD: $ccd \t $ccd_im \n"; 
			if (abs($dra)<2  and abs($ddec)<2 and $ccd_im eq $ccd) {
			  #ccd match is important, as if the field is overlapping onto another CCD, a different ZP might be found.
				print "4\n";
				my @foo = split ' ', `grep $ut exptime.list`;
        			my $exptime = $foo[1];
;
				my $expcor = 2.5 * log($exptime/60.0)/log(10.0);
       				my $corr_mag = $mag_im + $expcor;
				my $zp = ($mag_sloan - $corr_mag);
#				print "$exptime\t$mag_im\t$corr_mag\n";
#				print "This is the QUEST error foreach Field Star: $mag_error\n";
				print $zpsfile "$zp\t$sloan_error\t$mag_sloan\t$corr_mag\t$mag_error\t$ra_sloan\t$dec_sloan\t$dra\t$ddec\t$g_r\n";
				#printf( $testfile "%10.3f %10.3f %6.3f %5.3f %6.3f %5.3f\n", $ra_im, $dec_im, $mag_im,$corr_mag, $mag_sloan, $zp);
				


#				print "$zp\t$g_r\t@line_im\n";
			}
		seek($sloan, 0, 0);	
		}

	}
}	

sub clip_zp {
        # This section of code is designed to take out any outliers from the zero point list.
	my ($id) = @_;
	my $clippedfile = new FileHandle "> $id/$id.clipped" or die "Couldn't open clippedfile:$!\n";
        my $inputfile = new FileHandle "< $id/$id.zps" or die "Couldn't open inputfile:$!\n";
	my @column= ();
	my @columna = ();

	foreach (<$inputfile>) {
      		my @line = split ' ', $_;
      		my $zp = $line[0];
		if ($zp != -1) {
			push @column, $zp;	
			}

	}

       		my $total1 = 0;
		foreach my $zp (@column) {
		$total1 += $zp;
		}
		my $mean1 = $total1 / (scalar @column);

		# find the mean of the squares of the differences between each number and mean1
		my $total2 = 0;
		foreach my $zp (@column) {
		$total2 += ($mean1-$zp)**2;
		}
		my $mean2 = $total2 / (scalar @column);
		my $std_dev = sqrt($mean2);

#		print "$mean1\t$std_dev\n";

		foreach my $zp(@column) {
			my $diff = $mean1 - abs($zp);
			my $max_diff = 1*$std_dev;
			print "$max_diff\t$diff\t$mean1\t$std_dev\n";
			if (abs($diff) < $max_diff) {
				push @columna, $zp;
				}
		}		
				foreach my $zp (@columna) {
				print $clippedfile "$zp\n";
				}
}

sub weighted_zp {
        #This section of code calculates a weighted zp and error.
	my ($id) = @_;
	my @column = ();   
	my @row = ();
	my @part = ();
	my $clippedoutput = new FileHandle "< $id/$id.clipped" or die "Couldn't open clippedfile:$!\n";
	my $vcal = new FileHandle "> $id/$id.vcal" or die "Couldn't open zpsfile:$!\n"; 
	while (<$clippedoutput>){
		my @line = split ' ', $_;
		my $zp = $line[0];
		my $zpsfile = new FileHandle "< $id/$id.zps" or die "Couldn't open zpsfile:$!\n";
		foreach (<$zpsfile>) {
			my @lcline = split ' ', $_;
			my $lczp = $lcline[0];
			my $lcerr = $lcline[1];
			my $lcerr2 = ($lcerr)**2;
			my $lcerr_inverse = (1/$lcerr2);
			my $zp_over_err2 = ($zp/$lcerr2);
			my $g_r = $lcline[9];
#			print "$lczp\t$zp\n";
			if ($lczp eq $zp) {
				print $vcal "$lczp\t$g_r\n";
				push @column, $lcerr_inverse;
				push @row, $zp_over_err2;
				push @part, $lcerr_inverse;
			}
			seek($zpsfile, 0, 0); 
		}
	}
#	print "@column\t@row\n";
	my $total1 = 0;
	foreach my $err (@column) {
	$total1 += $err;
	}
	my $mean1 = $total1 / (scalar @column);

	# find the mean of the squares of the differences between each number and mean1
	my $total2 = 0;
	foreach my $zp (@row) {
	$total2 += $zp;
	}
	my $weighted_zp = $total2 / $total1;
	
	my $total3 = 0;
	foreach my $err2d (@part) {
	$total3 += $err2d;
	}
	my $zp_err = 1/sqrt($total3);
#	print "This is the weighted zeropoint error for the field stars: $zp_err\n";
	return($weighted_zp, $zp_err);
	
#	print "$total2\t$total1\n";
		
		
#	my $total1 = 0; 
#	foreach my $lcerr2(@column) {
#	my $total1 += ($lcerr2);
#	}
#	
#	my $total2 = 0;
#	foreach my $zp_over_err2(@row) {
#	my $total2 += ($zp_over_err2);
#	}
#	my $weighted_zp = ($total2/$total1);
#	return($weighted_zp);
}			



sub get_ut {
        # this section extracts the Universal Time from the id list.
	my ($id) = @_;

	my $lctarget = "calib190_195.dat";
	my $lcfile = new FileHandle "< $lctarget" or die "Couldn't open $lctarget:$!\n";
	foreach (<$lcfile>) {
		if(/^#/){
			next;
			print "chomp";
		}
  		
		my @result_line = split ' ', $_;
		my $id1 = $result_line[0];
		my $jd = $result_line[3];
		my $ut = $result_line[4];
		my $lcmag = $result_line[5];
#		print "$maxmag\t$result_mag\n";
		if ($id eq $id1) {
#			print "$ut\t$jd\n";
			return ($ut, $jd, $lcmag);
			$lcfile -> close;
		}
	}
}

sub get_maxmag {
        # this section of code finds the matching magnitude for the id from the id list. "maxmag" is a bit of misnomer. 
	my ($id) = @_;
#	print "$id\n";
	my @column= ();
	my $lctarget = "calib190_195.dat";
	my $lcfile = new FileHandle "< $lctarget" or die "Couldn't open lcfile:$!\n";
        while (<$lcfile>) {
		if(/^#/){
			next;
			print "chomp";
		}
		chomp;
                my @line = split ' ', $_;
                my $mag = $line[2];
                if ($mag != 'nan') {
                        push @column, $mag;
                }
	}
	#my @sorted = sort @column; 
	my $maxmag = $column[0];
	print "$maxmag\n";
	return($maxmag, $lcfile);
}	

sub mag_SDSS_selection {
        #this section further cuts the possible sloan field stars.
	my ($id) = @_;
	my $sdss_output = new FileHandle "<$id/$id.sdss" or die "couldn't open SDSSoutput:$!\n";
	my $sloanfile = new FileHandle ">$id/$id.sloan" or die "couldn't create sloanfile:$!\n";
	foreach (<$sdss_output>) {
		my @line = split ' ', $_;
		my $ra = $line[0];
		my $dec = $line[1];
		my $g_mag = $line[4];
		my $g_err = $line[5];
		my $r_mag = $line[6];
		my $r_err = $line[7];
		my $v_err = sqrt($r_err**2 + $g_err**2);
		my $g_r = $g_mag - $r_mag;
		my $calc_v = ($g_mag - .57844*($g_mag - $r_mag) - .0038);
		#color and magnitude cuts are used to ensure good quality field stars.
		if ($g_mag > 10.0 and $g_mag < 24.0 and $g_r > .1 and $g_r < 1.1) {
#			print "Foreach fieldstar: g error: $g_err\t r error: $r_err\t V error: $v_err\n"; 
			print $sloanfile "$ra\t$dec\t$calc_v\t$v_err\t$g_mag\t$g_err\t$g_r\n";
			}
	}
}


sub SDSS_query {
        # this section queries the sloan database through the internet, and downloads acceptable stars into a buffer (out.csv).
	my ($ra1, $dec1, $id) = @_;
	my $ra_low = ($ra1 - 0.083333);
	my $ra_high = ($ra1 + 0.083333);
	my $dec_low = ($dec1 - 0.083333);
	my $dec_high = ($dec1 + 0.083333);
	my $outfilename = "$id/$id.sdss";
	
	#these are for converting to log based magnitudes
	my $bu = 1.4E-10;
	my $bg = 0.9E-10;
	my $br = 1.2E-10;
	my $bi = 1.8E-10;
	my $bz = 7.4E-10;

	#first, download sloan data in the appropriate area.
	#download all stars with mags brighter than that where the flux = 0.
	#(fainter than that we can't convert it to log based fluxes.)

#	print "Sending SDSS query...\n";

	system "wget -q -O \"out.csv\" \"http://cas.sdss.org/astro/en/tools/search/x_sql.asp?format=csv&cmd=SELECT%20ra%2Cdec%2Cu%2CErr_u%2Cflags_u%2Cg%2CErr_g%2Cflags_g%2Cr%2CErr_r%2Cflags_r%2Ci%2CErr_i%2Cflags_i%2Cz%2CErr_z%2Cflags_z%20FROM%20Star%20WHERE%20u%3C24.63%20AND%20g%3C25.11%20AND%20r%3C24.8%20AND%20i%3C24.36%20AND%20z%3C22.83%20AND%20ra%3E$ra_low%20AND%20ra%3C$ra_high%20AND%20dec%3E$dec_low%20AND%20dec%3C$dec_high\"";

#	print "Done fetching SDSS\n";

	#these are the flags that we want to cut on.
	#if the object fails in any filter (r,i,z) then we don't want it totally.
	my @badpowers;
	push @badpowers, 1; #BRIGHT
	push @badpowers, 2; #EDGE
	#push @badpowers, 3; #BLENDED
	push @badpowers, 5; #PEAKCENTER
	push @badpowers, 7; #NOPROFILE
	push @badpowers, 18; #SATURATED
	push @badpowers, 19; #NOTCHECKED
	push @badpowers, 40; #BAD_COUNTS_ERROR
	push @badpowers, 47; #PSF_FLUX_INTERP

	#this is the flag we want to make sure is set
	my $goodpower = 28; #BINNED1

	#these are for keeping track of whether an object is passing the cuts
	my $acceptable = 0;
	my $skip = 0;
	my $badobj = 0;

	my $numtot = 0;
	my $numkept = 0;

	open( INFILE, "out.csv" ) or die "Can't open out.csv: $!\n";
	open( OUTFILE, ">$outfilename" ) or die "Can't open $outfilename: $!\n";
	
	while( my $line = <INFILE> ){

	    if( $line =~ "ra" ){
	    }
	    else{
		++$numtot;
		$badobj = 0;
		$acceptable = 0;
		$skip = 0;
		my( $ra, $dec, $u, $eu, $uflags, $g, $eg, $gflags, $r, $er, $rflags, $i, $ei, $iflags, $z, $ez, $zflags ) = split( ",", $line );

		#for now, we only care about the g, r, i, and z flags.

		my %colorflags;
	#	$colorflags{$u} = {
	#	    f => $uflags,
	#	    e => $eu,
	#	};
		$colorflags{$g} = {
		    f => $gflags,
		    e => $eg,
		};
		$colorflags{$r} = {
		    f => $rflags,
		    e => $er,
		};
		$colorflags{$i} = {
		    f => $iflags,
		    e => $ei,
		};
		$colorflags{$z} = {
		    f => $zflags,
		    e => $ez,
		};

		foreach my $color ( sort keys %colorflags ){
		    my $flags = $colorflags{$color}{f};
		    my $error = $colorflags{$color}{e};

		    my $power = 64;
		    if( $flags > 2 ** 64 ){
#			print "Problem, flags number is too high\n";
			exit 0;
		    }
		    else{
	#	    print "Flags: ";
			while( ($power >= 0) && (!$skip) ){
			    if( $flags >= 2 ** $power ){
				$flags -= 2 ** $power;
	#		    print "2^$power ";
				foreach my $bad( @badpowers ){
				    if( $power == $bad ){
					$skip = 1;
	#			    print "BAD ";
				    }
				}
				if( $power == $goodpower ){
				    $acceptable = 1;
				}
				#if INTERP_CENTER (44) is set, then if COSMIC_RAY (12) 
				#is also set then we want to throw out the object.
				if( $power == 44 ){
				    push @badpowers, 12;
				}
				#if DEBLEND_NOPEAK (46) is set, then if the psf error 
				# > 0.2 then we want to throw out the object
				if( $power == 46 ){
				    if( $error > 0.2 ){
					$skip = 1;
				    }
				}
	#	      	    print "$power ";
			    }
			    --$power;
			}
		    }
		    if( $skip || !$acceptable ){
			$badobj = 1;
		    }
	#	print "\n";
		}
		if( !$badobj ){
		    ++$numkept;
		    #convert from hyperbolic trig magnitudes to log-based magnitudes
		    my ($u_c, $eu_c) = convertmags($u, $eu, $bu);
		    my ($g_c, $eg_c) = convertmags($g, $eg, $bg);
		    my ($r_c, $er_c) = convertmags($r, $er, $br);
		    my ($i_c, $ei_c) = convertmags($i, $ei, $bi);
		    my ($z_c, $ez_c) = convertmags($z, $ez, $bz);
		my $flag=0;
		    #these two cuts ensure that the field doesn't include the target star
		   unless(abs($ra1-$ra)<0.000555)
			{
			unless(abs($dec1-$dec)<0.000555)
			{
		    printf( OUTFILE "%10.6f %10.6f %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f\n", 
			    $ra, $dec, 
			    $u_c, $eu_c,
			    $g_c, $eg_c,
			    $r_c, $er_c, 
			    $i_c, $ei_c, 
			    $z_c, $ez_c );
			$flag = 10;
			}
			}
		if($flag>5.0){
		        #this just saves the sloan star information for later use if needed, note that there are possible multiples, depending on how
		        #good the astrometry is.
			my $starfile = new FileHandle "> $id/$id.star" or die "Couldn't create $id.quest:$!\n";
			print $starfile "$id \t $ra \t $dec \t $u_c \t $eu_c \t $g_c \t $eg_c \t $r_c \t $er_c \t $i_c \t $ei_c \t $z_c \t $ez_c \n";
		}
}
	    } # else if not header line
	} # while objects

#	print "Kept $numkept objects out of $numtot\n";

	# system "rm out.csv";

	close( INFILE );
	close( OUTFILE );
}

sub convertmags{
    #imperical (?) mag conversion
    my( $mag, $err, $const ) = @_;
    my $ratio = 2*$const*sinh( -1*( log($const) + 0.92103404*$mag ) );

    return( -2.5*log10($ratio), $err*2*$const*sqrt( 1+($ratio/(2*$const))**2 ) / $ratio  );

}

sub get_ccd {
        #this section gets the CCD from the target object 
	my ($id, $catalog, $ra, $dec) = @_;
	my $logfile2 = new FileHandle "> $id/$id.log" or die "Couldn't create $id.quest:$!\n";
	unlink "tmpcat";
    	system("zcat $catalog>tmpcat");
	my $im = new FileHandle "< tmpcat"  or die "Couldn't open tmpcat:$!\n";
#	print "in get_ccd sub\n";
	while (<$im>){
		chomp;
		my @line_im = split ' ', $_;
		my $ccd = $line_im[0];
		my $ra_im = $line_im[1];
		my $dec_im = $line_im[2];
		
		my $dra = ($ra - $ra_im)*3600;
		my $ddec = ($dec - $dec_im)*3600;
		print $logfile2 "RA: $ra \t $ra_im \t DEC: $dec \t $dec_im \t Coordinate: $dra \t $ddec \t  CCD: $ccd \n"; 
		if (abs($dra)<2  and abs($ddec)<2) {
			print "This is the ccd $ccd\n";
			return ($ccd);		
		}
	}
	print "no ccd?! /n";
}

sub check_sdss {
        # checks existance of .sdss file
        my ($id) = @_;
        my $sdssfile = new FileHandle "< $id/$id.sdss" or die "Couldn't open sdss file:$!\n";
        my @sdssfile = <$sdssfile>;
        my $lines = @sdssfile;
        if ($lines eq 0) {
            next;
         }
}

sub check_sloan {
        #checks existance of .sloan file
        my ($id) = @_;
        my $sloanfile = new FileHandle "< $id/$id.sloan" or die "Couldn't open sloan file:$!\n";
        my @sloanfile = <$sloanfile>;
        my $lines = @sloanfile;
        if ($lines eq 0) {
            next;
         }
}

sub check_zps {
        #checks existance of .zps file
        my ($id) = @_;
        my $zpsfile = new FileHandle "< $id/$id.zps" or die "Couldn't open zps file:$!\n";
        my @zpsfile = <$zpsfile>;
        my $lines = @zpsfile;
        if ($lines eq 0) {
            next;
         }
}

sub check_clipped {
        #checks existance of .clipped file
        my ($id) = @_;
        my $clippedfile = new FileHandle "< $id/$id.clipped" or die "Couldn't open clipped file:$!\n";
        my @clippedfile = <$clippedfile>;
        my $lines = @clippedfile;
        if ($lines eq 0) {
            next;
         }
}
sub starfile{
	my ($id, $weighted_zp, $zp_err) = @_;

	my $star = new FileHandle "< $id/$id.star"  or die "Couldn't open star!\n";
#	print "in get_ccd sub\n";
	while (<$star>){
		chomp;
		my @line_star = split ' ', $_;
		
		my $Gc = $line_star[3];
		my $Rc = $line_star[4];
		my $g_r = $Gc-$Rc;
		my $Pweighted_zp = $weighted_zp + (0.1213*($g_r)-0.0681);
	return($Pweighted_zp, $zp_err);
}
}

#______________________________________________
#main program

#my $candlist = new FileHandle "<RRLcands.idradecgr" or die "g:$!\n";

open ( MYFILE, '>> RRLYRAEmaster1');
open ( MYFILE2, '>> RRLYRAEmaster_magcorrection1');
while (<$idlist>) {
	my @line = split ' ', $_;
	my $id = $line[0];
#	my $ra = $line[1];
#	my $dec = $line[2];
	if (-e "calib190_195.dat") {
		mkdir ($id);
		print "$id\n";
		print "control \n";
		my ($ra, $dec) = get_coords($id);
	
		SDSS_query($ra, $dec, $id);
	 	#run the SDSS Query

	 	check_sdss($id);
	
		mag_SDSS_selection($id);
		#select those with a G mag between 14-19, print ra, dec, mag, err to .sloan file
			
		check_sloan($id);
			
		#my ($maxmag, $lcfile) = get_maxmag($id);
		#extract the maximum uncalibrated magnitude from lc. files
		#print "$maxmag\t$lcfile\n";
		my ($ut, $jd, $lcmag) = get_ut($id); 
		# find the ut for the max. calibrated magnitude from .result file
		print "\n get_catalog \n";
		my ($catalog) = get_catalog($ut, $jd);
		#get the frame catalog path from the ut
		print "\n get_ccd \n";
		my ($ccd) = get_ccd($id, $catalog, $ra, $dec);
	
		my $sloanfile = "$id/$id.sloan";
		my $sloan = new FileHandle "< $sloanfile" or die "Couldn't create $id.sloan:$!\n";
		print "\n get quest \n";
		get_quest_fs_zps($id, $sloan, $ut, $ccd); 
		# compare the list of sloan fs to Quest fs, take the zp (difference of mags), print to .zps file
	
		check_zps($id);
	
		clip_zp($id);		
		# take the clipped mean of the zps to knock out any variables or different chips
	
		check_clipped($id);
	
		my ($weighted_zp, $zp_err) = weighted_zp($id);

		my $newmag = $lcmag+$weighted_zp;
		#my ($Pweighted_zp, $Pzp_err) = starfile($id, $weighted_zp, $zp_err);
		my $newzp = $weighted_zp + (-0.01892*($newmag)+0.34622);
		my $newzp_err = sqrt(($zp_err)**2+(0.01892*$zp_err)**2);
		my $flag=1;
		# this is a long and ugly way to flag chips that are along crack, note that post 2010ish, the crack was "sealed" with,
		# as David Rabinowitz puts it, a piece of scotch tape, and supposedly it is much better now. Something to check out
		unless($ccd eq 14){
		unless($ccd eq 15){
		unless($ccd eq 42){
		unless($ccd eq 43){
		unless($ccd eq 70){
		unless($ccd eq 71){
		unless($ccd eq 98){
		unless($ccd eq 99){
		$flag=0;
		}}}}}}}}
		print MYFILE "$id \t $weighted_zp \t $zp_err \t $ccd \t $flag\n";
		#the second, magnitude corrected file has time and again been shown to be more accurate than without the magnitude correction
		print MYFILE2 "$id \t $newzp \t $newzp_err \t $ccd \t $flag\n";
	

	}
}

# right now, the script assumes exptime.list, coordinate, and the lc files are in the same directory as this program. 
# for each rrl there will return 4 files: .sdss (querried SDSS) .sloan (sloan field stars and mags within 14-19), 
#.zps (zeropoints of the sloan fs also found in quest's frame), 
# .clipped (clipped zeropoitns) and .abs (the absolute calibration for the rrl)
