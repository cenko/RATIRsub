#!/usr/bin/perl

# Packages
use LWP::Simple;
use URI::URL;

# Defaults
$TRUE=1;
$FALSE=0;
$CLOBBER=$TRUE;
$DTOR=3.14159265359/180;
%MAGKEYS=("GSC","FJVN",
          "USNOA","BR",
	  "USNOB","BRI",
          "2MASS","JHK",
          "TYCHO","BV");

## Not doing anything with DEBUG right now
#$DEBUG=$TRUE;
#$DEBUG=$FALSE;

$addnums=$FALSE;
$addmags=$FALSE;
$radius=12;
$catalog="GSC";
$color="cyan";
$regfile="";
$fullcat="";
@now=gmtime(time());
# Following formula is correct only to ~1day
$epoch=(1900+$now[5]+$now[7]/365.0);

($exec)=($0=~/\/?([^\/]+)$/);

################
# Command Line #
################

&usage if (@ARGV<3);

while(@ARGV) {
    $_=shift(@ARGV);
    if (/^[-\+]\D/) {
	# Options
	if (/^[-\+]d/i) {
	    # Database to use
	    $code=shift(@ARGV);
	    if ($code=~/^g/i) {
		$catalog="GSC";
	    }elsif ($code=~/^u.+a/i) {
		$catalog="USNOA";
	    }elsif ($code=~/^u.+b/i) {
		$catalog="USNOB";
	    }elsif ($code=~/^2/) {
		$catalog="2MASS";
	    }elsif ($code=~/^t/i) {
		$catalog="TYCHO";
	    }else{
		print STDERR "# Unknown catalog code, defaulting to GSC\n";
		$catalog="GSC";
	    }
	    print STDERR "# Retrieving the $catalog catalog\n";
	}elsif(/^[-\+]c/i) {
	    # Color for regions
	    $color=shift(@ARGV);
	    print STDERR "# Using color $color for regions\n";
	}elsif(/^[-\+]e/i) {
	    # Epoch for Proper Motions
	    $epoch=shift(@ARGV);
	    print STDERR "# Using Epoch $epoch for proper motions\n";
	}elsif(/^[-\+]r/i) {
	    # Search radius
	    $radius=shift(@ARGV);
	    print STDERR "# Setting search radius to $radius arcmin\n";
	}elsif (/^[-\+]n/i) {
	    # Number lines
	    $addnums=$TRUE;
	}elsif (/^[-\+]m/i) {
	    # Output magnitudes
	    $addmags=$TRUE;
	}elsif (/^[-\+]f/i) {
	    # Write a region file
	    $regfile=shift(@ARGV);
	    print STDERR "# Writing region file $regfile\n";
	}elsif (/^[-\+]s/i) {
	    # Save full catalog
	    $fullcat=shift(@ARGV);
	    print STDERR "# Saving full catalog to $fullcat\n";
	}else{
	    # Unrecognized options
	    print STDERR "# Unrecognized option $_... ignoring.\n";
	}
    }else{
	# Arguments
	push(@args,$_);
    }
}

&usage if (@args<3);
($ra,$dec,$out)=@args;

######################
# Parse Input Coords #
######################

# Parse RA
if ($ra=~/[:\s]/) {
    # Sexagesimal RA
    ($rah,$ram,$ras)=($ra=~/(\d+)[\s:]+(\d+)[\s:]+([\d\.]+)/);
    $rahd=ten($rah,$ram,$ras);
    $rad=$rahd*15;
}else{
    $rad=$ra;
    $rahd=$ra/15;
    ($rah,$ram,$ras)=&sixty($rahd);
}

# Parse Dec
if ($dec=~/[:\s]/) {
    # Sexagesimal Dec
    ($dcd,$dcm,$dcs)=($dec=~/([\+\-]?\d+)[\s:]+(\d+)[\s:]+([\d\.]+)/);
    $decd=ten($dcd,$dcm,$dcs);
}else{
    $decd=$dec;
    ($dcd,$dcm,$dcs)=&sixty($decd);
}


#####################
# Execute Retrieval #    
#####################

$astrom=&getastrom;

exit($astrom);


######################################################################

sub usage {
    my $diestr=(@_ ? $_[0] : "");
    $diestr .= "Usage:  $exec [-d <catalog>] [-r <radius/arcmin>] [-e <epoch>] [-s <savefile>] [-f <regfile>] [-m(ags)] [-n(umber)] [-c <color>] <ra> <dec> <outfile>\n";
    die $diestr;
}

sub sixty {
    my $dcml=$_[0];
    my $sign = ($dcml<0 ? -1 : 1);
    $dcml*=$sign;
    my $hexi=int($dcml);
    my $hexmd=60*($dcml-$hexi);
    my $hexm=int($hexmd);
    my $hexs=sprintf("%6.3f",60*($hexmd-$hexm));
    return ($sign*$hexi,$hexm,$hexs);
}

sub ten {
    my ($din,$min,$sec)=@_[0..2];
    my $sign = ($din=~/^-/ ? -1 : 1);
    $din *= $sign;
    my $deg=$sign*($din+($min+$sec/60)/60);
    return $deg;
}

sub max {
    my @a=@_;
    my $max=shift(@a);
    while (@a) {
	$el=shift(@a);
	$max=$el if ($el>$max);
    }
    return $max;
}

sub getastrom {
    # gets astrometry catalog for USNO-A2.0, 2MASS, GSC2.2
    # puts data into file $out, which contains:
    #     index ra DEC distance [mag magu]
    # where distance = output distance in arcsec
    #            mag = catalog magnitude if requested
    #           magu = catalog magnitude uncertainty

    # Size in kb of HTML doc's returned for failed requests
    my $failsize=120;
    my $failsafe;

    my $cmd, $url;
    my $RA,$DEC,$dra,$ddec,@text,$i,$r,$content;

    # For USNO forms: Supply a rectangular window width in arcsec
    if ($catalog=~/USNO/) {
	# convert to width in arcsec
	$radius*=2*60;
    }

    # USNO-A2.0 Catalog Retrieval
    if ($catalog eq "USNOA") {
	$url=url('http://asteroid.lowell.edu/cgi-bin/koehn/webnet');
	$url->query_form("rawidth"=>$radius,"decwidth"=>$radius,
			 "ra"=>"$rah $ram $ras","dec"=>"$dcd $dcm $dcs",
			 "redlo"=>"-5.0","redhi"=>"30.0","bluelo"=>"-5.0",
			 "bluehi"=>"30.0","cilo"=>"-10.0","cihi"=>"10.0",
			 "density"=>"10000000.0","shape"=>"Circular",
			 "sortby"=>"Radius","sortdirection"=>"Ascending",
			 "outstyle"=>"ASCII","focalength"=>"5.0",
			 "pixelsize"=>"15.0","catalog"=>"USNO-A2.0");

    # USNO-B1.0 Catalog Retrieval
    }elsif ($catalog eq "USNOB") {

	$url=url('http://www.nofs.navy.mil/cgi-bin/tfch3tH.cgi');
	$url->query_form("ra"=>"$rah $ram $ras","dec"=>"$dcd $dcm $dcs",
			 "equinox"=>"J2000","epoch"=>"2000.0",
			 "cextract"=>"rect","rawid"=>$radius,"decwid"=>$radius,
			 "wunits"=>"Seconds","cat"=>"USNO B1.0",
			 "surims"=>"None","getcat"=>"yes",
			 "colbits"=>"cb_id","colbits"=>"cb_ra",
			 "slf"=>"hh/dd mm ss","colbits"=>"cb_sigra",
			 "colbits"=>"cb_mep","colbits"=>"cb_mura",
			 "colbits"=>"cb_smura",
			 "colbits"=>"cb_flg","colbits"=>"cb_mag",
			 "clr"=>"R2","bri"=>"0","fai"=>"100",
			 "colbits"=>"cb_mflg",
			 "colbits"=>"cb_dstctr","skey"=>"dstctr",
			 "gzf"=>"No","cftype"=>"ASCII");

    # 2MASS Catalog Retrieval
    }elsif ($catalog eq "2MASS") {
        #$url=url('http://vizier.hia.nrc.ca/viz-bin/VizieR-4');
        $url=url('http://vizier.u-strasbg.fr/viz-bin/VizieR');
	$url->query_form("-source"=>"II/246/out","-from"=>"-3","-to"=>"4",
			 "-this"=>"-3","-oc.form"=>"dec","-order"=>"+",
			 "-meta"=>2,"-file"=>".","-out.max"=>999999,
			 "-out.form"=>"ascii table","-c"=>"$rad $decd",
			 "-c.eq"=>"J2000",
			 "-c.r"=>$radius,"-c.u"=>"arcmin","-c.geom"=>"r",
			 "-out.add"=>"_r","-sort"=>"_r",
			 "-out"=>"RAJ2000","-out"=>"DEJ2000",
			 "-out"=>"errMaj","-out"=>"errMin","-out"=>"errPA",
			 "-out"=>"Jmag","-out"=>"e_Jmag",
			 "-out"=>"Hmag","-out"=>"e_Hmag",
			 "-out"=>"Kmag","-out"=>"e_Kmag");

    # Tycho Catalog Retrieval
    }elsif ($catalog eq "TYCHO") {
	$url=url('http://vizier.u-strasbg.fr/cgi-bin/VizieR-4');
	$url->query_form("-source",=>"I/259/tyc2","-from"=>"-3","-to"=>"4",
			 "-this"=>"-3",
			 "-meta"=>"1","-file"=>".","-out.max"=>"999999",
			 "-out.form"=>"ascii table","-c"=>"$rad $decd",
			 "-c.eq"=>"J2000",
			 "-c.r"=>$radius,"-c.u"=>"arcmin",
			 "-out.add"=>"_r","-sort"=>"_r",
			 "-out"=>"RAmdeg","-out"=>"DEmdeg",
			 "-out"=>"pmRA","-out"=>"pmDE",
			 "-out"=>"e_RAmdeg","-out"=>"e_DEmdeg",
			 "-out"=>"e_pmRA","-out"=>"e_pmDE",
			 "-out"=>"BTmag","-out"=>"e_BTmag",
			 "-out"=>"VTmag","-out"=>"e_VTmag");

    # GSC Catalog Retrieval
    }elsif ($catalog eq "GSC") {
	$url=url('http://www-gsss.stsci.edu/cgi-bin/gsc22query.exe');
	$url->query_form("ra"=>$rahd,"dec"=>$decd,"r1"=>0,"r2"=>$radius,
			 "m1"=>0,"m2"=>19.5,"n"=>5000);
    }

    # Report on progress
    printout("\# Fetching $url\n");
    unless (defined ($content=get($url))) {
	killoff("Could not get $url $!\n");
    }

    # For USNO-B1.0 retrieval, need to retrieve next page:
    if ($catalog eq "USNOB") {
	@lines=split(/\n/,$content);
	@refresh=grep(/refresh/i,@lines);
	($newurl)=($refresh[0]=~/=(http.+)\"/); # "
	$newurl=~s/_fch\.html/_usnob/;
	$newcontent="";
	print "\# Checking for catalog data at $newurl\n";
	until ($newcontent=~/^FIELD/) {
	    sleep(15);
	    $newcontent=get($newurl);
	}
	$content=$newcontent;
    }

    # Save all catalog data to file if requested
    if ($fullcat) {
	open(FCAT,"> $fullcat");
	print FCAT $content;
	close(FCAT);
    }

    if ($regfile) {
	if(!open(REG,"> $regfile")){
	    print STDERR "# Failed to open region file $regfile\n";
	    $regfile="";
	}else{
	    print REG "\# Region file format: DS9 version 3.0\n";
	    print REG "global color=green font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0\n";
	}
    }

    # Reformat output
    @text=split(/\n/,$content);
    
    if (-e $out && !($CLOBBER)) {
	killoff("Clobber set to no and file $out exists\n");
    }

    if (!open(OUT,"> $out")) {
	killoff("Couldn't open $out for writing: $!",__LINE__);
    }

    my $i=0;
    foreach $_ (@text) {

	if ($catalog eq "USNOA") {
	    if (/\d{2} \d{2} [\d\.]+ +[\+\-]\d{2} \d{2} [\d\.]+/) {
		# Data line
		next if (/[me]/i); # Skip stars w/ bad magnitudes!
		# Good data
		s/^\s+//;
		@data=split(/\s+/,$_);
		$RA=15*($data[0]+$data[1]/60+$data[2]/3600);
		$decsn=($data[3]=~/^-/ ? "-" : "+");
		$DEC=abs($data[3])+$data[4]/60+$data[5]/3600;
		$DEC=$decsn.$DEC;
		($RMAG,$BMAG)=@data[6..7];
		++$i;
		$rdist{$i}=$data[8];
		$mags{$i}="$BMAG,0.1,$RMAG,0.1";
		$dataline{$i}="$RA $DEC 0.3 0.3";
	    }	

	}elsif ($catalog eq "USNOB") {
	    if (/\d{4}-\d{7} \d{2} \d{2} [\d\.]+ [+\-]\d{2} \d{2} [\d\.]+/) {
		# Data line
		# next if (/[me]/i); # Skip stars w/ bad magnitudes!
		# Good data
		s/^\s+//;
		@data=split(/\s+/,$_);
		$RA=15*($data[1]+$data[2]/60+$data[3]/3600);
		$decsn=($data[4]=~/^-/ ? "-" : "+");
		$DEC=abs($data[4])+$data[5]/60+$data[6]/3600;
		$DEC=$decsn.$DEC;
		($RAU,$DECU)=(0.001*$data[7],0.001*$data[8]);
		($GSCEP,$PMRA,$PMDEC,$PMRAU,$PMDECU)=@data[9..13];
		if (abs($PMRA)+abs($PMDEC)>0) {
		    $rapmm=($epoch-$GSCEP)/cos($DEC*$DTOR);
		    $decpmm=($epoch-$GSCEP);
		    $RA=$RA+$PMRA*$rapmm/3.6e6;
		    $DEC=$DEC+$PMDEC*$decpmm/3.6e6;
		    $RAU=sqrt($RAU**2 + ($PMRAU*$rapmm*1e-3)**2);
		    $DECU=sqrt($DECU**2 + ($PMDECU*$decpmm*1e-3)**2);
		}
		($BMAG,$RMAG,$IMAG)=@data[19,21,23];
		++$i;
		if (abs($PMRA)+abs($PMDEC)>0) {
		    # Recalculate distance to field center
		    $dra=($RA-$rad)*cos($decd*$DTOR);
		    $ddec=($DEC-$decd);
		    $rdist{$i}=sqrt(($dra**2)+$ddec**2)*3600;
		}else{
		    $rdist{$i}=$data[25];
		}
		$mags{$i}="$BMAG,0.1,$RMAG,0.1,$IMAG,0.1";
		$dataline{$i}=sprintf("%s %s %.3f %.3f",$RA,$DEC,$RAU,$DECU);
	    }

	}elsif ($catalog eq "2MASS") {
	    return -1 if (/No object found/i);
	    if (/VizieR-5/) {
		# data line
		s/<\/?FONT[^>]*>//g;
		s/<[^>]*>//g;
		s/^\s+//;
		@data=split(/\s+/,$_);
		($RA,$DEC,$MAJX,$MINX,$THETA,
		 $JMAG,$JMAGU,$HMAG,$HMAGU,$KMAG,$KMAGU)=@data[2..12];
                ++$i;
		$rdist{$i}=$data[1]*60;
		$mags{$i}="$JMAG,$JMAGU,$HMAG,$HMAGU,$KMAG,$KMAGU";
		$RAU=max(abs($MAJX*sin($THETA*$DTOR)),
			 abs($MINX*cos($THETA*$DTOR)));
		$DECU=max(abs($MAJX*cos($THETA*$DTOR)),
			  abs($MINX*sin($THETA*$DTOR)));
		$dataline{$i}=sprintf("%s %s %.3f %.3f",$RA,$DEC,$RAU,$DECU);
	    }

	}elsif ($catalog eq "TYCHO") {
	    return -1 if (/No object found/i);
	    if (/VizieR-5/) {
		# data line
		s/<\/?FONT[^>]*>//g;
		s/<[^>]*>//g;
		s/^\s+//;
		@data=split(/\s+/,$_);
		($RA,$DEC,$PMRA,$PMDEC,$RAU,$DECU,$PMRAU,$PMDECU,
		 $BMAG,$BMAGU,$VMAG,$VMAGU)=@data[2..13];
                ++$i;
		if (abs($PMRA)+abs($PMDEC)>0) {
		    $rapmm=($epoch-2000)/cos($DEC*$DTOR);
		    $decpmm=($epoch-2000);
		    $RA=$RA+$PMRA*$rapmm/3.6e6;
		    $DEC=$DEC+$PMDEC*$decpmm/3.6e6;
		    $RAU=sqrt($RAU**2 + ($PMRAU*$rapmm)**2);
		    $DECU=sqrt($DECU**2 + ($PMDECU*$decpmm)**2);
		}
		$mags{$i}="$BMAG,$BMAGU,$VMAG,$VMAGU";
		$rdist{$i}=$data[1]*60;
		$RAU*=1e-3;
		$DECU*=1e-3;
		$dataline{$i}=sprintf("%s %s %.3f %.3f",$RA,$DEC,$RAU,$DECU);
	    }

	}elsif ($catalog eq "GSC") {
	    if (!(/---/ || /GSC2ID/) && /\d+/) {
		# data line
		@data=split(/\s+/,$_);
		($RA,$DEC,$RAU,$DECU,$GSCEP,$PMRA,$PMDEC,$PMRAU,$PMDECU,
		 $FMAG,$FMAGU,$JMAG,$JMAGU,$VMAG,$VMAGU,$NMAG,$NMAGU)=
		    @data[1..17];
		# apply proper motions
		if (abs($PMRA)+abs($PMDEC)>0) {
		    $rapmm=($epoch-$GSCEP)/cos($DEC*$DTOR);
		    $decpmm=($epoch-$GSCEP);
		    $RA=$RA+$PMRA*$rapmm/3.6e6;
		    $DEC=$DEC+$PMDEC*$decpmm/3.6e6;
		    $RAU=sqrt($RAU**2 + ($PMRAU*$rapmm*1e-3)**2);
		    $DECU=sqrt($DECU**2 + ($PMDECU*$decpmm*1e-3)**2);
		}
		$dra=($RA-$rad)*cos($decd*$DTOR);
		$ddec=($DEC-$decd);
		++$i;
		$mags{$i}="$FMAG,$FMAGU,$JMAG,$JMAGU,$VMAG,$VMAGU,$NMAG,$NMAGU";
		$rdist{$i}=sqrt(($dra**2)+$ddec**2)*3600;
		$dataline{$i}=sprintf("%s %s %.3f %.3f",$RA,$DEC,$RAU,$DECU);
	    }
	}
    }
    my @lines=(1..$i);
    @lines = sort { $rdist{$a} <=> $rdist{$b} } @lines;
    $ndig=1+int(log(0.1+@lines)/log(10));
    my $n=0;
    @magkeys=split(//,$MAGKEYS{$catalog});
    foreach $i (@lines) {
	++$n;
	$lab=sprintf("%s-%0${ndig}d",$catalog,$n);
	@mags=split(/,/,$mags{$i});
	# Catalog line
	$catline="";
	$catline=sprintf("%${ndig}d ",$n) if $addnums;
	$catline.=sprintf("%s-%0${ndig}d %s %.2f",
			  $catalog,$n,$dataline{$i},$rdist{$i});
	$catline.=" @mags" if $addmags;
	print OUT "$catline\n";
	# Region file line
	if ($regfile) {
	    @d=split(/\s+/,$dataline{$i});
	    $regline="fk5;ellipse($d[0],$d[1],$d[2]\",$d[3]\",0) \# text=\{$lab\} color=$color";
	    for ($ii=0;$ii<@magkeys;++$ii) {
		$jj=2*$ii;
		$regline.=" ${magkeys[$ii]}MAG=\{$mags[$jj]\}";
		++$jj;
		$regline.=" ${magkeys[$ii]}MAGU=\{$mags[$jj]\}";
	    }
	    print REG "$regline\n";
	}
    }
    close(OUT);
    close(REG) if ($regfile);

    return (-s $out > 0);
    }

sub printout {
    my $inp=shift;
    print STDERR $inp;
}

sub catline {
}

sub regline {
    my $ii,$jj;
}

sub killoff {
    my $message=shift;
    my $line=0;
    $line=shift if ($#_ >= 0);
    my $s;

    $s=$message;
    $s.= " at " . __FILE__ . " line " . $line . ".\n" if ($message !~ /\n$/ && $line > 0);    
    die $s;
}


__END__

