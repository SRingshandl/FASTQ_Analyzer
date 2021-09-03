#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use File::Basename;
use Getopt::Std;
use Data::Dumper;

system("clear");

## Get User Input
##############################################################################################

my %opts;
getopts('f:o:m:t:l:s:dr:',\%opts);
sub usage {
    my $d=$_[0];
    print "Usage:\n\n";
    print "Options -f,-m,-t,-l have to be given.\n";
    print "-f\tPath to the Fastq file.\n";
    print "-o\tPath where to write the output (Without the programm will write into a subfolder of the Fastq file).\n";
    print "-m\tSelect minimal mean quality for reads (Range: 0-40).\n";
    print "-t\tSelect cutoff to trim reads (Range: 0-40).\n";
    print "-l\tSelect minimal length of reads after trimming (Range: >=0).\n";
    print "-s\tGive sequence to search for. Write ^ before or \$ after sequence to search only at the beginning or end of a read.\n";
    print "-d\tWrite -d to perform a direct sequence search without interpretation of the letter code.\n";
    print "-r\tGive sequence that is removed BEFORE subsequent steps. Write ^ before or \$ after sequence to search only at the beginning or end of a read.\n";
    
    die "\n".$d."\n\n";
}
my $Fastqpath=$opts{f};
my $Resultpath=0;
if($opts{o} && -d $opts{o}){$Resultpath=$opts{o}};
unless($Resultpath=~/\/$/ || $Resultpath eq "0"){$Resultpath=$Resultpath."\/"}
my $Readcutoff=$opts{m};
my $Trimmingcutoff=$opts{t};
my $Minimallength=$opts{l};
my $Seqsearch="";
my $Seqremoval="";
if($opts{s}){$Seqsearch=uc($opts{s})};
if($opts{r}){$Seqremoval=uc($opts{r})};

my $first=0;
unless($opts{d}){
    my @Seqsearch;
    my @Seqsearch2=split(//,$Seqsearch);
    foreach my $letter (@Seqsearch2) {
        if ($letter=~/\^/ && $first==0){push(@Seqsearch,"^");$first=1}
        if ($letter=~/\$/){push(@Seqsearch,"\$")}
        if ($letter=~/[ATGC]/){push(@Seqsearch,$letter)}
        if ($letter=~/N/){push(@Seqsearch,"[AGTCN]")}
        if ($letter=~/R/){push(@Seqsearch,"[AGR]")}
        if ($letter=~/Y/){push(@Seqsearch,"[CTY]")}
        if ($letter=~/S/){push(@Seqsearch,"[GCS]")}
        if ($letter=~/W/){push(@Seqsearch,"[ATW]")}
        if ($letter=~/K/){push(@Seqsearch,"[GTK]")}
        if ($letter=~/M/){push(@Seqsearch,"[ACM]")}
        if ($letter=~/B/){push(@Seqsearch,"[CGTB]")}
        if ($letter=~/D/){push(@Seqsearch,"[AGTD]")}
        if ($letter=~/H/){push(@Seqsearch,"[ACTH]")}
        if ($letter=~/V/){push(@Seqsearch,"[ACGV]")}
        if ($letter=~/\.|\-/){push(@Seqsearch,".*")}    
          
    }
    $Seqsearch=join('',@Seqsearch);
}

my ($sec,$min,$hour,$day,$month,$year,$wday,$yday,$isdst) = localtime();
$year+=1900;$month+=1;

unless ($Resultpath) {
    $Resultpath = "Output_";
    if($Fastqpath=~/^\.\//){
        $Fastqpath=cwd()."\/".basename($Fastqpath);
    } elsif($Fastqpath=~/^\.\.\//){
        $Fastqpath=dirname(cwd())."\/".basename($Fastqpath);
    } else {
        unless ($Fastqpath=~/^\/Users\//){
            $Fastqpath=cwd()."\/".$Fastqpath;
        }
    }
    if ($Fastqpath=~/^\/Users\//) { 
        unless(-d dirname($Fastqpath)."\/$year"."\-".$month."\-".$day."\_Results\/") {
            mkdir dirname($Fastqpath)."\/$year"."\-".$month."\-".$day."\_Results\/";
        }
        unless(-d dirname($Fastqpath)."\/$year"."\-".$month."\-".$day."\_Results\/".basename($Fastqpath)."\/") {
            mkdir dirname($Fastqpath)."\/$year"."\-".$month."\-".$day."\_Results\/".basename($Fastqpath)."\/";
        }
        $Resultpath=dirname($Fastqpath)."\/$year"."\-".$month."\-".$day."\_Results\/".basename($Fastqpath)."\/";
    }
}

usage("ERROR: (-f) Fastq file is not valid! It doesn't exist!") unless ($opts{f} && -f $opts{f});
if($opts{o}){
    usage("ERROR: (-o) Output path is not valid! Does the folder already exist?") unless (-d $opts{o});
}
usage("ERROR: (-m) No or invalid minimal mean quality given!") unless ($opts{m} && $opts{m}<=40 && $opts{m}>=0);
usage("ERROR: (-t) No or invalid cutoff for trimming given!") unless ($opts{t} && $opts{t}<=40 && $opts{t}>=0);
usage("ERROR: (-l) No or invalid minimal length of final read given!") unless ($opts{l} && $opts{l}>=0);
if ($opts{s}){
    usage("ERROR: (-s) No valid sequence given! Unknown letter(s)!") unless ($opts{s}=~/^[AGCTNRYSWKMBDHV\.\-^\$]*$/);
}
if ($opts{r}){
    usage("ERROR: (-r) No valid sequence given! Unknown letter(s)!") unless ($opts{r}=~/^[AGCTN^\$]*$/);
}

print "File is: \t\t$Fastqpath\n";
print "Resultpath is: \t\t$Resultpath\n";
print "Minimal mean quality: \t$Readcutoff\n";
print "Trimming-Cutoff: \t$Trimmingcutoff\n";
print "Minimal readlength: \t$Minimallength\n";
if($opts{s}){print "Search sequence: \t$Seqsearch\n"};
if($opts{r}){print "Removal sequence: \t$Seqremoval\n"};

## Read Input and define Lines
##############################################################################################

open(FASTQ, "<", $Fastqpath) || die "Program failed to open the Fastq file. Reason: $!\n";
open(FASTQR, ">>", $Resultpath."Trimmed\_".basename($Fastqpath)) || die "Program failed to generate $Resultpath the Fastq that should contain modified reads. Reason: $!\n";
open(FASTQS, ">>", $Resultpath."Statistics\_".basename($Fastqpath).".txt") || die "Program failed to generate the Fastq that should contain statistics to the reads. Reason: $!\n";
print FASTQS "Name\tLength\tGC content\tQuality Mean\tQuality Median\tQuality Stddev\n";

my ($Line1,$Line2,$Line3,$Line4);
my $Linevar = 1;
my @Qualities;
my @Qualitiesnum;
my $Alldefined = 0;
my $phred=0;
my %Qualitypos;
my $Resultcount=0;
my $Startstring="";

while (<FASTQ>){
    chomp($_);
    next unless ($_);
    
    if ($Linevar==4){
        $Line4=$_; 
        
        @Qualities=split(//,$Line4);
        foreach my $Qualityvalue (@Qualities){
            if(ord($Qualityvalue)<64){
                $phred=33;
                print "Program detected phred score system -33\n";
                last;
            } elsif (ord($Qualityvalue)>73){
                $phred=64;
                print "Program detected phred score system -64\n";
                last;
            }
        }
        $Linevar=0;
    }
    if ($phred !=0){last}
    $Linevar++;
}
close(FASTQ);

if ($phred==0){print "Program couldn't detect phred score. Taking phred 33 for calculation of Qualities.\n";$phred=33}

open(FASTQ, "<", "$Fastqpath") || die "Program failed to open the Fastq file. Reason: $!\n";
$Linevar=1;
while (<FASTQ>){
    
    chomp($_);
    next unless ($_);
    $Alldefined = 0;    
    
    if ($Linevar==1){$Line1=$_}
    if ($Linevar==2){$Line2=$_}
    if ($Linevar==3){$Line3=$_}
    if ($Linevar==4){
        $Line4=$_; 
        $Alldefined=1;
        $Linevar=0;
    }

    if ($Alldefined == 1) {                                 # Start with checks and writing
        
        ##REMOVE PIECES
        if($opts{r}){        
            while ($Line2 =~ m/$Seqremoval/g) {
                $Startstring = pos($Line2);
                $Line2=substr($Line2,0,$Startstring-length($Seqremoval)).substr($Line2,$Startstring);
                $Line4=substr($Line4,0,$Startstring-length($Seqremoval)).substr($Line4,$Startstring);
            }
        }
        ##END OF REMOVAL
        
        @Qualities=split(//,$Line4);
        @Qualitiesnum=();
        foreach my $Qualityvalue (@Qualities){
            push(@Qualitiesnum,ord($Qualityvalue)-$phred);
        }
   
        my $pos=0;
        for($pos=scalar(@Qualitiesnum)-1;$pos>0;$pos--){
            if ($Qualitiesnum[$pos] >= $Trimmingcutoff){last}            #last ends only small loop!
        }
        $Line2=substr($Line2,0,$pos+1);
        $Line4=substr($Line4,0,$pos+1);
        @Qualitiesnum=splice(@Qualitiesnum,0,$pos+1);
        
        if (mean(\@Qualitiesnum)> $Readcutoff && length($Line2)>=$Minimallength){
            if ($Seqsearch) {
                unless ($Line2=~/$Seqsearch/){$Linevar++;next}
            }
                        
            ## Block schreibt Werte für Datei Qpos
            ###########################################################################################################
            my $Position=1;
            foreach my $Qvalue (@Qualitiesnum) {
                if ($Qualitypos{$Position}->{$Qvalue}) {
                    $Qualitypos{$Position}->{$Qvalue}+=1;
                } else {
                    $Qualitypos{$Position}->{$Qvalue}=1;
                }
                $Position++;
            }            
            ###########################################################################################################
            
            print FASTQR $Line1."\n".$Line2."\n".$Line3."\n".$Line4."\n";
            print FASTQS $Line1."\t".length($Line2)."\t".GCcont($Line2)."\t".mean(\@Qualitiesnum)."\t".median(\@Qualitiesnum)."\t".stddev(\@Qualitiesnum)."\n";
            $Resultcount++;
        }
    }
    if ($. % 5000 == 0){print "Currently working on line $.. Please wait!\n"}
    $Linevar++;
}                                               ## This sign ends the while loop!!!

close(FASTQS);
close(FASTQR);
close(FASTQ);


if (-z $Resultpath."Trimmed\_".basename($Fastqpath)) {die "No suitable results found! Aborting and not generating graphical output!\n"}
print "$Resultcount suitable reads found!\n";

##Quantengenerator
##############################################################################################################
foreach my $Qval (0..40) {
    foreach my $Pos (6..9){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{5}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (11..14){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{10}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (16..19){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{15}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (21..24){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{20}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (26..29){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{25}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (31..39){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{30}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (41..49){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{40}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (51..74){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{50}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (76..99){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{75}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (101..149){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{100}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (151..199){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{150}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (201..249){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{200}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (251..299){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{250}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (301..399){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{300}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (401..499){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{400}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (501..599){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{500}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (601..699){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{600}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (701..799){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{700}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (801..899){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{800}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (901..999){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{900}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
    foreach my $Pos (1001..3000){if($Qualitypos{$Pos}->{$Qval}){$Qualitypos{1000}->{$Qval}+=$Qualitypos{$Pos}->{$Qval};delete $Qualitypos{$Pos}}}
}
##############################################################################################################



open(QPOS,">>", $Resultpath."Qpos.txt") || die "Couldn't generate file!\n";
foreach my $Qpos (sort({$a<=>$b} keys(%Qualitypos))) {
    foreach my $Qval (keys(%{$Qualitypos{$Qpos}})) {
     print QPOS $Qpos."\t".$Qval."\t".$Qualitypos{$Qpos}->{$Qval}."\n";
    }
}
close(QPOS);

my $Rpath=$Resultpath;
my $Rstatpath=$Resultpath."Statistics\_".basename($Fastqpath).".txt";
my $RQpospath=$Resultpath."Qpos.txt";

system("Rscript RStatistics.R $Rpath $Rstatpath $RQpospath");                                       #Start R-Script

##DELETE QPOS FILE
##############################################################################################################
while (!-e "Graphical-Results.pdf"){
    print "running";
}
if (-e "Graphical-Results.pdf"){
    unlink $Resultpath."Qpos.txt";
}

##SUBS
##############################################################################################################
sub mean {
    #calculate Sum
    my $sum=0;
    foreach my $num (@{$_[0]}) {$sum+=$num}
    #calculate Mean
    my $Mean=$sum/scalar(@{$_[0]});
    # return values
    return($Mean);
}

sub median {
    my @values = sort ({$a <=> $b} @{$_[0]});
    my $length = scalar(@values);
    if($length %2 == 1){                    # Number of elements is odd
        return $values[int($length/2)];
    } else {                                # Number of elements is even
        return ($values[int($length/2)-1] + $values[int($length/2)])/2;
    }
}

sub stddev {
    #calculate Sum
    my $sum=0;
    foreach my $num (@{$_[0]}) {$sum+=$num}
    #calculate Mean
    my $Mean=$sum/scalar(@{$_[0]});
    #calculate Numerator for standard deviation
    my $numerator=0;
    foreach my $num (@{$_[0]}){$numerator+=($num-$Mean)**2}
    #calculate rest for standard deviation
    my $Std=sqrt($numerator/(scalar(@{$_[0]})-1));
    #return value
    return($Std);
}

sub GCcont {
    my $GCcount = 0;
    my $ATcount = 0;
    my @Bases=split(//,$_[0]);
    foreach my $Base (@Bases){
        if ($Base eq "G" || $Base eq "C"){$GCcount++}
        if ($Base eq "A" || $Base eq "T"){$ATcount++}
    }
    my $GCcont=($GCcount/($GCcount+$ATcount));
    #return value
    return($GCcont);
}