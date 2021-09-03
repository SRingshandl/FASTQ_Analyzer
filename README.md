# FASTQ_Analyzer

<img src="https://github.com/SRingshandl/FASTQ_Analyzer/blob/main/Example_Image.png" width="100%">

Final Unix/Perl Course Project 2016.  
Reads a FASTQ file, trims and filters sequences according to user specification. Additional graphics output!  

Besides graphics two files are generated:  
Trimmed contains the FASTQ sequences that correspond to the entered filter sprecifications.  
Statistics contains the Name and the corresponding statistics of the filtered sequences.  
For filter options see below.

Usage:  
Options -f,-m,-t,-l have to be given.  
-f      Path to the Fastq file.  
-o      Path where to write the output (Without the programm will write into a subfolder of the Fastq file).  
-m      Select minimal mean quality for reads (Range: 0-40).  
-t      Select cutoff to trim reads (Range: 0-40).  
-l      Select minimal length of reads after trimming (Range: >=0).  
-s      Give sequence to search for. Write ^ before or $ after sequence to search only at the beginning or end of a read.  
-d      Write -d to perform a direct sequence search without interpretation of the letter code.  
-r      Give sequence that is removed BEFORE subsequent steps. Write ^ before or $ after sequence to search only at the beginning or end of a read.  

Usage example on command line (UNIX systems):  
./Fastq_analysis.pl -f 10k_reads.fastq -m 20 -t 35 -l 200

Usage example on command line (on Windows cmd with installed Perl 5):  
perl Fastq_analysis.pl -f 10k_reads.fastq -m 20 -t 35 -l 200

**The software was created on Mac OS X and optimization for usage on Windows systems was sparsely conducted.** 
**It runs and output is properly generated but output folders are not timestamped, naming of output files is suboptimal and other small issues might appear.**
