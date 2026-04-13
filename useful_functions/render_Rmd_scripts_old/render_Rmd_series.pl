#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

####### render_Rmd.pl my_script_1.Rmd my_script_2.Rmd
## runs those scripts one after the other

my $series_name = "zzz_Rmd_series"; ## this will be used as the start of the output file names for those files that relate to all Rmd files in the series

my $R_module = "fhR/4.4.1-foss-2023b-R-4.4.1";
my $pandoc_module = "Pandoc/2.13";
my $keep_html = 0;

my $use_sbatch = 1;
my $numThreads = 4;
my $walltime = "3-0";
my $debug = 0;

GetOptions("series=s"      => \$series_name,
           "R=s"           => \$R_module,
           "pandoc=s"      => \$pandoc_module,
           "html=i"        => \$keep_html,
           "t=i"           => \$numThreads,           # '--t 4' to use 4 threads
           "sbatch=i"      => \$use_sbatch,
           "wall=s"        => \$walltime,             # '--wall 0-6' to specify 6 hrs
           "debug"         => \$debug                 # '--debug' to just test
           ) or die "\n\nterminating - unknown option(s) specified on command line\n\n"; 


# module load fhR/4.4.1-foss-2023b-R-4.4.1
# module load Pandoc/2.13
# Rscript -e 'rmarkdown::render("test_randomCodeBits.Rmd", output_format="github_document", clean=TRUE)' > test_randomCodeBits.Rrender.Rout.txt 2> test_randomCodeBits.Rrender.Rerr.txt
# rm test_randomCodeBits.html
# module purge

#### notes: 
# - we get TWO output files for some reason (md and html)
# - running R this way, I do not have /home/jayoung/R/x86_64-pc-linux-gnu-library/4.4 in .libPaths(), whereas when I run it from Hutch RStudio-server I do. I don't know where that comes from - I suspect something to do with Rstudio server setup:
# Rscript -e '.libPaths()'


### test - it works from same directory as the Rmd file:
# cd ~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes/Rscripts/
# ~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes/useful_functions/render_Rmd.pl test_randomCodeBits.Rmd

### test from project top dir - works
# cd ~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes/
# ~/FH_fast_storage/git_more_repos/Rtest_and_Rnotes/useful_functions/render_Rmd.pl Rscripts/test_randomCodeBits.Rmd

#####################
if ($use_sbatch == 1) {print "\n\nUsing sbatch to parallelize\n\n";}

my $shellScript = "$series_name.Rrender.sh";
my $logfile = "$series_name.Rrender.log.txt";

open (SH, "> $shellScript");
open (LOG, "> $logfile");

print SH "#!/bin/bash\n";
print SH "source /app/lmod/lmod/init/profile\n";
print SH "module purge\n";
print SH "module load $R_module\n";
print SH "module load $pandoc_module\n";

foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nterminating - Rmd file $file does not exist\n\n";
    }
    my $fileStem = $file; $fileStem =~ s/\.Rmd$//;
    my $Rout = "$fileStem.Rrender.Rout.txt";
    my $Rerr = "$fileStem.Rrender.Rerr.txt";

    print SH "\n\n";
    print SH "echo 'Running script $file' >> $logfile\n";
    # print SH "Rscript -e 'rmarkdown::render(\"$file\", output_format=\"github_document\", clean=TRUE)' > $Rout 2> $Rerr || echo '    FAILED!' >> $logfile && exit 1\n";
    print SH "Rscript -e 'rmarkdown::render(\"$file\", output_format=\"github_document\", clean=TRUE)' > $Rout 2> $Rerr\n\n";

    print SH "STATUS=\$?\n";
    ## Need double-quote not single to successfully echo a variable
    # print SH "echo \"  exit code was \$STATUS\" >> $logfile\n";

    print SH "if [ \$STATUS -eq 0 ]\n";
    print SH "then\n";
    print SH "  echo '  SUCCEEDED' >> $logfile\n";
    print SH "else\n";
    print SH "  echo '  FAILED' >> $logfile\n";
    print SH "  exit 1\n";
    print SH "fi\n";

    if (!$keep_html) { print SH "rm $fileStem.html 2>&1 >> $logfile\n"; }

}
print SH "\nmodule purge\n";

## print a message to the log file if everything worked
print SH "\necho 'All scripts ran OK' >> $logfile\n\n";

close SH;

my $command;
if ($use_sbatch == 1) {
    $command = "sbatch --job-name=Rmd -t $walltime --cpus-per-task=$numThreads --wrap=\"bash ./$shellScript 2>&1 >> $logfile\"";
} else {
    $command = "bash ./$shellScript 2>&1 >> $logfile";
}
print "command $command\n\n";
print LOG "running command:\n$command\n\n";
close LOG;
if ($debug == 0) { system($command); }
