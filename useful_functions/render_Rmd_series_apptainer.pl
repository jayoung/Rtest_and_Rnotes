#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

####### render_Rmd_series_apptainer.pl my_script_1.Rmd my_script_2.Rmd
## runs those scripts one after the other

## see https://github.com/FredHutch/wiki/blob/r-apptainer-cmdline/_scicomputing/software_R.md#containers

## questions - if I want to use multiple threads inside the apptainer, does it work any differently?

my $series_name = "zzz_Rmd_series"; ## this will be used as the start of the output file names for those files that relate to all Rmd files in the series

my $apptainer_module = "Apptainer/1.1.6";
my $sif = "https://sif-registry.fredhutch.org/bioconductor_docker_RELEASE_3_22-R-4.5.2.sif";

my $keep_html = 0;

my $use_sbatch = 1;
my $numThreads = 1;
my $walltime = "3-0";
my $debug = 0;

GetOptions("series=s"      => \$series_name,
           "apptainer=s"   => \$apptainer_module,
           "sif=s"         => \$sif,
           "html=i"        => \$keep_html,
           "t=i"           => \$numThreads,
           "sbatch=i"      => \$use_sbatch,
           "wall=s"        => \$walltime,
           "debug"         => \$debug # '--debug' to just test
           ) or die "\n\nterminating - unknown option(s) specified on command line\n\n"; 


if ($numThreads > 1) {
    die "\n\nTerminating - don't know how to make it work with >1 CPU yet\n\n";
}


#####################
if ($use_sbatch == 1) {print "\n\nUsing sbatch to parallelize\n\n";}

## we'll have one script for the actual R commands and another for the external apptainer wrapper
my $shellScript = "$series_name.Rrender.Rcode.sh";
my $apptainer_wrapper = "$series_name.Rrender.apptainerWrap.sh";
my $logfile = "$series_name.Rrender.log.txt";

open (LOG, "> $logfile");


#### write the shell script containing the R commands
open (SH, "> $shellScript");
print SH  "#!/bin/bash\n";
foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nterminating - Rmd file $file does not exist\n\n";
    }
    my $fileStem = $file; $fileStem =~ s/\.Rmd$//;
    my $Rout = "$fileStem.Rrender.Rout.txt";
    my $Rerr = "$fileStem.Rrender.Rerr.txt";

    print SH "\n\n";
    print SH "echo 'Running script $file' >> $logfile\n";

    print SH "    Rscript -e 'rmarkdown::render(\"$file\", output_format=\"github_document\", clean=TRUE)' > $Rout 2> $Rerr\n\n";

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

## print a message to the log file if everything worked
print SH "\necho 'All scripts ran OK' >> $logfile\n\n";

close SH;

### write the apptainer shell script that calls the R-command shell script
open (APP, "> $apptainer_wrapper");

print APP "#!/bin/bash\n";
print APP "source /app/lmod/lmod/init/profile\n";
print APP "module purge\n";
print APP "module load $apptainer_module\n";

### figure out .libPaths()
### .libPaths() should be figured out by .Rprofile
## because I am setting this using ~/gizmo_login_items.bash to 
# export R_LIBS_USER=$HOME/malik_lab_shared/linux_gizmo/R_packages/v_fh_4.4.0-foss-2023b
## see also https://github.com/FredHutch/wiki/blob/r-apptainer-cmdline/_scicomputing/software_R.md#notes-about-these-r-containers
## I could instead add this to my ~/.Rprofile - that will set the libPaths appropriately
# source("https://sciwiki.fredhutch.org/assets/apptainer_config.R")
# print APP "unset R_LIBS_USER\n";
# print APP "export R_LIBS_USER=/home/jayoung/R/x86_64-pc-linux-gnu-library/4.5-apptainer\n";

print APP "\necho 'Running apptainer' >> $logfile\n\n";

print APP "apptainer run \\\n";
# print APP "    --cpus $numThreads \\\n";
# FATAL:   container creation failed: while applying cgroups config: rootless cgroups require a D-Bus session - check that XDG_RUNTIME_DIR and DBUS_SESSION_BUS_ADDRESS are set

print APP "    --bind /fh/fast:/fh/fast \\\n";
print APP "    $sif \\\n";
print APP "    bash $shellScript >> $logfile\n";

print APP "\necho 'Finished running apptainer' >> $logfile\n\n";
print APP "\nmodule purge\n";

close APP;

my $command;
if ($use_sbatch == 1) {
    $command = "sbatch --job-name=Rmd -t $walltime --cpus-per-task=$numThreads --wrap=\"bash ./$apptainer_wrapper 2>&1 >> $logfile\"";
} else {
    $command = "bash ./$apptainer_wrapper 2>&1 >> $logfile";
}
print "command $command\n\n";
print LOG "running command:\n$command\n\n";
close LOG;
if ($debug == 0) { system($command); }
