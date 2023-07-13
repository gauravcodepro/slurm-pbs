#!/bin/bash
PS3 = "Altair PBS Pro Manual at hand (#)? "; \
select i in qsub qdel qhold qrls qmove qstat qselect qalter; do
if [[ $i == qsub ]]; then
        echo "this is a qsub option and this is used for the 
        submission of the job queue and the job queue should
        follow the following format for the submission to the
        #PBS -o path_out_file
        #PBS -e path_error_file
        #PBS -l walltime=1:00:00
        #PBS -l mem=100GB
        #PBS -l ncpus=64
        #PBS -j jobname
        submission should be made as 
        qsub script_name.sh
        qsub -M user@pbspro.com"
        break; fi; 
if [[ $i == qdel ]]; then
        echo " this is used for the continous monitoring of the 
        job queues on the system and is used for the deletion of
        the jobs from the system and this can be used only by the
        administrator of the system"; 
        break; fi;
if [[ $i == qhold ]]; then
        echo "The qhold is used to hold the jobs on the
        system when the requested resources are not available
        qhold [ -h hold_list ] job_identifier
        n:none, u:user,o:operator,s:system"; 
        break; fi;
if [[ $i == qrls ]]; then
        echo "The qrls is used to release the jobs holded
        due to the system resource unavailability 
        qrls [ -h hold_list ] job_identifier
        qstat -a jobid # get the job
        qhold jobid # hold the job
        qrls -h u jobid # release the job"; 
        break; fi;
if [[ $i == qmove ]]; then
        echo "this allows for the moving of the jobs to
        the new allocated server and to allow for the new
        server to be used.
        qmove destination job_identifier(s)"; 
        break; fi;
if [[ $i == qstat ]]; then
        echo "qstat is the most user friendly command
        to be used on the PBS and it is used to check 
        your server status and also check your job status
        Available options are qstat -a or qstat -Bf or
        qstat -Q or qstat -Qf qstat -f jobname or qstat -u username
        qstat -G (GB mem list), qstat -M (megaword listing),
        qstat -n (display all the nodes), qstat -s (display job comments)
        qstat -q (display limits) "; 
        break; fi;
if [[ $i == qselect ]]; then
        echo "The qselect is used to alter the allocated resources and 
        also to allow for the change. List all the nodes used by the 
        user which require such number of nodes
        qselect -u user -l ncpus.gt.cpu # where cpu is the cpu number
        qstat -a ' qselect -u user -l ncpus.gt.cpu"; 
        break; fi;
if [[ $i == qalter ]]; then
        echo "qalter is used to alter the allocated resources and 
        it should be used to change the resources on the jobs 
        already allocated. 
        qalter -l walltime=110:00 -N bladename username"; 
        break; fi;
done