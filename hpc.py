import subprocess
import re

# Available clusters on the Stevin infrastructure
clusters = ['raichu', 'delcatty', 'haunter', 'gastly', 'golett']


def running(clusters=clusters):
    '''
    Get all running calculations.

    Arguments:
        clusters:       clusters to check for running calculations
                        (default: ['haunter', 'gastly', 'raichu', 'delcatty'])
    Return:
        {'job_name':{job information}}
    '''

    calculations = dict()

    for cluster in clusters:

        # Retrieve job numbers for every cluster.
        qstat = subprocess.Popen(['module swap cluster/%s ; qstat' % cluster],
                                 stdout=subprocess.PIPE, shell=True)
        qstatout, qstaterr = qstat.communicate()
        jobids = re.findall('([0-9]*)\.master', qstatout)

        for jobid in jobids:

            # Retrieve the output of qstat -f $JOB_ID for every job.
            qstatf = subprocess.Popen(['module swap cluster/%s ; qstat -f %s'
                                       % (cluster, jobid)],
                                      stdout=subprocess.PIPE, shell=True)
            qstatf_output, qstatf_error = qstatf.communicate()

            calculation = dict()
            # Extract various options from qstat -f $JOB_ID output.
            for option in (('job_name', 'Job_Name'),
                           ('job_state', 'job_state'),
                           ('walltime_used', 'resources_used.walltime'),
                           ('walltime_remaining', 'Walltime.Remaining')):
                try:
                    calculation[option[0]] = re.search(option[1]
                                                       + ' = (.*)',
                                                       qstatf_output).group(1)
                except:
                    pass

            calculation['cluster'] = cluster
            calculations[calculation['job_name']] = calculation

    # Return a dictionary of running calculations
    # {'job_name':{job information}}
    return calculations


def set_status(*logfiles):
    r = running()
    for lf in logfiles:
        if lf.file in r:
            print r[lf.file]['job_state']
            lf.hpc = r[lf.file]['job_state']
