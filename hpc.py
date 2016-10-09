import subprocess
import re

# Available clusters on the Stevin infrastructure
clusters = ['raichu', 'delcatty', 'golett']

def running():
    '''
    Get all running calculations on the current cluster.

    Return:
        {'job_name':{job information}}
    '''

    calculations = dict()

    qstat = subprocess.Popen(['qstat'], stdout=subprocess.PIPE, shell=True)
    qstatout, qstaterr = qstat.communicate()
    jobids = re.findall('([0-9]*)\.master', qstatout)

    for jobid in jobids:

        # Retrieve the output of qstat -f $JOB_ID for every job.
        qstatf = subprocess.Popen(['qstat -f %s' % jobid],
                stdout=subprocess.PIPE, shell=True)
        qstatf_output, qstatf_error = qstatf.communicate()

        calculation = dict()
        # Extract various options from qstat -f $JOB_ID output.
        for option in (('job_name', 'Job_Name'),
                       ('job_state', 'job_state'),
                       ('walltime_used', 'resources_used.walltime'),
                       ('walltime_requested', 'Resource_List.walltime'),
                       ('walltime_remaining', 'Walltime.Remaining')):
            try:
                calculation[option[0]] = re.search(option[1]
                                                   + ' = (.*)',
                                                   qstatf_output).group(1)
            except:
                calculation[option[0]] = ''

        calculation['id'] = jobid
        calculations[calculation['job_name']] = calculation

    # Return a dictionary of running calculations
    # {'job_name':{job information}}
    return calculations




def runningold(clusters=clusters, current=False):
    '''
    Get all running calculations.

    Arguments:
        clusters:       clusters to check for running calculations
                        (default: ['haunter', 'gastly', 'raichu', 'delcatty', 'golett'])
        current:        only check current cluster
    Return:
        {'job_name':{job information}}
    '''

    calculations = dict()

    if current:
        clusters = [subprocess.Popen(['echo $VSC_INSTITUTE_CLUSTER'], stdout=subprocess.PIPE, shell=True).communicate()[0].strip()]

    print clusters

    for cluster in clusters:

        # Retrieve job numbers for every cluster.
        qstat = subprocess.Popen(['module swap cluster/%s  > /dev/null 2>&1 ; qstat' % cluster],
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
                           ('walltime_requested', 'Resource_List.walltime'),
                           ('walltime_remaining', 'Walltime.Remaining')):
                try:
                    calculation[option[0]] = re.search(option[1]
                                                       + ' = (.*)',
                                                       qstatf_output).group(1)
                except:
                    calculation[option[0]] = ''

            calculation['cluster'] = cluster
            calculation['id'] = jobid
            calculations[calculation['job_name']] = calculation

    # Return a dictionary of running calculations
    # {'job_name':{job information}}
    return calculations


def set_status(*logfiles):
    r = running()
    for lf in logfiles:
        if lf.file in r:
            lf.hpc = r[lf.file]['job_state']
            lf.cluster = r[lf.file]['cluster']
            lf.jobid = r[lf.file]['id']
