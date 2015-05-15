import subprocess
import re
import fnmatch

# Available clusters on the Stevin infrastructure
clusters = ['raichu', 'delcatty', 'haunter', 'gastly']


def running(clusters=clusters):
    '''
    Get all running calculations.

    Arguments:
        clusters:       clusters to check for running calculations
                        (default: ['haunter', 'gastly', 'raichu', 'delcatty'])
    Return:
        {'job_name':{job information}}
    '''

    running_calculations = {}

    for cluster in clusters:

        # Retrieve job numbers for every cluster.
        qstat = subprocess.Popen(['module swap cluster/%s ; qstat' % cluster],
                                 stdout=subprocess.PIPE, shell=True)
        qstat_output, qstat_error = qstat.communicate()
        jobnumbers = re.findall('([0-9]*)\.master', qstat_output)

        for jobnumber in jobnumbers:

            # Retrieve the output of qstat -f $JOB_ID for every job.
            qstatf = subprocess.Popen(['module swap cluster/%s ; qstat -f %s'
                                       % (cluster, jobnumber)],
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
                try:
                    # First node on which the calculation is running.
                    calculation['node'] = re.search(r'exec_host = ([^/]*)',
                                                    qstatf_output).group(1)
                except:
                    pass
                try:
                    # Retrieve path.
                    start = qstatf_output.find('Output_Path')
                    start = qstatf_output.find(':', start) + 1
                    end = qstatf_output.find('Priority')
                    end = qstatf_output.rfind('/', 0, end)
                    path = ''.join([x.strip() for x
                                    in qstatf_output[start:end].split('\n')])
                    calculation['path'] = path
                except:
                    pass

            calculation['cluster'] = cluster
            running_calculations[calculation['job_name']] = calculation

    # Return a dictionary of running calculations
    # {'job_name':{job information}}
    return running_calculations

def set_status(*logfiles):
    r = running()
    for lf in logfiles:
        if lf.file in r:
            lf.hpc = r[lf.file]['job_state']
