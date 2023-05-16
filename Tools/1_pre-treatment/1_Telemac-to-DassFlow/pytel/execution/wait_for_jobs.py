""" Scripts to get info after cluster submission """
import time
import subprocess


def check_job_slurm(job_id, call_count=0):
    """
    Return a job status

    @param job_id (str) Job id of the job to check
    @param call_count (int) Recursive number of call if 5 is reached returns -1

    @returns (int) 'success' if completed
                   'failed' if job crashed
                   'timeout' if job timed out
                   'running' otherwise (job still running or pending)
    """
    cmd = "sacct -j {} -o JobID,State -P".format(job_id)

    # Scontrol will give us information on the job
    # TODO: Extract more than just the path ? We could get the listing...
    try:
        tmp = subprocess.check_output(cmd, shell=True)
    except Exception as e:
        time.sleep(1)
        tmp = subprocess.check_output(cmd, shell=True)
    output = tmp.decode('utf-8')
    # In case the return is not as it should trying again after 5 tries
    # returning failed
    if len(output.split('\n')) < 3:
        if call_count == 5:
            return 'failed'

        return check_job_slurm(job_id, call_count=call_count+1)

    # Extract line containing status (second one)
    line = output.split('\n')[1]
    # Line should be something like that: 30346693|COMPLETED
    state = line.split('|')[1]

    if state == 'COMPLETED':
        return 'success'
    if state in ['CANCELLED', 'FAILED', 'DEADLINE']:
        return 'failed'
    if state == 'TIMEOUT':
        return 'timeout'

    return 'running'

def get_job_time_slurm(job_id):
    """
    Return the time that took a job

    @param job_id (str) Job id of the job to check

    @returns (float) time in seconds
    """
    cmd = "sacct -j {} -o ElapsedRaw -P".format(job_id)
    try:
        tmp = subprocess.check_output(cmd, shell=True)
    except:
        time.sleep(1)
        tmp = subprocess.check_output(cmd, shell=True)
    output = tmp.decode('utf-8')

    # Extract line containing status (second one)
    line = output.split('\n')[1]
    # Line should be something like that: 666
    run_time = float(line)

    return run_time
