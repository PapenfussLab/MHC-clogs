#!/usr/bin/env python
###############################################################################
#
# execution_engine.py
#
# Copyright (C) 2011 Christopher Davoren
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################
""" 
Job execution engine - schedules and executes jobs based on a given job dependency graph 
"""

import sys
import os
import threading
import multiprocessing
import subprocess
import time
import logging
import traceback

from graph import Graph
from job import Job
from log import get_logger
from log import get_output_logger
from log import get_multiprocess_logger
from utilities import which
from threadpool import ThreadPool
from threadpool import WorkRequest
from exceptions import NoHeadNodesError
from exceptions import NonZeroExitCodeError
from exceptions import MissingInputFileError
from config import get_configuration

def detect_number_of_processors():
    """
    Returns the number of available processors on the system.

    :returns: The number of available processors on the current system.
    """
    return multiprocessing.cpu_count()

class ExecutionEngine:
    """
    The most basic execution engine - takes a job graph and executes one job at a time until all jobs are complete."
    """
    def __init__(self, job_graph=None, terminate_on_nonzero_exit=True, stop_all_on_failure=True):
        self.completed_jobs = []
        self.logger = self.get_logger()
        self.set_graph(job_graph)
        self.terminate_on_nonzero_exit = terminate_on_nonzero_exit
        self.stop_all_on_failure = stop_all_on_failure

        self.verbose = 1
        self.max_parallel = detect_number_of_processors()

    def set_graph(self, job_graph):
        """
        Sets the job graph that the execution should use when executing.

        :param job_graph: The job graph that should be used.
        """
        self.graph = job_graph

    @staticmethod
    def get_logger():
        """
        Gets the debug logger used by all execution engines.

        :returns: The debug logger used by all execution engines.
        """
        return get_logger()

    def set_max_parallel(self, max_parallel="autodetect", buffer=0):
        if isinstance(max_parallel, str) and max_parallel == "autodetect":
            self.max_parallel = detect_number_of_processors()
        else:
            self.max_parallel = max_parallel

        self.max_parallel -= buffer

    def set_max_parallel_buffer(self, buffer):
        self.max_parallel = max(0, self.max_parallel - buffer)

    def _output(self, msg, verbose=1):
        if verbose >= self.verbose:
            # print >> sys.stderr, msg
            get_output_logger().info(msg)

    def _output_job_complete(self, job):
        self._output("Job [%s]: Complete." % str(job))

    def _output_job_executing(self, job):
        self._output("Job [%s]: Executing..." % str(job))

    def _output_job_skipped(self, job, reason):
        self._output("Job [%s]: Skipping, %s." % (str(job), reason))

    def execute(self):
        """
        Executes the jobs in the execution engine's job graph, honouring its dependencies.  Jobs are only executed if they are "out-of-date".
        """
        if not self.graph:
            raise NoJobGraphError("The job graph has not been set")

        if len(self.graph.nodes) == 0:
            return

        if self.graph.is_cyclic():
            raise CyclicDependencyError("The job dependency graph contains cyclic dependencies")

        job_queue = self.graph.get_head_nodes()

        if len(job_queue) == 0:
            raise NoHeadNodesError("The job dependency graph has no head nodes - unable to determine starting jobs")

        while len(job_queue) > 0:
            job = job_queue[0]
            if not job in self.completed_jobs:
                if not job.check_inputs():
                    self.logger.error("Missing input files for job '%s', terminating." % str(job))
                    return False

                if job.out_of_date():
                    print >> sys.stdout, "\tEXECUTING [Job '%s'] (%s)" % (str(job), job.command)
                    self.logger.info("Executing job " + str(job))
                    exit_status = job.execute()
                    
                    if exit_status != 0:
                        if self.terminate_on_nonzero_exit:
                            self.logger.error("Job " + str(job) + " returned a non-zero exit code")
                            return False
                        else:
                            self.logger.warning("Job " + str(job) + " returned a non-zero exit code")

                    self.logger.info("Job completed: " + str(job))
                else:
                    pass
                    # self.logger.info("Job " + str(job) + " is not out of date, skipping")
                    print >> sys.stdout, "\tSKIPPED: NOT OUT OF DATE [Job '%s'] (%s) " % (str(job), job.command)

                self.completed_jobs.append(job)

            job_queue.remove(job)

            self._append_available_jobs(job_queue)

        return True

    def _append_available_jobs(self, job_list):
        for job in self.graph.nodes:
            if job in job_list or job in self.completed_jobs:
                continue

            if len(self.graph.node_parents[job]) == 0:
                job_list.append(job)
                continue

            all_parents_completed = True
            for parent_job in self.graph.node_parents[job]:
                if not parent_job in self.completed_jobs:
                    self.logger.debug("Found uncompleted parent for job " + str(job) + " (" + str(parent_job) + ") - job not added to queue")
                    all_parents_completed = False
                    break

            if all_parents_completed:
                self.logger.info("New job available: " + str(job))
                job_list.append(job)

    def _remove_dependent_jobs(self, job, job_graph, remove_self=True):
        for child in job_graph.node_children[job]:
            self._remove_dependent_jobs(child, job_graph, remove_self=True)

        if remove_self:
            job_graph.remove_node(job)


class MultiprocessExecutionEngine(ExecutionEngine):
    COMPLETION_POLLING_INTERVAL = 0.1

    def __init__(self, job_graph=None, max_processes=-1, terminate_on_nonzero_exit=True):
        ExecutionEngine.__init__(self, job_graph, terminate_on_nonzero_exit)

        if max_processes > 0:
            self.max_processes = max_processes
        elif max_processes != 0:
            nprocs = detect_number_of_processors
            if nprocs > 0:
                self.max_processes = nprocs
            else:
                self.max_processes = 1

        self.active_processes = 0
        self.active_jobs = []
        self.scheduler_lock = threading.Lock()
        self.error_encountered = False
        self.error_job = None

    def _job_process(self, job):
        # job.logger = logging.getLogger("massivetools.JobMultiprocess")
        job.logger = get_multiprocess_logger()
        self.logger = get_multiprocess_logger()

        try:
            exit_status = job.execute()
            sys.exit(exit_status)
        except Exception, e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print >> sys.stderr, "Exception encountered in process for job " + str(job) + " (" + exc_type.__name__ + "): " + str(e) + "\n" + "\n".join(traceback.format_tb(exc_traceback))
            sys.exit(1)

    def _job_thread(self, job):
        try:
            if job.out_of_date():
                job_process = multiprocessing.Process(target=self.__class__._job_process, args=[self,job])
                # job_process.daemon = False
                job_process.start()
                job_process.join()

                exit_status = job_process.exitcode

                if exit_status != 0:
                    if self.terminate_on_nonzero_exit:
                        self.logger.error("Job " + str(job) + " returned a non-zero exit code")
                        self.logger.debug("_job_thread: acquiring scheduler lock")
                        self.scheduler_lock.acquire()
                        self.logger.debug("_job_thread: scheduler lock acquired")
                        self.error_encountered = True
                        self.error_job = job
                        self.logger.debug("_job_thread: releasing scheduler lock")
                        self.scheduler_lock.release()
                        
                        self._schedule()
                        return
                    else:
                        self.logger.warning("Job " + str(job) + " returned a non-zero exit code")

                self.logger.info("Job completed: " + str(job))
            else:
                self.logger.info("Job " + str(job) + " is not out of date, skipping")
        except Exception, e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            self.logger.debug("Exception encountered in starter thread for job " + str(job) + " (" + exc_type.__name__ + "): " + str(e) + "\n" + "\n".join(traceback.format_tb(exc_traceback)))

            self.logger.debug("_job_thread: acquiring scheduler lock")
            self.scheduler_lock.acquire()
            self.logger.debug("_job_thread: scheduler lock acquired")
            self.error_encountered = True
            self.error_job = job
            self.logger.debug("_job_thread: releasing scheduler lock")
            self.scheduler_lock.release()
            
            self._schedule()
            return


        self.logger.debug("_job_thread: acquiring scheduler lock")
        self.scheduler_lock.acquire()
        self.logger.debug("_job_thread: scheduler lock acquired")

        self.job_queue.remove(job)
        self.active_jobs.remove(job)

        self.completed_jobs.append(job)
        self.active_processes -= 1

        self.logger.info("Number of active processes is now " + str(self.active_processes))
        self.logger.debug("_job_thread: releasing scheduler lock")
        self.scheduler_lock.release()

        self._schedule() 

    def _start_job(self, job):
        job_thread = threading.Thread(target=MultiprocessExecutionEngine._job_thread, args=(self, job))
        job_thread.start()

    def _schedule(self):
        self.logger.debug("_schedule: acquiring scheduler lock")
        self.scheduler_lock.acquire()
        self.logger.debug("_schedule: scheduler lock acquired")
        self._append_available_jobs(self.job_queue)

        if self.finished:
            self.logger.info("Nothing to schedule - already finished")
            self.logger.debug("_schedule: releasing scheduler lock")
            self.scheduler_lock.release()
            return

        if self.error_encountered:
            self.logger.error("Terminating after encountering a non-zero exit code from process for job " + str(self.error_job))
            self.logger.debug("_schedule: releasing scheduler lock")
            self.scheduler_lock.release()
            self.finished = True
            return

        for job in self.job_queue:
            if self.max_processes > 0 and self.active_processes >= self.max_processes:
                break

            if job in self.active_jobs:
                continue

            self.active_jobs.append(job)
            self.logger.info("Starting process for job " + str(job))
            self._start_job(job)
            self.logger.info("Process for job " + str(job) + " started")
            self.active_processes += 1
            self.logger.info("Number of active processes is now " + str(self.active_processes))

        if len(self.job_queue) == 0 and len(self.active_jobs) == 0:
            self.finished = True

        self.logger.debug("_schedule: releasing scheduler lock")
        self.scheduler_lock.release()

    def execute(self):
        if not self.graph:
            raise NoJobGraphError("No job graph has been set")

        if len(self.graph.nodes) == 0:
            self.logger.debug("Job graph contains no nodes, exiting")
            return

        if self.graph.is_cyclic():
            raise CyclicDependencyError("The job dependency graph contains cyclic dependencies")

        self.finished = False
        self.job_queue = []
        self._schedule()

        # self.logger.debug("execute: acquring scheduler lock")
        self.scheduler_lock.acquire()
        while not self.finished:
            # self.logger.debug("execute: scheduler lock acquired")
            # self.logger.debug("execute: releasing scheduler lock")
            self.scheduler_lock.release()
            time.sleep(MultiprocessExecutionEngine.COMPLETION_POLLING_INTERVAL)
            # self.logger.debug("execute: acquiring scheduler lock")
            self.scheduler_lock.acquire()

        # self.logger.debug("execute: releasing scheduler lock")
        self.scheduler_lock.release()

        return not self.error_encountered


class ThreadPooledExecutionEngine(ExecutionEngine):
    """
    This execution engine spawns jobs from a central thread pool.  The maximum number of parallel jobs can be configured by changing the number of worker threads in the pool.  See the constructor documentation for details
    """

    def __init__(self, job_graph=None, max_parallel="autodetect", terminate_on_nonzero_exit=True, stop_all_on_failure=True):
        """
        Constructs a new execution engine that uses worker threads to execute jobs concurrently.  Workers fetch jobs from a job queue.

        :param job_graph:  The job graph to use during execution.
        :param max_parallel:  The number of worker threads to use.  If this argument is set to the string value "autodetect" (default), then the number of worker threads that will be set to the number of processors available on the host machine.
        :param terminate_on_non_zero_exit:  Indicates whether the execution should terminate if a command job returns a non-zero exit code.  If this occurs, jobs that are already executing jobs will be allowed to finish.
        """
        ExecutionEngine.__init__(self, job_graph, terminate_on_nonzero_exit, stop_all_on_failure)

        if isinstance(max_parallel, str) and max_parallel == "autodetect":
            self.max_parallel = detect_number_of_processors()
        else:
            self.max_parallel = max_parallel

        self._schedule_lock = threading.Lock()

        self.terminated = False
        self.error_encountered = False
        self.results = {}
        self.pool = None

        self.jobs_to_execute = []
        self.jobs_executing = []

    def _execute_job(self, job):
        result = JobResult()

        try:
            if not job.out_of_date():
                self._output_job_skipped(job, "not out-of-date")
                result.exit_code = 0
                return result

            self._output_job_executing(job)
            result.exit_code = job.execute()
        except Exception, e:
            result.exception = e
            result.exc_info = sys.exc_info()

        return result

    def _exception_callback(self, job, request, exc, exc_info):
        error_str = "Exception detected in job '%s':\n" % str(job)
        error_str += str(exc)
        error_str += "\n".join(traceback.format_tb(exc_info[2])) + "\n"
        self.logger.error(error_str)

        result = JobResult()
        result.exception = exc
        result.exc_info = exc
        result.exit_code = None
        self.results[job] = result

        self._handle_error(job)

    def _remove_dependents(self, job):
        for child in self.jobs_to_execute:
            if child in self.graph.node_children[job]:
                self._remove_dependents(child)

        if job in self.jobs_to_execute:
            self.jobs_to_execute.remove(job)

    def _handle_error(self, job):
        self.error_encountered = True
        if self.stop_all_on_failure:
            self._set_finished(True)
        else:
            self._remove_dependents(job)

        self._schedule_lock.acquire()
        self.jobs_executing.remove(job)
        self._schedule_lock.release()

        self._schedule()

    def _completed_callback(self, job, result):
        self.logger.debug("Job '%s' completion callback" % str(job))
        self.results[job] = result

        if not isinstance(result, JobResult):
            self.logger.error("Job '%s' returned a result of unknown type - must be JobResult")
            self._handle_error(job)
            return

        if result.exception is not None:
            error_str = "Exception detected in job '%s':\n" % str(job)
            error_str += str(result.exception)
            error_str += "\n".join(traceback.format_tb(result.exc_info[2]))
            self.logger.error(error_str)
            self._handle_error(job)
            return

        if result.exit_code != 0:
            error_str = "Job '%s' returned non-zero exit code" % str(job)
            self.logger.error(error_str)
            self._handle_error(job)
            return

        self._schedule_lock.acquire()
        self.jobs_executing.remove(job)
        self._output_job_complete(job)
        self.completed_jobs.append(job)
        self._schedule_lock.release()

        self._schedule()

    """
    def _exception_callback(self, job, worker_request, exception, exc_info):
        self._schedule_lock.acquire()
        self.jobs_executing.remove(job)
        self._schedule_lock.release()

        ThreadPooledExecutionEngine._exception_callback(self, job, worker_request, exception, exc_info)
    """

    def _schedule(self):
        self._schedule_lock.acquire()
        self.logger.debug("_scheduler() entered.")

        if self.terminated:
            self.logger.debug("_scheduler() exiting.")
            self._schedule_lock.release()
            return

        self.logger.debug("Current jobs to execute size: %d" % len(self.jobs_to_execute))
        self.logger.debug("Current job queue size (should be zero): %d" % len(self.job_queue))
        self.logger.debug("Current jobs executing: %d" % len(self.jobs_executing))
        for job in list(self.jobs_to_execute):
            self.logger.debug("Checking job '%s'" % str(job))
            if self._job_available(job):
                self.job_queue.append(job)
                self.jobs_to_execute.remove(job)

        # for job in self.job_queue:
            # self.logger.debug("queueing %s" % job.name)

        while len(self.job_queue) > 0:
            job = self.job_queue.pop()
            self.logger.debug("Sending job '%s' to pool" % str(job))

            self.jobs_executing.append(job)

            self.pool.make_request(
                func = self.__class__._execute_job,
                fargs = [self, job], 
                callback = self.__class__._completed_callback,
                cargs = [self, job],
                exc_callback = self.__class__._exception_callback,
                exc_cargs = [self, job])
        
        if not len(self.job_queue) and not len(self.jobs_to_execute) and not len(self.jobs_executing):
            self.logger.debug("All jobs finished, signalling threadpool to end.")
            self._set_finished()

        self.logger.debug("_scheduler() exiting.")
        self._schedule_lock.release()

    def _set_finished(self, error_encountered = False):
        self.pool.set_finished()
        self.terminated = True

    def _job_available(self, job):
        self.logger.debug("Checking availability of job '%s'" % str(job))
        for parent in self.graph.get_parents(job):
            if not parent in self.completed_jobs:
                self.logger.debug("Parent not completed for job '%s' (%s) - not adding to queue" % (str(job), str(parent)))
                return False

        return True

    def execute(self):
        if not self.graph:
            raise NoJobGraphError("The job graph has not been set")

        if self.graph.is_cyclic():
            raise CyclicDependencyError("The job graph contains a cyclic dependency")

        self.job_queue = []

        self.jobs_to_execute = self.graph.topological_sort()

        for i, job in enumerate(self.jobs_to_execute):
            self.logger.debug("%d : %s" % (i, str(job)))

        self.pool = ThreadPool(num_workers = self.max_parallel, end_condition = ThreadPool.END_SIGNALLED)

        self.terminated = False
        self.logger.debug("Calling scheduler for the first time...")
        self._schedule()
        self.logger.debug("Waiting on pool to finish...")
        
        self.pool.wait()
        self.logger.debug("Execute all done.")

        if self.error_encountered:
            self._output("\nErrors were encountered:\n")
            
            i = 1
            for job, result in self.results.items():
                if result.exception:
                    self._output("%d : %s (Exception):" % (i+1, str(job)))
                    self._output(str(result.exception))
                    self._output("\n".join(traceback.format_tb(result.exc_info[2])) + "\n")
                    i += 1
                elif result.exit_code != 0:
                    if result.exit_code == -2:
                        self._output("%d : %s (non-zero exit code: %d) [keyboard interrupt]\n" % (i, str(job), result.exit_code))
                    else:
                        self._output("%d : %s (non-zero exit code: %d)\n" % (i, str(job), result.exit_code))
                    i += 1

            if not self.stop_all_on_failure:
                self._output("NOTICE: The option 'stop_all_on_failure' was set to FALSE, so only dependents of the listed failed jobs were cancelled.  To terminate the whole pipeline when an error is encountered, ensure this option is set to TRUE.")
        else:
            self._output("Execution complete.")

        return not self.error_encountered


class QSubTPExecutionEngine(ThreadPooledExecutionEngine):
    def __init__(self, job_graph=None, qsub_options="-cwd -b y -q all.q", max_parallel="autodetect", terminate_on_nonzero_exit=True, stop_all_on_failure=True):
        ThreadPooledExecutionEngine.__init__(self, job_graph=job_graph, max_parallel=max_parallel, terminate_on_nonzero_exit=terminate_on_nonzero_exit, stop_all_on_failure=stop_all_on_failure)
        if isinstance(max_parallel, str) and max_parallel == "autodetect":
            self.max_parallel = get_configuration()["max_cluster_queue_size"]

        self.qsub_options = qsub_options

    def set_max_parallel(self, max_parallel, buffer=0):
        if isinstance(max_parallel, str) and max_parallel == "autodetect":
            ThreadPooledExecutionEngine.set_max_parallel(self, get_configuration.max_cluster_queue_size, buffer)
        else:
            ThreadPooledExecutionEngine.set_max_parallel(self, max_parallel, buffer)

    def _execute_job(self, job):
        result = JobResult()

        try:
            if not job.out_of_date(comp_in_out=False):
                self._output_job_skipped(job, "not out-of-date")
                result.exit_code = 0
                return result

            for input in job.file_inputs:
                self.logger.debug("Job '%s' checking for existence of input file '%s'" % (str(job), input))

                tries = 0
                try_limit = get_configuration()["qsub_file_tries"]
                try_delay = get_configuration()["qsub_file_try_delay"]
                input_found = False
                while tries < try_limit and not input_found:
                    if not os.path.exists(input) or not os.path.isfile(input):
                        tries += 1
                        if tries < try_limit:
                            self.logger.debug("Attempt to find input file '%s' for job '%s' (try %d/%d) failed, will try again in %.2f seconds..." % (input, str(job), tries, try_limit, try_delay))
                            time.sleep(try_delay)
                        else:
                            self.logger.debug("Attempt to find input file '%s' for job '%s' (try %d/%d) failed, no more retries." % (input, str(job), tries, try_limit))
                    else:
                        self.logger.debug("Input file '%s' for job '%s' found." % (input, str(job)))
                        input_found = True

                if not input_found:
                    raise MissingInputFileError("The required input file '%s' for job '%s' does not exist" % (input, str(job)))

            self._output_job_executing(job)

            if job.type == Job.COMMAND_JOB:
                use_sync = get_configuration()["qsub_use_sync"]

                if use_sync:
                    qsub_command = "qsub -terse -V -sync y -N '%s' %s '%s'" % (str(job), self.qsub_options, job.command)
                    self.logger.debug(qsub_command)
                    job_process = subprocess.Popen(qsub_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    self.logger.debug("qsub pid: " + str(job_process.pid))
                    output, error = job_process.communicate()
                    result.exit_code = job_process.returncode

                    if len(output) > 0:
                        self.logger.debug("Output for qsub job '%s':\n%s" % (str(job), output))

                    if len(error) > 0:
                        self.logger.debug("Error for qsub job '%s':\n%s" % (str(job), error))

                    self.logger.debug("Exit code for qsub job '%s': %d" % (str(job), result.exit_code))
                else:
                    qsub_command = "qsub -terse -V -N '%s' %s '%s'" % (str(job), self.qsub_options, job.command)

                    self.logger.debug(qsub_command)
                    job_process = subprocess.Popen(qsub_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    self.logger.debug("qsub pid: " + str(job_process.pid))
                    output, error = job_process.communicate()
                    exit_code = job_process.returncode

                    if exit_code != 0:
                        result.exit_code = exit_code
                        self.logger.debug("Failed qsub submission for job '%s' with exit code %d\nOutput:\n%s\nError:\n%s\n" % (str(job), exit_code, output, error))

                    else:
                        job_id = int(output)

                        job_in_queue = True

                        poll_interval = get_configuration()["qsub_poll_interval"]

                        while job_in_queue:
                            time.sleep(poll_interval)

                            self.logger.debug("Polling status of job id %d (for job '%s')" % (job_id, str(job)))
                            qstat_check_command = "qstat -j %d" % job_id
                            job_process = subprocess.Popen(qstat_check_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            output, error = job_process.communicate()
                            exit_code = job_process.returncode

                            job_in_queue = exit_code == 0
                        
                        # Use qacct to determine the return code of the finished job
                        # We wait a few seconds so the log file can catch up
                        time.sleep(2.0)
                        qacct_command = "qacct -j %d" % job_id
                        job_process = subprocess.Popen(qacct_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        output, error = job_process.communicate()
                        exit_code = job_process.returncode

                        if exit_code != 0:
                            self.logger.debug("Error retrieving exit code for job '%s': qacct finished with non-zero exit code %d\nOutput:\n%s\nError:\n%s\n" % (str(job), exit_code, output, error))
                            result.exit_code = exit_code

                        else:
                            for line in output.split("\n"):
                                if line.strip().startswith("exit_status"):
                                    result.exit_code = int(line.strip().split()[1])

                            self.logger.debug("Exit code for job '%s': %d" % (str(job), result.exit_code))

            else:
                result.exit_code = job.execute()

        except Exception, e:
            result.exception = e
            result.exc_info = sys.exc_info()
            result.exit_code = 1

        finally:
            return result 


class ThreadedExecutionEngine(ThreadPooledExecutionEngine):
    """
    A threaded execution engine.
    """
    pass


class QSubExecutionEngine(QSubTPExecutionEngine):
    pass


class JobResult:
    def __init__(self):
        self.exit_code = -1
        self.exception = None
        self.exc_info = None
        pass


class ExecutionEngineFactory:

    @staticmethod
    def get_best_engine():
        qsub_path = which("qsub")

        if qsub_path:
            return QSubExecutionEngine()
        else:
            return ThreadPooledExecutionEngine()

