#!/usr/bin/env python
###############################################################################
#
# threadpool.py
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

import sys
import threading
import Queue
import traceback
import signal

class ThreadPool:
    POLLING_INTERVAL = 0.100

    END_EMPTY = 0
    END_SIGNALLED = 1

    def __init__(self, num_workers = 10, end_condition = END_EMPTY):
        self.num_workers = num_workers
        self.end_condition = end_condition

        self.finished_event = threading.Event()
        self.request_queue = Queue.Queue()
        self.results_queue = Queue.Queue()

        self.finished_event.clear()

        self.workers = []
        self.requests_pending = []

        self.make_workers(num_workers)

    def make_workers(self, num_workers):
        for i in range(0, num_workers):
            worker = WorkerThread(self.request_queue, self.results_queue)
            self.workers.append(worker)
            worker.start()

    def exc_callback(self, request, exc_info):
        pass

    def make_request(self, func = None, fargs = [], fkwargs = {}, callback = None, cargs = [], ckwargs = {}, exc_callback = None, exc_cargs = [], exc_ckwargs = {}, request_obj = None):
        if not request_obj:
            request = WorkRequest(func, fargs, fkwargs, callback, cargs, ckwargs, exc_callback, exc_cargs, exc_ckwargs)
        else:
            request = request_obj
            
        self.requests_pending.append(request)
        self.request_queue.put(request)

    def set_finished(self):
        self.finished_event.set()

    def wait(self):
        try:
            if self.end_condition == ThreadPool.END_EMPTY:
                while self.requests_pending:
                    # print >> sys.stderr, "More requests pending..."
                    request, result = self.results_queue.get()
                    self.requests_pending.remove(request)

                    if request.exception:
                        print "Exception detected: %s" % str(request.exception)
                        traceback.print_tb(result[2])

            elif self.end_condition == ThreadPool.END_SIGNALLED:
                # print >> sys.stderr, "Waiting for event object signal"
                self.finished_event.wait()

        except Exception, e:
            print >> sys.stderr, "Exception while waiting for thread pool end condition: %s" % e
            traceback.print_tb(sys.exc_info()[2])

        except KeyboardInterrupt, ki:
            print >> sys.stderr, "User keyboard interrupt caught, terminating thread pool."

        finally:
            # print >> sys.stderr, "End condition satisfied, dismissing workers"
            for worker in self.workers:
                worker.dismiss()

            for worker in self.workers:
                worker.join()
            
            # print >> sys.stderr, "Dismissed workers, all done."

class WorkRequest:
    def __init__(self, func, fargs = [], fkwargs = {}, callback = None, cargs = [], ckwargs = {}, exc_callback = None, exc_cargs = [], exc_ckwargs = {}): 

        self.func = func
        self.fargs = fargs
        self.fkwargs = fkwargs
        
        self.callback = callback
        self.cargs = cargs
        self.ckwargs = ckwargs

        self.exc_callback = exc_callback
        self.exc_cargs = exc_cargs
        self.exc_ckwargs = exc_ckwargs

        self.exception = None

class WorkerThread(threading.Thread):
    GET_TIMEOUT = 0.1

    def __init__(self, request_queue, results_queue):
        threading.Thread.__init__(self)
        self.request_queue = request_queue
        self.results_queue = results_queue

        self.dismissed = threading.Event()
        self.dismissed.clear()

    def run(self):
        while not self.dismissed.is_set():
            try:
                # print >> sys.stderr, "Worker is attempting to get from requests queue..."
                request = self.request_queue.get(True, WorkerThread.GET_TIMEOUT)
            except Queue.Empty:
                continue

            try:
                if request.func:
                    result = request.func(*request.fargs, **request.fkwargs)
                    self.results_queue.put((request, result))

                if request.callback:
                    request.callback(*(request.cargs + [result]), **(request.ckwargs))

            except Exception, exc:
                exc_info = sys.exc_info()
                request.exception = exc
                self.results_queue.put((request, exc_info))
                request.exc_callback(*(request.exc_cargs + [request, exc, exc_info]), **request.exc_ckwargs)

    def dismiss(self):
        # print >> sys.stderr, "Worker received dismiss signal..."
        self.dismissed.set()


