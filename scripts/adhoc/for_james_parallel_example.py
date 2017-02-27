import multiprocessing as mp
import time


def function_to_parallelise(*args, **kwargs):
    """
    Here we have the main function - in your case running the classification on a single image.
    Your input arguments might be the image filename or data and the trained classifier?
    :param args:
    :return:
    """
    time.sleep(3)
    return "All done, stupid head."


def the_main_code_body():
    """
    Here is the main body of code. Are you running windows? If so, I think it has to be defined in a function, not just
    sitting in a script. On Linux (and Mac?) it doesn't matter.
    :return:
    """
    # this controls the number of processes (cores) to use - None means all.
    # try setting it to 1 - the time takes should be ~ the time to do it in serial
    n_proc = None

    # we can kill jobs after a certain amount of time has elapsed
    max_time_secs = 20
    poll_every = 1  # secs between checking up on jobs

    my_big_array_of_inputs = [
        'foo', 'bar', 'gilbey', 'gabs'
    ]

    # create a pool of workers, by default using all available cores
    # you can specify processes=n if you want to do something different
    pool = mp.Pool(processes=n_proc)

    # here we're going to keep track of the jobs
    # I usually make this a dictionary and use an index that is handy for tracking (e.g. image number)
    jobs = {}

    # we can setup kwargs now if necessary
    kwds = {
        'pointless_kwd1': 'Ben',
        'pointless_kwd2': 'Folds',
    }

    # now loop time
    print "Start loop"
    tic = time.time()

    for i, inp in enumerate(my_big_array_of_inputs):
        # do some prep if you need to, but this will NOT be parallel
        j = pool.apply_async(function_to_parallelise, args=[inp], kwds=kwds)
        # we'll store a pointer to the task AND the start time
        jobs[i] = (j, time.time())

    # we should arrive here almost instantly - but the processes are going on in a different thread now
    print "Submitted all the parallel jobs in %.2f secs" % (time.time() - tic)

    # let's tell the pool that it shouldn't accept any more jobs
    pool.close()

    tic_all = time.time()

    # now we need to gather up the results, waiting for them to be ready
    my_results = []
    while len(jobs):
        for i, (j, tic) in jobs.items():
            if j.ready():
                # a job is 'ready' either because it has finished successfully...
                if j.successful():
                    print "Job %d finished successfully in %.2f secs. Rejoice!" % (i, time.time() - tic)
                    # in practice, you may want to do something more sophisticated when gathering results
                    # i'm just going to add it to an array
                    my_results.append(j.get())
                    jobs.pop(i)
                else:
                    time_running = (time.time() - tic)
                    if time_running > max_time_secs:
                        print "TOO SLOW - we're going to kill job %d" % i
                        # you probably want to log this more carefully
                        jobs.pop(i)
        time.sleep(poll_every)

    # we fall out of the loop when all jobs are finished

    print "Completed all jobs in %.2f seconds" % (time.time() - tic_all)

if __name__ == "__main__":
    the_main_code_body()