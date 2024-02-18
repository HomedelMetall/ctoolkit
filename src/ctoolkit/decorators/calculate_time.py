import time
from ctoolkit.global_vars.timers_dict import *

def calculate_time(func):

    # added arguments inside the inner1,
    # if function takes any arguments,
    # can be added like this.
    def inner1(*args, **kwargs):

        # storing time before function execution
        begin = time.time()

        val = func(*args, **kwargs)

        # storing time after function execution
        end = time.time()
        timer_name = func.__name__
        timer_time = end-begin
        if timer_name in timers_dict:
            timers_dict[timer_name] += timer_time
        else:
            timers_dict[timer_name] = timer_time

        return val

    return inner1
