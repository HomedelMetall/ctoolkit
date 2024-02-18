#from toolkit.tools.randomTools import *
from ctoolkit.tools import tools
#tk = tools.tools()
#tools()
#tools = tk.tools.tools()
#tools
# edit this to win.
import functools
def p2f(*args2, **kwargs2):
    def decorator(func):
        @functools.wraps(func)
        # added arguments inside the inner1,
        # if function takes any arguments,
        # can be added like this.
        def inner(*args, **kwargs):

            val = func(*args, **kwargs)
            tk = tools.tools()
            newargs2 = args2 + ([val],)
            tk.print_to_file(*newargs2, **kwargs2)
            #print(*args2, **kwargs2)
            return val

        return inner
    return decorator

