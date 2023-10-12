# DISCLAIMER
# THIS PROGRAM WILL UPDATE THE VARIABLES THAT DEFINE THE VERSION OF TOOLKIT.
import sys
from datetime import datetime

# RUN WITH:
#   $ python update_version.py version_name
s = 'def print_version():\n'
s += '\tprint("Toolkit v_%s"%tk_version)\n\n'
s += "tk_version = '%s'\n" % (sys.argv[1])
s += "tk_date = '%s'\n" % (datetime.today().strftime('%Y-%m-%d @ %H:%M:%S'))
s += "if __name__=='__main__':\n"
s += "\tprint('Version:', tk_version)\n"
s += "\tprint('Date:', tk_date)\n"
with open("version.py", 'w') as file:
    file.write(s)

print("Version updated!\n%s" % sys.argv[1])
