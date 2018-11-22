from subprocess import call
import sys
import logging

def run_cmd(name, cmd):
    """ Run a command and stdout print in console.
    ::
        from baseq.mgt import run_cmd
        run_cmd("list files", "ls -l")
    """
    print("[run] Command: {}".format(cmd))
    try:
        call(cmd, shell=True)
        print("[info] {} complete.".format(name, cmd))
    except:
        sys.exit("[error] {} exit with error.".format(name))