"""
Temporary, until it's merged in astropy: https://github.com/astropy/astropy/pull/2793
"""
def in_ipynb_kernel():
    """
    Determine whether the current code is being executed from within an IPython
    notebook kernel.  If ``True``, the code may be executed from a notebook or
    from a console connected to a notebook kernel, but (we believe) it's not
    being executed in a console that was initialized from the command line.
    """
    try:
        cfg = get_ipython().config
        app = cfg['IPKernelApp']
        # ipython 1.0 console has no 'parent_appname',
        # but ipynb does
        if ('parent_appname' in app and
            app['parent_appname'] == 'ipython-notebook'):
            return True
        else:
            return False
    except NameError:
        # NameError will occur if this is called from python (not ipython)
        return False
