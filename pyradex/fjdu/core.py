try:
    from myradex import myradex_wrapper

class fjdu(object):
    def __init__(self, **kwargs):
        pass

    def load_datfile(self, filename=None):
        if filename is not None:
            self.dpath = os.path.dirname(filename) or self.dpath
            self.fname = os.path.basename(filename)
        myradex_wrapper.config_basic(self.dpath,
                                     self.fname,
                                     True)
        
