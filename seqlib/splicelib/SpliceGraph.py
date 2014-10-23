

class SpliceGraph(object):

    def __init__(self):
        
        """
          5'      3'
        + GT------AG
        
          3'      5'
        - CT------AC
        """
        self.F_3p_5p_ss = {}
        self.F_5p_3p_ss = {}

        self.R_3p_5p_ss = {}
        self.R_5p_3p_ss = {}
        

