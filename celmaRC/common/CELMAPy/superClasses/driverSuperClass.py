#!/usr/bin/env python

"""
Contains the super class for the drivers
"""

#{{{DriverSuperClass
class DriverSuperClass(object):
    """
    The parent driver class
    """

    #{{{Constructor
    def __init__(self                ,\
                 dmp_folders         ,\
                 useSubProcess = True,\
                 collectPaths  = None,\
                 ):
        #{{{docstring
        """
        This constructor

        * Sets the collect path.
        * Toggles the plotting backend

        Parameters
        ----------
        dmp_folders : tuple
            Tuple of the dmp_folder (output from bout_runners).
        collectPaths : [None|tuple]
            Tuple containing the collect paths
            If None, the dmp_folders will be used
        useSubProcess : bool
            Whether each job will be made by a new sub process, if not,
            the jobs will be done in series.
        """
        #}}}

        # Set collect paths
        if collectPaths == None:
            collectPaths = dmp_folders

        # Set the member data
        self._collectPaths  = collectPaths
        self._useSubProcess = useSubProcess

        if self._useSubProcess:
            #{{{ The multiprocess currently only works with the Agg backend
            # Qt4Agg currently throws
            # [xcb] Unknown sequence number while processing queue
            # [xcb] Most likely this is a multi-threaded client and XInitThreads has not been called
            # [xcb] Aborting, sorry about that.
            # python: ../../src/xcb_io.c:274: poll_for_event: Assertion
            # `!xcb_xlib_threads_sequence_lost' failed.
            #}}}
            plt.switch_backend("Agg")
    #}}}
#}}}
