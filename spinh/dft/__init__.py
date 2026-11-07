"""Set up, submit and harvest DFT jobs, one directory per spin configuration.

A job directory holds main.py, which builds a job set and runs the action named on the command
line: ``python main.py setup``; ``python main.py help`` lists the actions.

    CollinearVaspJobs, CollinearOrcaJobs    'u'/'d' configurations (isotropic exchange)
    SiteJobs                                one moment direction per site (VASP)
    SaxisJobs                               one spin axis for the whole system (VASP)
    recipes                                 ready-made directions and configurations
"""

from spinh.dft import recipes
from spinh.dft.collinear import CollinearVaspJobs
from spinh.dft.noncollinear import SaxisJobs, SiteJobs
from spinh.dft.orca import CollinearOrcaJobs

__all__ = ["CollinearOrcaJobs", "CollinearVaspJobs", "SaxisJobs", "SiteJobs", "recipes"]
