"""Noncollinear VASP jobs with one MAGMOM direction per site.

Set up the configurations below, then run an action:

    python main.py print          # configurations (local frames) and their directories
    python main.py setup          # create the directories and write the inputs
    python main.py submit
    python main.py check          # convergence;  resubmit: resubmit the unconverged jobs
    python main.py restart-test   # restart plan; restart: do it (myjob.restart_options)
    python main.py energies       # energies.dat;  dipoles: dipoles.dat
    python main.py help           # all actions

do() lists directions of the magnetic moments; do_spin() lists directions of the spins (opposite
to the moments). Pick one at the bottom.
"""

import sys

import data
from spinh.dft.noncollinear import SiteJobs
from spinh.dft.vasp import VaspJob


def check_and_resubmit():

    dirnames = [ "4319101_0_0", "3238801_0_0" ]
    for dirname in dirnames:
        VaspJob(data, dirname).check_convergence(restart=False)
        # VaspJob(data, dirname).check_convergence(restart=True)

def find_global_directions():
    direction_emin = (42.9,  102.5)
    myjob = SiteJobs(data, n_site=1, magnetic_moments=[3], n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(sites=[0], z_directions=[direction_emin], from_file=False)
    myjob.add_configurations([[[2,90]]], local_ref_frame=True)
    myjob.print_dirs_and_configs(local_ref_frame=False)

def transform_reference_frame(emat):

    myjob = SiteJobs(data, n_site=1, magnetic_moments=[1], n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(sites=[0], from_file=True)
    myjob.add_configurations([[[  0, 0]]], local_ref_frame=False)
    direction_emin = (0.0, 0.0)
    myjob.setup_local_ref_frames(sites=[0], z_directions=[direction_emin], emats=[emat], from_z_direction=False, from_file=False)
    myjob.transform_configurations(global2local = True)
    myjob.print_dirs_and_configs(local_ref_frame=True)

def do(direction_emin = (0.0, 0.0)):

    myjob = SiteJobs(data, n_site=10, magnetic_moments=[7, 7, 4, 4, 4, 4, 4, 4, 4, 4], directions_of="magmom", n_theta=1801, n_phi=3600)
    # myjob.setup_local_ref_frames(sites=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], from_file=True)

    ### Sample magnetic configurations

    # uuuuuuuuuu
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0]]], local_ref_frame=False)
    # dduuuuuuuu
    myjob.add_configurations([[[180,0], [180,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0]]], local_ref_frame=False)
    # uduuuudddd
    myjob.add_configurations([[[  0,0], [180,0], [  0,0], [  0,0], [  0,0], [  0,0], [180,0], [180,0], [180,0], [180,0]]], local_ref_frame=False)
    # uuuudduudd
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [  0,0], [180,0], [180,0], [  0,0], [  0,0], [180,0], [180,0]]], local_ref_frame=False)
    # uuudduuddu
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [180,0], [180,0], [  0,0], [  0,0], [180,0], [180,0], [  0,0]]], local_ref_frame=False)
    # uuududdudu
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [180,0], [  0,0], [180,0], [180,0], [  0,0], [180,0], [  0,0]]], local_ref_frame=False)
    # uuudududud
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [180,0], [  0,0], [180,0], [  0,0], [180,0], [  0,0], [180,0]]], local_ref_frame=False)
    # uddddduuuu
    myjob.add_configurations([[[  0,0], [180,0], [180,0], [180,0], [180,0], [180,0], [  0,0], [  0,0], [  0,0], [  0,0]]], local_ref_frame=False)
    # uuuuuudddd
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [180,0], [180,0], [180,0], [180,0]]], local_ref_frame=False)
    # uuudududdd
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [180,0], [  0,0], [180,0], [  0,0], [180,0], [180,0], [180,0]]], local_ref_frame=False)
    # uuududdudd
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [180,0], [  0,0], [180,0], [180,0], [  0,0], [180,0], [180,0]]], local_ref_frame=False)
    # uuuuududud
    myjob.add_configurations([[[  0,0], [  0,0], [  0,0], [  0,0], [  0,0], [180,0], [  0,0], [180,0], [  0,0], [180,0]]], local_ref_frame=False)

    ### Fix or flip one spin for all configurations.

    # myjob.set_direction_for_one_site(i_site=0, direction=[180, 0], local_ref_frame=True)
    # myjob.set_direction_for_one_site(i_site=1, direction=[180, 0], local_ref_frame=True)
    # myjob.set_direction_for_one_site(i_site=2, direction=[180, 0], local_ref_frame=True)
    # myjob.flip_direction_for_one_site(i_site=2)

    ### If more than one local reference frames are used, it is necessary to update the local coordinates with respect to the last local reference frame.

    # myjob.transform_configurations(global2local = True)

    myjob.restart_options = dict(max_angle=45, from_neighbor=True, file_name="WAVECAR", copy_it=False)
    return myjob

def do_spin(direction_emin = (0.0, 0.0)):

    myjob = SiteJobs(data, n_site=3, magnetic_moments=[2, 2, 3], directions_of="spin", n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(sites=[0, 1, 2], from_file=True)

    ### Get initial CHGCAR or WAVECAR.

    # myjob.add_configurations([[[  0,0], [  0,0], [  0,0]]], local_ref_frame=False)
    # myjob.add_configurations([[[180,0], [180,0], [180,0]]], local_ref_frame=False)

    ### Search for the local z-axis.

    # myjob.setup_local_ref_frames(sites=[0, 1, 2], z_directions=[direction_emin, direction_emin, direction_emin], from_file=False)

    # myjob.add_thetas_colinear_spin(phi=180, theta_min=6, theta_max=0, ntheta=4, local_ref_frame=True)
    # myjob.add_thetas_colinear_spin(  phi=0, theta_min=0, theta_max=6, ntheta=4, local_ref_frame=True)
    # myjob.add_thetas_colinear_spin(phi=270, theta_min=6, theta_max=0, ntheta=4, local_ref_frame=True)
    # myjob.add_thetas_colinear_spin( phi=90, theta_min=0, theta_max=6, ntheta=4, local_ref_frame=True)
    # myjob.add_thetas_colinear_spin(phi=180, theta_min=174, theta_max=180, ntheta=4, local_ref_frame=True)
    # myjob.add_thetas_colinear_spin(  phi=0, theta_min=180, theta_max=174, ntheta=4, local_ref_frame=True)
    # myjob.add_thetas_colinear_spin(phi=270, theta_min=174, theta_max=180, ntheta=4, local_ref_frame=True)
    # myjob.add_thetas_colinear_spin( phi=90, theta_min=180, theta_max=174, ntheta=4, local_ref_frame=True)

    ### Sample three great circles.

    myjob.add_thetas_colinear_spin(  phi=0, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    myjob.add_thetas_colinear_spin(phi=180, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    myjob.add_thetas_colinear_spin( phi=90, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    myjob.add_thetas_colinear_spin(phi=270, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    myjob.add_phis_colinear_spin( theta=90,   phi_min= 0,   phi_max=330,   nphi=12, local_ref_frame=True)

    ### Exchange interaction for two sites.

    # myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[90,90], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    # myjob.add_configurations([[[90,90], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    # myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    # myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[90,90], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    # myjob.flip_direction_for_one_site(i_site=0)                     
                                                                    
    # myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[90,90], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    # myjob.add_configurations([[[90,90], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    # myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    # myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[90,90], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)

    ### Exchange interaction for more than two sites. Conquer pair by pair.

    # myjob.add_configurations([[[  0, 0], [ 90, 0], [ 90, 180], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90, 0], [ 90, 270], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90, 0], [180,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [ 90, 180], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [ 90, 270], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [180,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [ 90, 180], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [ 90, 270], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [180,   0], [180, 0]]], local_ref_frame=True)

    # myjob.add_configurations([[[  0, 0], [ 90, 0], [ 90,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90, 0], [ 90,  90], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90, 0], [  0,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [ 90,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [ 90,  90], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [  0,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [ 90,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [ 90,  90], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [  0,   0], [180, 0]]], local_ref_frame=True)

    # myjob.flip_direction_for_one_site(i_site=0)
    # myjob.flip_direction_for_one_site(i_site=3)

    # myjob.add_configurations([[[  0, 0], [ 90, 0], [ 90, 180], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90, 0], [ 90, 270], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90, 0], [180,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [ 90, 180], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [ 90, 270], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [180,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [ 90, 180], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [ 90, 270], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [180,   0], [180, 0]]], local_ref_frame=True)

    # myjob.add_configurations([[[  0, 0], [ 90, 0], [ 90,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90, 0], [ 90,  90], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90, 0], [  0,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [ 90,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [ 90,  90], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [ 90,90], [  0,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [ 90,   0], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [ 90,  90], [180, 0]]], local_ref_frame=True)
    # myjob.add_configurations([[[  0, 0], [180, 0], [  0,   0], [180, 0]]], local_ref_frame=True)

    ### Fix or flip one spin for all configurations.

    # myjob.set_direction_for_one_site(i_site=0, direction=[180, 0], local_ref_frame=True)
    # myjob.set_direction_for_one_site(i_site=1, direction=[180, 0], local_ref_frame=True)
    # myjob.set_direction_for_one_site(i_site=2, direction=[180, 0], local_ref_frame=True)
    # myjob.flip_direction_for_one_site(i_site=2)

    ### If more than one local reference frames are used, it is necessary to update the local coordinates with respect to the last local reference frame.

    # myjob.transform_configurations(global2local = True)

    myjob.restart_options = dict(max_angle=45, from_neighbor=True, file_name="WAVECAR", copy_it=False)
    return myjob

if __name__ == "__main__":

    # check_and_resubmit()

    # find_global_directions()

    # transform_reference_frame(some_emat)

    myjob = do(direction_emin = (0.0, 0.0))
    # myjob = do_spin(direction_emin = (0.0, 0.0))

    myjob.run(sys.argv[1:])
