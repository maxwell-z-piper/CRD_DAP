"""Scientific invariants of the global model, independent of the test galaxy."""
import numpy as np
from crd_utils.global_disk import (DiskProblem,curvature_penalty,disk_coordinates,
                                   fit_disk,profile_cube,radial_nodes)


def test_interpolate_before_fraction_profile_does_not_mix_incompatible_fractions():
    cube=np.zeros((1,2,2,2))
    cube[0,:,:,0]=[[0,100],[100,0]]
    cube[0,:,:,1]=[[100,0],[0,100]]
    score,fraction=profile_cube(cube,np.array([0.,1.]),np.array([0.,1.]),[.2,.8],[.5],[.5])
    assert np.allclose(score,50.)
    assert fraction[0] in (.2,.8)
    # Profiling the corners first would incorrectly give zero.
    assert np.min(cube,axis=-1).mean()==0


def test_exact_node_ignores_failed_unused_corners_and_outside_is_rejected():
    cube=np.full((1,2,2,2),np.inf)
    cube[0,0,0]=[3.,4.]
    axes=np.array([0.,1.])
    assert profile_cube(cube,axes,axes,[.2,.8],[0.],[0.])[0][0]==3
    assert profile_cube(cube,axes,axes,[.2,.8],[.5],[.5]) is None
    assert profile_cube(cube,axes,axes,[.2,.8],[1.01],[0.]) is None


def test_component_exchange_preserves_score_and_complements_fraction():
    cube=np.arange(3*3*3,dtype=float).reshape(1,3,3,3)
    axes=np.array([-10.,0.,10.]);fractions=np.array([.2,.5,.8])
    first=profile_cube(cube,axes,axes,fractions,[-4.],[7.])
    exchanged=cube.transpose(0,2,1,3)[...,::-1]
    second=profile_cube(exchanged,axes,axes,fractions,[7.],[-4.])
    assert np.allclose(first[0],second[0])
    assert np.allclose(first[1],1-second[1])


def test_central_anchor_is_not_in_curvature_and_extent_scaling_keeps_cost():
    assert curvature_penalty([150.,160.,170.],[100.,110.,120.],.002)==(0.,0.)
    inner_A=(6/1.)**4/200.**2
    outer_A=(12/1.)**4/200.**2
    c=.002
    assert np.isclose((c/inner_A)/16,c/outer_A)
    assert curvature_penalty([100,200,100],[100,100,100],c)[0]==80.


def test_sky_pa_and_outer_grid_coverage():
    r,cosine=disk_coordinates(np.array([0.,0.,1.]),np.array([1.,-1.,0.]),0.,60.)
    assert np.allclose(cosine[:2],[1.,-1.])
    assert np.isclose(cosine[2],0.)
    assert radial_nodes(1.,1.,3.1,outward=True)[-1]==4
    assert radial_nodes(1.,1.,3.1,outward=False)[-1]==3


def problem_fixture():
    north=np.array([[-3.,-2.,-1.,1.,2.,3.]])
    east=np.zeros_like(north)
    geometry=np.array([0.,0.,0.,30.,0.])
    bounds=np.vstack((np.column_stack([geometry,geometry]),np.tile([0.,200.],(6,1))))
    axis=np.linspace(-100,100,21)
    problem=DiskProblem(east,north,np.ones_like(north),np.zeros(6,int),np.arange(6),
        np.arange(6),np.arange(6),np.arange(6),np.zeros((6,21,21,3)),axis,axis,np.array([.2,.5,.8]),
        np.array([1.,2.,3.]),bounds,geometry,100.,1)
    truth=np.r_[geometry,[80.,120.,140.],[60.,100.,120.]]
    va,vb,*_=problem.predict(truth)
    aa,bb,ff=np.meshgrid(axis,axis,problem.fractions,indexing='ij')
    problem.chi2=np.array([((aa-a)/5)**2+((bb-b)/5)**2+((ff-.5)/.1)**2 for a,b in zip(va,vb)])
    return problem,truth


def test_known_coherent_rotation_curves_are_recovered():
    problem,truth=problem_fixture()
    result=fit_disk(problem,0.,seed=4,population_factor=3,global_iterations=12,
                    local_starts=3,local_iterations=120,tolerance=1e-7)
    assert np.allclose(result['theta'][5:],truth[5:],atol=.05)
    assert result['data_score']<.02
    assert result['bin_ids']==list(range(6))


def test_psf_prediction_weights_photons_rather_than_averaging_velocities():
    problem,_=problem_fixture()
    problem.east=np.array([[0.,1.,2.]])
    problem.north=np.zeros((1,3))
    problem.light=np.array([[100.,10.,0.]])
    problem.member_y=np.array([0,0]);problem.member_x=np.array([0,1])
    problem.member_bin=np.array([0,1]);problem.bin_ids=np.array([0,1]);problem.selected=np.array([0,1])
    problem.kernel=np.ones((1,3))/3
    theta=np.r_[[0,0,90,30,0],[200,200,200],[200,200,200]]
    va,vb,*_=problem.predict(theta)
    assert np.isclose(va[0],1000/110)
    assert np.isclose(vb[0],-1000/110)


def test_last_node_continues_without_changing_node_count():
    problem,truth=problem_fixture()
    problem.north*=2
    before=problem.nnode
    va,*_=problem.predict(truth)
    assert np.isclose(va[-1],70.)
    assert problem.nnode==before
