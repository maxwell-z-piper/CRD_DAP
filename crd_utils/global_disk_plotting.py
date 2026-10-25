"""Script 04 review figures. Descriptive curve distances are not error bars."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from .global_disk import disk_coordinates


COLORS=('tab:blue','tab:orange')


def _save(fig,path):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=180)
    fig.savefig(path.with_suffix('.pdf'))
    plt.close(fig)


def _map(problem,values):
    image=np.full(problem.light.shape,np.nan)
    all_values=np.full(len(problem.bin_ids),np.nan)
    all_values[problem.selected]=values
    image[problem.member_y,problem.member_x]=all_values[problem.member_bin]
    return image


def _draw_map(ax,problem,values,**kwargs):
    artist=ax.pcolormesh(problem.east,problem.north,_map(problem,values),shading='nearest',**kwargs)
    ax.set(xlabel='East offset (arcsec)',ylabel='North offset (arcsec)',aspect='equal')
    return artist


def _curves(ax,result,sign_a,*,reference=None):
    nodes=np.asarray(result['nodes'])
    k=len(nodes)
    theta=np.asarray(result['theta'])
    for d,sign in enumerate((sign_a,-sign_a)):
        if reference is not None:
            ref=np.asarray(reference['theta'])[5+d*k:5+(d+1)*k]
            ax.plot(np.r_[0,nodes],sign*np.r_[0,ref],color=COLORS[d],alpha=.25,lw=1)
        speed=theta[5+d*k:5+(d+1)*k]
        ax.plot(np.r_[0,nodes],sign*np.r_[0,speed],'-o',ms=3,color=COLORS[d],label='Disk '+('A','B')[d])
    ax.axvspan(0,nodes[0],color='gray',alpha=.08)
    ax.axhline(0,color='gray',lw=.5)
    ax.set(xlabel='Disk-plane radius (arcsec)',ylabel='Signed intrinsic speed (km/s)')


def calibration_summary(problem,results,path,metadata):
    """One wide figure, TWO rows by SIX entries, with two plots/text each.

    The zero-centered map is trial data cost minus the SAME control bin cost.
    Never color by the penalized objective, and never rescale colors per tile.
    """
    if len(results)!=12 or results[0]['coefficient']!=0:
        raise ValueError('The review sheet requires a control and eleven ordered strength trials')
    if any(r['bin_ids']!=results[0]['bin_ids'] or r['nodes']!=results[0]['nodes'] for r in results):
        raise ValueError('Calibration plots require identical bins and nodes')
    control=results[0]
    k=problem.nnode
    reference=np.asarray(control['theta'])
    changes=[np.asarray(r['bin_score'])-control['bin_score'] for r in results]
    limit=max(1e-9,max(float(np.max(np.abs(d))) for d in changes))
    norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit)
    fig=plt.figure(figsize=(40,16),layout=None)
    grid=fig.add_gridspec(3,6,height_ratios=(.65,3,3),hspace=.4,wspace=.25,
                         left=.025,right=.985,bottom=.085,top=.91)
    strip=fig.add_subplot(grid[0,:])
    previous=control
    distances=[]
    for result in results[1:]:
        a=np.asarray(result['theta'])[5:]
        b=np.asarray(previous['theta'])[5:]
        distances.append([np.sqrt(np.mean((a[d*k:(d+1)*k]-b[d*k:(d+1)*k])**2)) for d in range(2)])
        previous=result
    for d in range(2):
        strip.plot(np.arange(11),np.asarray(distances)[:,d],'-o',color=COLORS[d],label='Disk '+('A','B')[d])
    strip.set(xticks=np.arange(11),xticklabels=['/32','/16','/8','/4','/2','baseline','×2','×4','×8','×16','×32'],
              ylabel='Adjacent-trial RMS\nchange (km/s)',title='Stability across increasing strength; first point compares weakest trial with control. Small change is not a quality verdict.')
    strip.legend(loc='upper right',ncol=2)
    map_axes=[]
    for index,result in enumerate(results):
        row,col=divmod(index,6)
        cell=grid[row+1,col].subgridspec(2,2,width_ratios=(1.65,1),hspace=.45,wspace=.12)
        curve=fig.add_subplot(cell[0,0])
        _curves(curve,result,problem.sign_a,reference=control)
        curve.set_ylim(-1.08*problem.bounds[5,1],1.08*problem.bounds[5,1])
        curve.set_xlim(0,problem.nodes[-1]*1.03)
        curve.legend(fontsize=6,loc='best')
        label='CONTROL' if index==0 else ('BASELINE' if index==6 else f'trial {index}')
        curve.set_title(f'{label}: C={result["coefficient"]:.6g}',fontsize=10)
        image_ax=fig.add_subplot(cell[1,0])
        artist=_draw_map(image_ax,problem,changes[index],cmap='RdBu_r',norm=norm)
        image_ax.set_title('Data cost change from control',fontsize=9)
        map_axes.append(image_ax)
        text=fig.add_subplot(cell[:,1]);text.axis('off')
        theta=np.asarray(result['theta']);difference=theta-reference
        rms=[np.sqrt(np.mean((theta[5+d*k:5+(d+1)*k]-reference[5+d*k:5+(d+1)*k])**2)) for d in range(2)]
        lines=[f"lambda = {result['lambda']:.6g}",f"Qdata = {result['data_score']:.5g}",
               f"ΔQdata = {result['data_score']-control['data_score']:+.4g}",f"Penalty = {result['penalty']:.4g}",
               f"PA = {theta[2]:.3f}° ({difference[2]:+.3f})",f"i = {theta[3]:.3f}° ({difference[3]:+.3f})",
               f"Vsys = {theta[4]:.2f} ({difference[4]:+.2f}) km/s",
               f"x0 = {theta[0]:.3f} ({difference[0]:+.3f}) arcsec",
               f"y0 = {theta[1]:.3f} ({difference[1]:+.3f}) arcsec",
               f"RMS from control:",f"A {rms[0]:.2f}; B {rms[1]:.2f} km/s",'Flags:']
        lines += [flag.replace('_',' ').lower() for flag in result['flags']] or ['none']
        text.text(0,1,'\n'.join(lines),va='top',fontsize=7.5,wrap=True)
        curve.tick_params(labelsize=7);image_ax.tick_params(labelsize=7)
    bar=fig.add_axes([.30,.035,.4,.009])
    fig.colorbar(artist,cax=bar,orientation='horizontal',label='Per-bin Δ data chi-square relative to control: blue = improved, red = worsened')
    fig.suptitle('Script 04 — empirical smoothing calibration\n'+
        f"Full aperture; PSF enabled; {len(problem.selected)} bins; {k} free nodes/disk; "
        f"spacing {metadata['spacing_arcsec']:.3f} arcsec; Rmax {metadata['Rmax_arcsec']:.3f} arcsec\n"
        'Curvature excludes the imposed origin. Faint curves are the unpenalized control. No uncertainty bands are implied.',fontsize=14)
    _save(fig,path)


def solution_diagnostics(problem,result,path):
    fig,axes=plt.subplots(2,3,figsize=(17,10))
    _curves(axes[0,0],result,problem.sign_a)
    axes[0,0].legend()
    for ax,key,label,cmap in ((axes[0,1],'va','Disk A LOS velocity (km/s)','RdBu_r'),
                              (axes[0,2],'vb','Disk B LOS velocity (km/s)','RdBu_r'),
                              (axes[1,0],'fraction','Coarse profiled RH3 fraction','viridis'),
                              (axes[1,1],'bin_score','Data Δchi-square from independent minima','magma')):
        im=_draw_map(ax,problem,result[key],cmap=cmap)
        fig.colorbar(im,ax=ax,label=label)
    history=result['history']
    for key in ('data_score','objective'):
        axes[1,2].plot([h['seconds'] for h in history],[h[key] for h in history],label=key)
    axes[1,2].set(xlabel='Search time (s)',ylabel='Objective contribution')
    axes[1,2].legend()
    fig.suptitle(f"{result.get('name','Global model')} | Qdata={result['data_score']:.5g}, P={result['penalty']:.5g}\n"+
                 ', '.join(result['flags']),fontsize=10)
    fig.tight_layout()
    _save(fig,path)


def geometry_preview(problem,info,path):
    fig,axes=plt.subplots(1,2,figsize=(13,6))
    axes[0].pcolormesh(problem.east,problem.north,problem.light,shading='nearest',cmap='gray')
    x0,y0,pa,inc,_=problem.initial_geometry
    t=np.linspace(0,2*np.pi,200)
    p=np.deg2rad(pa)
    for r in problem.nodes:
        major=r*np.cos(t);minor=r*np.sin(t)*np.cos(np.deg2rad(inc))
        axes[0].plot(x0+major*np.sin(p)+minor*np.cos(p),y0+major*np.cos(p)-minor*np.sin(p),lw=.5)
    axes[0].scatter([x0],[y0],marker='+',color='red')
    axes[0].set(aspect='equal',xlabel='East offset (arcsec)',ylabel='North offset (arcsec)',
                title='Fixed reference geometry / node coverage')
    radius,_=disk_coordinates(
        problem.east-x0,problem.north-y0,pa,inc)
    members=np.isin(problem.member_bin,problem.selected)
    axes[1].hist(radius[problem.member_y[members],problem.member_x[members]],
                 bins=np.r_[0.,(problem.nodes[:-1]+problem.nodes[1:])/2,problem.nodes[-1]+info['spacing_arcsec']/2])
    axes[1].set(xlabel='Disk-plane radius (arcsec)',ylabel='Selected member-spaxel count',
                title=f"{info['extent']}: {len(problem.selected)} bins\nDisk A sign {problem.sign_a:+d}; {info['label']['source']}")
    fig.tight_layout();_save(fig,path)


def extent_comparison(results,path):
    """Compare corresponding settings over common radii, without score ranking."""
    fig,axes=plt.subplots(2,2,figsize=(13,10),squeeze=False)
    for row,psf in enumerate((False,True)):
        for col,regularized in enumerate((False,True)):
            ax=axes[row,col]
            for extent,style in (('two_sigma','--'),('full','-')):
                result=results[f'{extent}_psf{int(psf)}_smooth{int(regularized)}']
                k=len(result['nodes']);theta=np.asarray(result['theta'])
                for d,sign in enumerate((result['sign_a'],-result['sign_a'])):
                    ax.plot(result['nodes'],sign*theta[5+d*k:5+(d+1)*k],style,color=COLORS[d],
                            label=f"{extent} / Disk {('A','B')[d]}")
            ax.set(title=f'PSF {psf}; regularity {regularized}',xlabel='Disk-plane radius (arcsec)',ylabel='Signed intrinsic speed (km/s)')
            ax.legend(fontsize=8)
    fig.tight_layout();_save(fig,path)
