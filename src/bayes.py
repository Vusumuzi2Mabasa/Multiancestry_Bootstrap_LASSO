"""Secondary comparators reused from Mabasa source study commit 03a5a078.
Published MultiSuSiE is an external GPL-3 dependency, not vendored here.
"""
from pathlib import Path
import os,sys,itertools,functools,math,contextlib,io
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['NUMBA_NUM_THREADS']='1'
from common import ROOT
import numpy as np
from scipy.special import ndtr,logsumexp
from scipy.optimize import brentq
from MultiSuSiE.susiepy_ss import susie_multi_ss
MASKS=np.array(list(itertools.product([0,1],repeat=8)),dtype=bool)
@functools.lru_cache(None)
def layout(rho,tau):
    groups=[]
    for s in range(1,9):
        ids=np.flatnonzero(MASKS.sum(1)==s);ix=np.array([np.flatnonzero(MASKS[i]) for i in ids])
        v=tau*tau*np.block([[np.eye(s),rho*np.eye(s)],[rho*np.eye(s),np.eye(s)]])
        groups.append((ids,ix,np.linalg.inv(v),np.linalg.slogdet(v)[1]))
    return groups
def mixture_quantiles(weights,means,variances,probs=(.025,.975)):
    """Equal-tailed quantiles including a genuine point mass at zero."""
    weights=np.asarray(weights);means=np.asarray(means);variances=np.maximum(variances,0)
    active=variances>1e-15;atom=float(weights[~active].sum())
    assert np.all(np.abs(means[~active])<1e-8)
    w=weights[active];m=means[active];sd=np.sqrt(variances[active])
    left=float(np.sum(w*ndtr(-m/sd))) if len(w) else 0.;right=left+atom
    bound=max(1.,float(np.max(np.abs(m)+12*sd)) if len(w) else 1.)
    def cdf(x):return float(np.sum(w*ndtr((x-m)/sd)))+atom*(x>=0)
    out=[]
    for p in probs:
        if left<=p<=right:out.append(0.)
        elif p<left:out.append(brentq(lambda x:float(np.sum(w*ndtr((x-m)/sd)))-p,-bound,0.,xtol=1e-10))
        else:out.append(brentq(lambda x:cdf(x)-p,0.,bound,xtol=1e-10))
    return out
def custom(Q,q,variance,rho=.8,tau=.15,pi=.2,intervals=False):
    means=np.zeros((256,16));vars=np.zeros_like(means)
    lw=MASKS.sum(1)*math.log(pi)+(8-MASKS.sum(1))*math.log1p(-pi)
    for ids,ix,pinv,ldv in layout(rho,tau):
        s=ix.shape[1];h=np.broadcast_to(pinv,(len(ids),2*s,2*s)).copy()
        h[:,:s,:s]+=Q[0][ix[:,:,None],ix[:,None,:]]/variance[0]
        h[:,s:,s:]+=Q[1][ix[:,:,None],ix[:,None,:]]/variance[1]
        u=np.concatenate([q[0][ix]/variance[0],q[1][ix]/variance[1]],axis=1)
        cov=np.linalg.inv(h);mu=np.einsum('bij,bj->bi',cov,u)
        lw[ids]+=-.5*ldv-.5*np.linalg.slogdet(h)[1]+.5*np.sum(u*mu,axis=1)
        target=np.concatenate([ix,ix+8],axis=1)
        means[ids[:,None],target]=mu;vars[ids[:,None],target]=np.diagonal(cov,axis1=1,axis2=2)
    weights=np.exp(lw-logsumexp(lw));pip=weights@MASKS;mean=weights@means
    sd=np.sqrt(np.maximum(0,weights@(vars+means**2)-mean**2));lo=hi=None
    if intervals:
        limits=np.array([mixture_quantiles(weights,means[:,j],vars[:,j]) for j in range(16)])
        lo=limits[:,0].reshape(2,8);hi=limits[:,1].reshape(2,8)
    return dict(pip=pip,mean=mean.reshape(2,8),sd=sd.reshape(2,8),lower=lo,upper=hi,converged=True)
def multisusie(Q,q,stats,variance,L=5,rho=.75,intervals=False):
    n=[s['n'] for s in stats];actual=sum((s['n']-1)*s['varY'] for s in stats)/(sum(n)-1)
    supplied=sum((s['n']-1)*v for s,v in zip(stats,variance))/(sum(n)-1)
    rs=[x/np.sqrt(np.diag(x)[:,None]*np.diag(x)[None,:]) for x in Q]
    with contextlib.redirect_stdout(io.StringIO()):
        f=susie_multi_ss(XTX_list=[x.copy() for x in Q],XTY_list=[x.copy() for x in q],
            YTY_list=[(s['n']-1)*v for s,v in zip(stats,variance)],population_sizes=n,
            rho=np.array([[1.,rho],[rho,1.]]),L=L,standardize=False,
            scaled_prior_variance=.2*actual/supplied,estimate_residual_variance=False,
            estimate_prior_variance=True,pop_spec_effect_priors=True,estimate_prior_method='early_EM',
            R_list=rs,float_type=np.float64,max_iter=200,tol=.001,min_abs_corr=.5)
    lo=hi=None
    if intervals:
        active=np.any(f.V>1e-9,axis=0);alpha=f.alpha.copy();alpha[~active]=0
        configs=np.array(list(itertools.product([0,1],repeat=L)),dtype=bool)
        lo=np.empty((2,8));hi=np.empty((2,8))
        for k in range(2):
            for j in range(8):
                a=alpha[:,j];w=np.prod(np.where(configs,a[None,:],1-a[None,:]),axis=1)
                m=configs@f.mu[k,:,j];v=configs@np.maximum(0,f.mu2[k,k,:,j]-f.mu[k,:,j]**2)
                lo[k,j],hi[k,j]=mixture_quantiles(w,m,v)
                assert np.isclose(w@m,f.coef[k,j],atol=1e-7)
    return dict(pip=f.pip,mean=f.coef,sd=f.coef_sd,lower=lo,upper=hi,converged=bool(f.converged),iterations=int(f.niter))
def match_recovery(selected,causal,R,threshold=.8):
    """Maximum cardinality matching; one selected proxy cannot recover two causes."""
    chosen=np.flatnonzero(selected);truth=np.flatnonzero(causal);best=0
    def rec(i,used,count):
        nonlocal best
        if i==len(truth):best=max(best,count);return
        rec(i+1,used,count)
        for j in chosen:
            if j not in used and (j==truth[i] or R[j,truth[i]]**2>=threshold):rec(i+1,used|{j},count+1)
    rec(0,set(),0)
    return best/max(1,len(truth))
