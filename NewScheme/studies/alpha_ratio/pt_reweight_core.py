"""Signed-yield pT reweighting with shared-input first-order covariance."""
from dataclasses import dataclass
import numpy as np

TAGS = ('btag','untag','vbftag')
CHANNELS = ('ee','mumu')
REGIONS = ('CR','SR','VR')


@dataclass
class Estimate:
    value: np.ndarray
    jac: np.ndarray

    def __add__(self, other):
        return Estimate(self.value+other.value, self.jac+other.jac)

    def __sub__(self, other):
        return Estimate(self.value-other.value, self.jac-other.jac)

    def __mul__(self, other):
        return Estimate(self.value*other.value,
                        self.jac*other.value[:,None]+other.jac*self.value[:,None])

    def __truediv__(self, other):
        if np.any(other.value<=0): raise ValueError('Nonpositive ratio denominator')
        return Estimate(self.value/other.value,
                        self.jac/other.value[:,None]-other.jac*(self.value/other.value**2)[:,None])

    def covariance(self, variance):
        return (self.jac*variance[None,:])@self.jac.T

    def total(self):
        return Estimate(self.value.sum(keepdims=True),self.jac.sum(axis=0,keepdims=True))


class Space:
    """Independent signed Z and DY cells, topology/tag/channel/region/mass/pT."""
    def __init__(self, arrays, cat):
        dy = np.stack([arrays[f'{cat}__{tag}__DY'] for tag in TAGS])
        data = np.stack([arrays[f'{cat}__{tag}__Data'] for tag in TAGS])
        tt = np.stack([arrays[f'{cat}__{tag}__TTplusWW'] for tag in TAGS])
        vz = np.stack([arrays[f'{cat}__{tag}__VZ'] for tag in TAGS])
        z = np.stack([data[...,0]-tt[...,0]-vz[...,0],
                      data[...,1]+tt[...,1]+vz[...,1]],axis=-1)
        self.mean = np.stack([z[...,0],dy[...,0]])
        self.var = np.stack([z[...,1],dy[...,1]])
        self.shape = self.mean.shape
        self.flat = self.mean.ravel()
        self.variance = self.var.ravel()
        self.index = np.arange(self.flat.size).reshape(self.shape)
        self.data = data[...,0]
        self.mass_edges = arrays['mass_edges__'+cat]
        self.pt_edges = arrays['pt_edges']
        assert self.shape[1:4] == (3,2,3)
        assert np.all(np.isfinite(self.flat)) and np.all(self.variance>=0)

    def indices(self, proc, tags, channels, regions, masses=None, pts=None):
        masses = np.arange(self.shape[-2]) if masses is None else masses
        pts = np.arange(self.shape[-1]) if pts is None else pts
        return self.index[np.ix_([proc],tags,channels,regions,masses,pts)].ravel()

    def linear(self, proc, tags, channels, regions, massgroups=None, ptgroups=None):
        if massgroups is None: massgroups=[np.arange(self.shape[-2])]
        if ptgroups is None: ptgroups=[np.arange(self.shape[-1])]
        rows=[]
        for mg in massgroups:
            for pg in ptgroups:
                row=np.zeros(self.flat.size)
                row[self.indices(proc,tags,channels,regions,mg,pg)]=1
                rows.append(row)
        jac=np.stack(rows)
        return Estimate(jac@self.flat,jac)

    def weighted(self,tags,channels,regions,massgroups,ptgroups,factors):
        values=[]; gradients=[]
        for mg in massgroups:
            dy=self.linear(1,tags,channels,regions,[mg],ptgroups)
            product=dy*factors
            total=product.total()
            values.append(total.value[0]);gradients.append(total.jac[0])
        return Estimate(np.array(values),np.stack(gradients))


def group_edges(edges, groups):
    return np.array([edges[group[0]] for group in groups]+[edges[groups[-1][-1]+1]])


def merge_bad_groups(groups, bad):
    i=int(np.flatnonzero(bad)[-1])
    j=i-1 if i==len(groups)-1 else i
    groups[j]=np.r_[groups[j],groups[j+1]]
    del groups[j+1]


def fit_pt(space, tags, regions, channels=(0,1), min_effective=25.):
    groups=[np.array([i]) for i in range(space.shape[-1])]
    while True:
        z=space.linear(0,tags,channels,regions,ptgroups=groups)
        dy=space.linear(1,tags,channels,regions,ptgroups=groups)
        vz=np.diag(z.covariance(space.variance))
        vd=np.diag(dy.covariance(space.variance))
        nz=np.divide(z.value**2,vz,out=np.zeros_like(vz),where=vz>0)
        nd=np.divide(dy.value**2,vd,out=np.zeros_like(vd),where=vd>0)
        bad=(z.value<=0)|(dy.value<=0)|(nz<min_effective)|(nd<min_effective)
        if not bad.any(): break
        if len(groups)==1:
            raise ValueError('Insufficient integrated calibration statistics')
        merge_bad_groups(groups,bad)
    return dict(groups=groups,edges=group_edges(space.pt_edges,groups),
                factor=z/dy,observation=z,dy=dy,effective_z=nz,effective_dy=nd)


def mass_groups(space, tag, fits):
    groups=[np.array([i]) for i in range(space.shape[-2])]
    while True:
        bad=np.zeros(len(groups),dtype=bool)
        for channel in (0,1):
            for region in (0,1,2):
                base=space.linear(1,[tag],[channel],[region],groups)
                bad |= base.value<=0
                for fit in fits:
                    weighted=space.weighted([tag],[channel],[region],groups,fit['groups'],fit['factor'])
                    bad |= weighted.value<=0
        if not bad.any(): return groups
        if len(groups)==1: raise ValueError('Nonpositive integrated DY prediction')
        merge_bad_groups(groups,bad)


def predictions(space,tag,channels,region,massgroups,fit,cr_regions=(0,)):
    if region in cr_regions:
        raise ValueError('Alpha source overlaps its target region')
    # Transfer each channel separately; common fitted pT factors share covariance.
    parts=[]
    for ch in channels:
        zc=space.linear(0,[tag],[ch],cr_regions,massgroups)
        zt=space.linear(0,[tag],[ch],[region],massgroups)
        dc=space.linear(1,[tag],[ch],cr_regions,massgroups)
        dt=space.linear(1,[tag],[ch],[region],massgroups)
        wc=space.weighted([tag],[ch],cr_regions,massgroups,fit['groups'],fit['factor'])
        wt=space.weighted([tag],[ch],[region],massgroups,fit['groups'],fit['factor'])
        norm=fit['observation'].total()/fit['dy'].total()
        norm=Estimate(np.repeat(norm.value,len(dt.value)),np.repeat(norm.jac,len(dt.value),axis=0))
        parts.append(dict(observed=zt,dy_nominal=dt,dy_reweighted=wt,
                          dy_normonly=dt*norm,
                          alpha_nominal=(dt/dc)*zc,alpha_reweighted=(wt/wc)*zc,
                          transfer_nominal=dt/dc,transfer_reweighted=wt/wc))
    result={key:parts[0][key] for key in parts[0]}
    for part in parts[1:]:
        for key in result: result[key]=result[key]+part[key]
    # Sum of transfer factors is not a combined alpha: do not expose it.
    if len(channels)>1:
        del result['transfer_nominal'];del result['transfer_reweighted']
    return result


def compare(observed,predicted,variance,valid=None):
    residual=observed-predicted
    cov=residual.covariance(variance)
    if valid is None: valid=np.ones(len(residual.value),dtype=bool)
    matrix=cov[np.ix_(valid,valid)]
    y=residual.value[valid]
    if matrix.size:
        eigen,rotation=np.linalg.eigh(matrix)
        keep=eigen>max(float(eigen.max())*1e-10,1e-12)
        projected=rotation[:,keep].T@y
        chi2=float(np.sum(projected**2/eigen[keep]))
        rank=int(keep.sum())
    else: chi2,rank=np.nan,0
    ot,pt=observed.total(),predicted.total()
    residual_total=residual.total()
    ratio=ot/pt if pt.value[0]>0 else None
    return dict(observed=float(ot.value[0]), observed_stat=float(np.sqrt(ot.covariance(variance)[0,0])),
                predicted=float(pt.value[0]), predicted_stat=float(np.sqrt(pt.covariance(variance)[0,0])),
                closure=float(ratio.value[0]) if ratio else np.nan,
                closure_stat=float(np.sqrt(max(0.,ratio.covariance(variance)[0,0]))) if ratio else np.nan,
                residual=float(residual_total.value[0]),
                residual_stat=float(np.sqrt(max(0.,residual_total.covariance(variance)[0,0]))),
                chi2=chi2,rank=rank, bins_used=int(valid.sum()))
