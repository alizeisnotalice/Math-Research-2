"""Independent frozen-mark quadrature, dyadic mesh comparison, and symbolic expansion."""
import json
import math
import platform
from pathlib import Path
import struct
import sys

import mpmath as mp
import scipy
from scipy.integrate import quad
import sympy as sp

mp.mp.dps=68

def f32(x):
    return struct.unpack('f',struct.pack('f',x))[0]

def primitives(d,R,phi):
    v=mp.mpf(1)/d
    c=mp.exp(-v/2)
    b=mp.mpf('.5') if d==1 else -mp.exp(-mp.mpf(1)/3)/2
    base=1+b*c*mp.cos(phi)
    skew=b*c*(5*v+v*v)*mp.sin(phi)/4
    tpol=v-v*v+mp.mpf(3)/8*(3*v*v-v**3)-(3*v*v-6*v**3+v**4)/32
    second=v+b*c*tpol*mp.cos(phi)
    return b,base,skew/R,second/R**2

def geometry(k,phi,d,panels):
    R=mp.pi*k+phi/2
    b,base,first,second=primitives(d,R,phi)
    peak=mp.sqrt(d/(2*mp.pi))
    def f(y):
        h=y/(2*R)
        v=mp.sqrt(1+h*h)
        u=v+h
        shift=y+y*y/(2*R*(v+1))
        return peak*mp.exp(-d*y*y/2)*(1+b*mp.cos(phi+shift))/(u**3*v)
    nodes=list(map(mp.mpf,(-22,-12,-8,-5,-2,0,2,5,8,12,22)))
    integral=mp.quad(f,nodes)
    def alternate(y):
        h=y/(2*R)
        v=mp.sqrt(1+h*h)
        eta=y*y/(2*R*(v+1))
        even=1+4*h*h
        odd=-h*(3+4*h*h)/v
        return peak*mp.exp(-d*y*y/2)*(even*(1+b*mp.cos(phi+eta)*mp.cos(y))-b*odd*mp.sin(phi+eta)*mp.sin(y))
    check=2*mp.quad(alternate,[0,2,5,8,12,22])
    assert abs(integral-check)<mp.mpf(10)**-61
    radius,angle,amp,norm=map(float,(R,phi,b,peak))
    def double_source(y):
        h=y/(2*radius)
        v=math.hypot(1,h)
        return norm*math.exp(-d*y*y/2)*(1+amp*math.cos(angle+y+y*y/(2*radius*(v+1))))/((v+h)**3*v)
    double=sum(quad(double_source,a,b,epsabs=2e-13,epsrel=2e-13)[0] for a,b in zip(nodes,nodes[1:]))
    double_error=abs(mp.mpf(double)-integral)
    assert double_error<mp.mpf('2e-12')
    step=mp.mpf(24)/panels
    # Same dyadic mesh as the GPU, but direct integrand and high-precision arithmetic.
    absMesh=step*mp.fsum(f(-12+(i+mp.mpf('.5'))*step) for i in range(panels))
    stableMesh=step*mp.fsum(f(-12+(i+mp.mpf('.5'))*step)
        -peak*mp.exp(-d*(-12+(i+mp.mpf('.5'))*step)**2/2)
        *(1+b*mp.cos(phi-12+(i+mp.mpf('.5'))*step)) for i in range(panels))
    residual=integral-base
    return {'total':integral,'leading':base,'residual':residual,'first':first,'second':second,
        'afterTwo':residual-first-second,'methodDifference':abs(integral-check),'quadpackAbsError':double_error,
        'stableMeshResidual':stableMesh,'absoluteMesh':absMesh,
        'stableMeshError':stableMesh-residual,'absoluteMeshError':absMesh-integral}

def symbolic():
    h,s,y,v=sp.symbols('h s y v',real=True)
    c=sp.symbols('c')
    eta=y*y*s/(2*(sp.sqrt(1+y*y*s*s/4)+1))
    weight=sp.exp(-3*sp.asinh(y*s/2))/sp.sqrt(1+y*y*s*s/4)
    assert sp.series(weight,s,0,3).removeO() == 1-sp.Rational(3,2)*s*y+s*s*y*y
    assert sp.series(eta,s,0,3).removeO()==s*y*y/4
    even=sp.simplify(((sp.sqrt(1+h*h)+h)**3+(sp.sqrt(1+h*h)-h)**3)/(2*sp.sqrt(1+h*h)))
    assert even==1+4*h*h
    correction=sp.expand(v-v*v+sp.Rational(3,8)*(3*v*v-v**3)-(3*v*v-6*v**3+v**4)/32)
    first=(5*v+v*v)/4
    return {'parityWeight':str(even),'firstVariancePolynomial':str(first),'secondVariancePolynomial':str(correction),
        'twoWellsFirstCoefficientInUnitsExpMinusHalf':str(sp.simplify((first.subs(v,1)-first.subs(v,sp.Rational(1,3)))/2)),
        'twoWellsSecondCosCoefficientInUnitsExpMinusHalf':str(sp.simplify((correction.subs(v,1)-correction.subs(v,sp.Rational(1,3)))/2))}

def encode(v):
    if isinstance(v,dict):return {k:encode(x) for k,x in v.items()}
    if isinstance(v,list):return [encode(x) for x in v]
    if isinstance(v,mp.mpf):return mp.nstr(v,42)
    return v

def main():
    source,out=map(Path,sys.argv[1:3])
    data=json.loads(source.read_text())
    results=[]
    for row in data['rows']:
        mode=row['mode']; phase_id=mode//6
        for label in ['small','large']:
            if label not in row:continue
            hit=row[label]
            k,phi,_,n=hit['point'];k,n=int(k),int(n)
            phi=mp.pi*phase_id/2 if phase_id<4 else mp.mpf(f32(phi))
            parts=[geometry(k,phi,d,n) for d in (1,3)]
            pair=sum(x['residual'] for x in parts)
            errors=[]
            for ix,part in enumerate(parts):
                errors.append({'stableGPUAbsoluteError':abs(mp.mpf(f32(hit['stable'][ix]))-part['residual']),
                    'rawGPUAbsoluteError':abs(mp.mpf(f32(hit['raw'][ix]))-part['total'])})
            results.append({'mode':mode,'kind':label,'k':k,'panels':n,'phi':phi,'parts':parts,'oneWellErrors':errors,
                'pairResidual':pair,'pairAfterTwo':sum(x['afterTwo'] for x in parts),
                'stableGPUPairAbsoluteError':abs(mp.mpf(f32(hit['stable'][2]))-pair),
                'rawGPUPairAbsoluteError':abs(mp.mpf(f32(hit['raw'][2]))-pair),
                'stableGPUPairResidual':hit['stable'][2],'rawGPUPairResidual':hit['raw'][2]})
    fixed=[]
    for phase_id,phi in [('axis',mp.mpf(0)),('intermediate',mp.mpf('1.15')),('quadrature',mp.pi/2)]:
        for k in (2,8,64,4096,65536):
            parts=[geometry(k,phi,d,128) for d in (1,3)]
            R=mp.pi*k+phi/2
            fixed.append({'phaseClass':phase_id,'phi':phi,'k':k,
                'pairResidual':sum(x['residual'] for x in parts),
                'scaledFirst':R*sum(x['residual'] for x in parts),
                'scaledSecond':R**2*sum(x['residual']-x['first'] for x in parts),
                'afterTwo':sum(x['afterTwo'] for x in parts),'R3AfterTwo':R**3*sum(x['afterTwo'] for x in parts),
                'predictedFirst':R*sum(x['first'] for x in parts),'predictedSecond':R**2*sum(x['second'] for x in parts)})
    out.write_text(json.dumps(encode({'python':platform.python_version(),'mpmath':mp.__version__,
        'scipy':scipy.__version__,'sympy':sp.__version__,'symbolic':symbolic(),'checked':results,'fixed':fixed}),
        indent=2)+'\n')
if __name__=='__main__':main()
