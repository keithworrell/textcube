import sys, os
H=os.path.dirname(os.path.abspath(__file__))
mode=sys.argv[1]
t=open(os.path.join(H,'template.html')).read(); w=open(os.path.join(H,'core.js')).read()
assert '</script' not in w
if mode=='local':
    t=t.replace('__THREE__','./node_modules/three/build/three.module.js').replace('__ADDONS__','./node_modules/three/examples/jsm/')
    out=os.path.join(H,'..','local.html')
else:
    t=t.replace('__THREE__','https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js').replace('__ADDONS__','https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/')
    out=os.path.join(H,'..','index.html')
open(out,'w').write(t.replace('__WORKER__',w)); print(out)
