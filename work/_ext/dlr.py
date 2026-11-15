import urllib.request, urllib.parse, os, sys, time
base="https://deb.debian.org/debian/pool/main/g/gmsh/"
budget=float(sys.argv[1]) if len(sys.argv)>1 else 100.0
files=["gmsh_4.15.2+ds1-1_loong64.deb","libgmsh-dev_4.15.2+ds1-1_loong64.deb","gmsh_4.15.2+ds1.orig.tar.xz"]
t0=time.time()
for f in files:
    url=base+urllib.parse.quote(f)
    size=int(urllib.request.urlopen(urllib.request.Request(url,method="HEAD"),timeout=20).headers["Content-Length"])
    have=os.path.getsize(f) if os.path.exists(f) else 0
    if have>=size: print("complete",f,size//1024,"KiB"); continue
    if time.time()-t0>budget: print("partial",f,have//1024,"/",size//1024,"KiB"); continue
    req=urllib.request.Request(url, headers={"Range":f"bytes={have}-"})
    with urllib.request.urlopen(req,timeout=60) as r, open(f,"ab") as o:
        n=have
        while time.time()-t0<=budget:
            b=r.read(1<<16)
            if not b: break
            o.write(b); n+=len(b)
    print(("complete" if n>=size else "partial"),f,n//1024,"/",size//1024,"KiB", flush=True)
