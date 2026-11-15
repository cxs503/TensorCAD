import urllib.request, urllib.parse, os, tarfile, io, lzma, sys
base="https://deb.debian.org/debian/pool/main/g/gmsh/"
def fetch(f):
    if not os.path.exists(f):
        d=urllib.request.urlopen(base+urllib.parse.quote(f),timeout=90).read(); open(f,"wb").write(d)
    return open(f,"rb").read()
def ar_members(buf):
    assert buf[:8]==b"!<arch>\n"; i=8; out=[]
    while i+60<=len(buf):
        hdr=buf[i:i+60]; name=hdr[:16].decode().strip(); size=int(hdr[48:58].decode().strip())
        out.append((name,buf[i+60:i+60+size])); i+=60+size+(size%2)
    return out
for deb in sys.argv[1:]:
    buf=fetch(deb)
    for name,data in ar_members(buf):
        if name.startswith("data.tar"):
            tf=tarfile.open(fileobj=io.BytesIO(data))
            tf.extractall("gmsh_pkg")
            print(deb,"->",name,len(data)//1024,"KiB")
