import os

cur = 32768
while cur <= 65536:
    source = f"{cur}_blob"
    print(source, flush=True)
    out = os.popen(f"powershell ../build/Release/Sim.exe {source}").read()
    out = out.split("\n")
    buh = out[-11:-3]
    for i in range(8):
        print(buh[i], flush=True)
    print("", flush=True)
    cur *= 2

cur = 32768
while cur <= 65536:
    source = f"{cur}_stream"
    print(source, flush=True)
    out = os.popen(f"powershell ../build/Release/Sim.exe {source}").read()
    out = out.split("\n")
    buh = out[-11:-3]
    for i in range(8):
        print(buh[i], flush=True)
    print("", flush=True)
    cur *= 2