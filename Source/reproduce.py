"""Compile the solver, run the six cases and regenerate the plotted results."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

here=Path(__file__).resolve().parent
sys.path.insert(0,str(here))
from postprocess import CASES

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=here.parent/"reproduced",
                        help="New directory, separate from the distributed project")
    parser.add_argument("--with-cad",action="store_true",help="Also build and verify STEP files (requires CadQuery)")
    args=parser.parse_args()
    out=args.output.resolve()
    if out.exists() and any(out.iterdir()):
        parser.error("Choose a new or empty output directory.")
    compiler=next((shutil.which(name) for name in ("c++","g++","clang++","cl") if shutil.which(name)),None)
    if not compiler: parser.error("Install a C++17 compiler or use a Visual Studio developer command prompt.")
    out.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,AIR_NOZZLE_OUTPUT=str(out))
    def run(cmd,cwd=None): subprocess.run([str(x) for x in cmd],env=env,cwd=cwd,check=True)
    with tempfile.TemporaryDirectory() as tmp:
        exe=Path(tmp)/("air_nozzle.exe" if os.name=="nt" else "air_nozzle")
        cpp=here/"nozzle_euler.cpp"
        if Path(compiler).stem.lower()=="cl": run([compiler,"/nologo","/O2","/EHsc","/std:c++17",cpp,f"/Fe:{exe}"],tmp)
        else: run([compiler,"-O3","-std=c++17",cpp,"-o",exe],tmp)
        run([sys.executable,here/"write_mesh.py"])
        for case,nx,nr in CASES:
            run([exe,nx,nr,case,out/"Results"/f"{case}_{nx}x{nr}",80000])
    run([sys.executable,here/"postprocess.py"])
    run([sys.executable,here/"figures.py"])
    if args.with_cad:
        run([sys.executable,here/"build_cad.py"])
        run([sys.executable,here/"render_cad.py"])
        run([sys.executable,here/"report.py"])
    print(f"Reproduced outputs: {out}")
if __name__=="__main__": main()
