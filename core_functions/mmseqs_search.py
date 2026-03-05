import subprocess
from core_functions.paths_and_parameters import *

def run_wsl(cmd, check=True):
    out = subprocess.run(["wsl", "bash", "-lc", cmd], capture_output=True, text=True)
    if check and out.returncode != 0:
        raise RuntimeError(
            f"WSL command failed (code {out.returncode}).\n"
            f"CMD:\n{cmd}\n\nSTDOUT:\n{out.stdout}\n\nSTDERR:\n{out.stderr}"
        )
    return out
  

def mmseqs_run_search(threads=8, search_params="-s 7.5"):
    run_wsl(f'mkdir -p "{TMP_DIR}"')
    cmd = (
        f'mmseqs search "{QUERY_DB}" "{TARGET_DB}" "{RESULT_DB}" "{TMP_DIR}" '
        f'{search_params} --threads {threads}'
    )
    return run_wsl(cmd)
  
# Convert MMseqs alignment DB to TSV including taxonomy + coverage
def convertalis_with_taxonomy():
    cmd = f'''
    mmseqs convertalis "{QUERY_DB}" "{TARGET_DB}" "{RESULT_DB}" "{TSV_ANNOT_WSL}" \
      --format-output "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,taxid,taxname,taxlineage"
    '''
    return run_wsl(cmd)
