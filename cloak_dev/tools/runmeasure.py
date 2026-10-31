"""Run a command with a wall-clock timeout and an address-space cap; report wall time and child peak RSS.
Usage: runmeasure.py <timeout_s> <mem_cap_gb> -- cmd args...
Prints one TSV line: status  wall_s  peak_rss_mb"""
import resource, subprocess, sys, time

timeout = float(sys.argv[1])
cap = int(float(sys.argv[2]) * 1024**3)
cmd = sys.argv[sys.argv.index('--') + 1:]

def limit():
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

t0 = time.time()
status = 'OK'
try:
    p = subprocess.run(cmd, preexec_fn=limit, timeout=timeout,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if p.returncode != 0:
        status = f'EXIT{p.returncode}'
        sys.stderr.write(p.stderr.decode(errors='replace')[-400:])
except subprocess.TimeoutExpired:
    status = 'TIMEOUT'
wall = time.time() - t0
rss_mb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
print(f'{status}\t{wall:.2f}\t{rss_mb:.1f}')
