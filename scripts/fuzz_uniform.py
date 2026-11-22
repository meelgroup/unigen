#!/usr/bin/env python3
# Fuzz unigen's samples for validity and uniformity. Usage: ./fuzz_uniform.py [-j 16] [--seed N]

import argparse
import itertools
import multiprocessing
import os
import random
import shutil
import subprocess
import sys
from scipy.stats import binomtest, chi2

HERE = os.path.dirname(os.path.abspath(__file__))
SAT_DIR = os.path.normpath(os.path.join(HERE, "..", ".."))
RED, GREEN, YELLOW, NC = "\033[31m", "\033[32m", "\033[33m", "\033[0m"


def physical_cores():
    try:
        out = subprocess.run(["lscpu", "-p=core"], capture_output=True, text=True).stdout
        return len({l for l in out.splitlines() if not l.startswith("#")}) or os.cpu_count()
    except OSError:
        return os.cpu_count()


def set_up_parser():
    p = argparse.ArgumentParser(description="Fuzz unigen for validity and uniformity of its samples")
    p.add_argument("--unigen", default=os.path.join(HERE, "..", "build-native", "unigen"))
    p.add_argument("--cms", default=os.path.join(SAT_DIR, "cryptominisat", "build", "cryptominisat5"),
                   help="Reference solver: enumerates the projected solutions, checks samples")
    p.add_argument("--cfuzz", default=os.path.join(SAT_DIR, "count_fuzzer"),
                   help="count_fuzzer dir, for its CNF generators")
    p.add_argument("--seed", type=int, help="Run only this one test")
    p.add_argument("--start", type=int, help="First seed. Default: random")
    p.add_argument("-n", "--num", type=int, default=0, help="Number of tests, 0 = forever")
    p.add_argument("-j", "--jobs", type=int, default=physical_cores(), help="Default: physical cores")
    p.add_argument("--outdir", default="fuzz_out")
    p.add_argument("--maxsols", type=int, default=1000,
                   help="Enumerate the projected solutions only up to this many. Above it, no chi-square")
    p.add_argument("--persol", type=int, default=15, help="Samples per projected solution for chi-square")
    p.add_argument("--maxsamples", type=int, default=5000)
    p.add_argument("--kernel", type=float, default=0.5, help="Probability of a Barbarik kernel test")
    p.add_argument("--pval", type=float, default=1e-7, help="Reject uniformity below this p-value")
    p.add_argument("--tout", type=int, default=300, help="Timeout of a unigen run")
    p.add_argument("--extra", default="", help="Extra unigen options, given to every run, instead of their random values")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


class Fail(Exception):
    pass


class Skip(Exception):
    pass


def read_cnf(fname):
    nvars, cls = 0, []
    with open(fname) as f:
        for line in f:
            t = line.split()
            if not t or t[0] in ("c", "w"):
                continue
            if t[0] == "p":
                nvars = int(t[2])
                continue
            lits = [int(x) for x in t]
            assert lits[-1] == 0
            cls.append(lits[:-1])
    return nvars, cls


def write_cnf(fname, nvars, cls, show):
    with open(fname, "w") as f:
        f.write(f"p cnf {nvars} {len(cls)}\n")
        if show is not None:
            f.write("c p show " + "".join(f"{v} " for v in show) + "0\n")
        for c in cls:
            f.write(" ".join(map(str, c)) + " 0\n")


def gen_cnf(t, fname):
    seed = t.rnd.randint(0, 10**9)
    size = t.rnd.choice(["small", "small", "normal"])
    if t.rnd.random() < 0.5:
        call = [os.path.join(t.opts.cfuzz, "biere-fuzz"), str(seed), "0"]
        if size == "small":
            call += ["5", "12", "1"]
    else:
        call = [os.path.join(t.opts.cfuzz, "cnf-fuzz-brummayer.py"), "-s", str(seed), "-T", "0"]
        if size == "small":
            call += ["-i", "3", "-I", "12"]
    t.log(" ".join(call))
    with open(fname, "w") as f:
        subprocess.run(call, stdout=f, check=True, cwd=t.opts.cfuzz)
    nvars, cls = read_cnf(fname)
    show = None
    if nvars > 0 and t.rnd.random() < 0.8:
        k = t.rnd.choice([t.rnd.randint(0, max(1, nvars // 5)), t.rnd.randint(nvars // 4, nvars // 2 + 1),
                          t.rnd.randint(1, min(nvars, 14))])
        show = sorted(t.rnd.sample(range(1, nvars + 1), min(k, nvars)))
    return nvars, cls, show


def run(cmd, tout, cwd=None):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=tout, cwd=cwd)
    except subprocess.TimeoutExpired:
        raise Skip(f"timeout: {' '.join(cmd)}")
    return r


def cms_enum(t, fname, maxsol):
    """Projected solutions, as tuples of signed lits over the show vars, or None if more than maxsol"""
    r = run([t.opts.cms, "--verb", "0", "--maxsol", str(maxsol + 1), fname], t.opts.tout)
    sols, cur = [], []
    for line in r.stdout.splitlines():
        if line.startswith("v "):
            cur += [int(x) for x in line.split()[1:]]
            if cur and cur[-1] == 0:
                sols.append(tuple(sorted(cur[:-1], key=abs)))
                cur = []
    if r.returncode not in (10, 20):
        raise Fail(f"cms exited with {r.returncode} on {fname}")
    if len(sols) != len(set(sols)):
        raise Fail("cms gave a solution twice")
    return None if len(sols) > maxsol else set(sols)


def brute_enum(nvars, cls, show):
    sols = set()
    for bits in itertools.product([False, True], repeat=nvars):
        if all(any(bits[abs(l) - 1] == (l > 0) for l in c) for c in cls):
            sols.add(tuple(v if bits[v - 1] else -v for v in show))
    return sols


def cms_sat(t, fname, units):
    nvars, cls = read_cnf(fname)
    tmp = fname + ".check.cnf"
    write_cnf(tmp, nvars, cls + [[u] for u in units], None)
    r = run([t.opts.cms, "--verb", "0", tmp], t.opts.tout)
    os.unlink(tmp)
    if r.returncode not in (10, 20):
        raise Fail(f"cms exited with {r.returncode}")
    return r.returncode == 10


def chain_cls(c, lits):
    """Clauses with exactly c solutions over lits, LSB first: the number they encode is < c"""
    assert 1 <= c < 2 ** len(lits)
    ones = [j for j in range(len(lits)) if (c >> j) & 1]
    cls = [[-lits[j] for j in ones]]
    for i in range(len(lits)):
        if not (c >> i) & 1:
            cls.append([-lits[i]] + [-lits[j] for j in ones if j > i])
    return cls


def unigen_opts(t):
    r = t.rnd
    o = ["--arjun", str(int(r.random() < 0.8))]
    o += ["--withe", str(int(r.random() < 0.8))]
    o += ["--backbone", str(r.randint(0, 1))]
    o += ["--sbva", r.choice(["0", "1000"])]
    o += ["--autarky", str(r.randint(0, 1))]
    o += ["--eappmc", str(r.randint(0, 1))]
    o += ["--evivif", str(r.randint(0, 1))]
    o += ["--eiter2", str(r.randint(0, 1))]
    o += ["--sparse", str(r.randint(0, 1))]
    o += ["--reusemodels", str(r.randint(0, 1))]
    # without simplification the banning clauses pile up, and sampling slows down a lot
    o += ["--simplify", r.choice(["0", "1", "1", "2", "2"])]
    o += ["--fullsol", str(r.randint(0, 1))]
    if r.random() < 0.2:
        o += ["--multisample", "1"]
    if r.random() < 0.2:
        o += ["--kappa", str(r.choice([0.3, 0.638, 0.9]))]
    extra = t.opts.extra.split()
    ret = []
    for name, val in zip(o[::2], o[1::2]):
        if name not in extra:
            ret += [name, val]
    return ret + extra


def max_samples(t, opts):
    return min(1000, t.opts.maxsamples) if opt_val(opts, "--simplify") == "0" else t.opts.maxsamples


def opt_val(opts, name):
    return opts[opts.index(name) + 1] if name in opts else None


def run_unigen(t, fname, opts, nsamples):
    sfile = fname + ".samples"
    cmd = [t.opts.unigen, "--verb", "0", "--seed", str(t.rnd.randint(0, 10**6)), "--samples", str(nsamples),
           "--sampleout", sfile] + opts + [fname]
    t.cmds.append(" ".join(cmd))
    t.log(t.cmds[-1])
    r = run(cmd, t.opts.tout)
    if r.returncode < 0:
        raise Fail(f"unigen died with signal {-r.returncode}\n{r.stdout[-3000:]}{r.stderr[-3000:]}")
    samples = []
    if os.path.exists(sfile):
        with open(sfile) as f:
            for line in f:
                lits = [int(x) for x in line.split()]
                if not lits or lits[-1] != 0:
                    raise Fail(f"malformed sample line: {line!r}")
                samples.append(lits[:-1])
    return r, samples


def check_samples(nvars, cls, show, opts, samples, nsamples):
    """Returns the samples projected on show, sorted by var"""
    if len(samples) != nsamples:
        raise Fail(f"asked for {nsamples} samples, got {len(samples)}")
    full = opt_val(opts, "--fullsol") == "1"
    want_vars = list(range(1, nvars + 1)) if full else show
    proj = []
    for s in samples:
        if sorted(abs(l) for l in s) != want_vars:
            raise Fail(f"sample is not over the expected vars {want_vars}: {s}")
        if full:
            val = {abs(l): l > 0 for l in s}
            for c in cls:
                if not any(val[abs(l)] == (l > 0) for l in c):
                    raise Fail(f"full sample {s} falsifies clause {c}")
        ss = set(show)
        proj.append(tuple(sorted((l for l in s if abs(l) in ss), key=abs)))
    return proj


def chi2_check(t, what, proj, sols, opts):
    if opt_val(opts, "--multisample") == "1":
        return ""
    cnt = {s: 0 for s in sols}
    for s in proj:
        cnt[s] += 1
    if len(sols) == 1:
        return ""
    e = len(proj) / len(sols)
    stat = sum((o - e) ** 2 / e for o in cnt.values())
    p = chi2.sf(stat, len(sols) - 1)
    lo, hi = min(cnt.values()), max(cnt.values())
    msg = f"{what}: sols {len(sols)} samples {len(proj)} min/max/exp {lo}/{hi}/{e:.1f} chi2-p {p:.3g}"
    t.log(msg)
    if p < t.opts.pval:
        raise Fail("NOT UNIFORM " + msg)
    return f" {what}-p {p:.2g}"


def kernel_test(t, fname, nvars, cls, show, opts, sigma1, sols):
    """Barbarik's idea: a kernel with exactly two projected solutions of the input, sigma1 from the
    sampler and sigma2 another one, each lifted to c solutions over fresh vars by different chain
    formulas. A uniform sampler picks each half the time"""
    others = [s for s in sols if s != sigma1] if sols is not None else []
    if sols is None:
        for _ in range(5):
            r = run([t.opts.cms, "--verb", "0", "--polar", "rnd", "--seed", str(t.rnd.randint(0, 10**6)), fname],
                    t.opts.tout)
            lits = [int(x) for line in r.stdout.splitlines() if line.startswith("v ") for x in line.split()[1:]]
            ss = set(show)
            s2 = tuple(sorted((l for l in lits if abs(l) in ss), key=abs))
            if s2 != sigma1:
                others = [s2]
                break
    if not others:
        return ""
    sigma2 = t.rnd.choice(others)
    s1, s2 = dict((abs(l), l) for l in sigma1), dict((abs(l), l) for l in sigma2)
    d = t.rnd.choice([v for v in show if s1[v] != s2[v]])
    sel = s1[d]  # true: sigma1, false: sigma2
    m = t.rnd.randint(3, 9)
    c = t.rnd.randint(2 ** (m - 1), 2 ** m - 1)
    z = list(range(nvars + 1, nvars + m + 1))
    flip = [t.rnd.choice([1, -1]) for _ in z]
    perm = z[:]
    t.rnd.shuffle(perm)
    za = z
    zb = [f * v for f, v in zip(flip, perm)]
    kcls = [c2[:] for c2 in cls]
    for v in show:
        if v != d:
            kcls.append([-sel, s1[v]])
            kcls.append([sel, s2[v]])
    kcls += [[-sel] + cl for cl in chain_cls(c, za)]
    kcls += [[sel] + cl for cl in chain_cls(c, zb)]
    kshow = sorted(show + z)
    kname = fname + ".kernel.cnf"
    write_cnf(kname, nvars + m, kcls, kshow)

    def zval(sample, lits):
        val = {abs(l): l > 0 for l in sample}
        return sum(1 << i for i, l in enumerate(lits) if val[abs(l)] == (l > 0))

    ksols = set()
    for a, b in ((sigma1, za), (sigma2, zb)):
        for n in range(c):
            asg = {abs(l): l for l in a}
            for i, l in enumerate(b):
                bit = (n >> i) & 1
                asg[abs(l)] = abs(l) if bool(bit) == (l > 0) else -abs(l)
            ksols.add(tuple(asg[v] for v in kshow))
    assert len(ksols) == 2 * c

    nks = min(max(300, t.opts.persol * 2 * c), max_samples(t, opts))
    r, samples = run_unigen(t, kname, opts, nks)
    if r.returncode != 0:
        raise Fail(f"unigen exited with {r.returncode} on the kernel\n{r.stdout[-3000:]}")
    kproj = check_samples(nvars + m, kcls, kshow, opts, samples, nks)
    n1 = 0
    for s in kproj:
        if s not in ksols:
            raise Fail(f"kernel sample {s} is not a solution of the kernel")
        sp = tuple(l for l in s if abs(l) <= nvars)
        n1 += sp == sigma1
        if zval(s, za if sp == sigma1 else zb) >= c:
            raise Fail(f"kernel sample {s} breaks its chain formula")
    msg = chi2_check(t, "kernel", kproj, ksols, opts)
    if opt_val(opts, "--multisample") != "1":
        p = binomtest(n1, len(kproj), 0.5).pvalue
        t.log(f"kernel: sigma1 {n1}/{len(kproj)} binom-p {p:.3g}")
        if p < t.opts.pval:
            raise Fail(f"NOT UNIFORM kernel: sigma1 picked {n1} of {len(kproj)} times, binom-p {p:.3g}")
    return msg


class Test:
    def __init__(self, opts, seed):
        self.opts, self.seed = opts, seed
        self.rnd = random.Random(seed)
        self.cmds, self.lines = [], []
        self.dir = None

    def log(self, s):
        self.lines.append(s)


def unique_dir(outdir, seed, max_tries=10000):
    """Claims a fresh dir: mkdir fails if it exists, so parallel fuzzers never share one"""
    for i in range(1, max_tries + 1):
        d = os.path.join(outdir, f"seed-{seed}" + ("" if i == 1 else f"-{i}"))
        try:
            os.mkdir(d)
            return d
        except FileExistsError:
            pass
    print(f"ERROR: cannot create a unique dir for seed {seed} in {outdir}")
    sys.exit(1)


def one_test(t):
    d = t.dir = unique_dir(t.opts.outdir, t.seed)
    fname = os.path.join(d, "fuzz.cnf")
    nvars, cls, show = gen_cnf(t, fname)
    write_cnf(fname, nvars, cls, show)
    if show is None:
        show = list(range(1, nvars + 1))
    sols = cms_enum(t, fname, t.opts.maxsols)
    if sols is not None and nvars <= 12 and brute_enum(nvars, cls, show) != sols:
        raise Fail("cms' projected solutions differ from brute force")
    opts = unigen_opts(t)
    if sols is not None and not sols:
        r, _ = run_unigen(t, fname, opts, 10)
        if r.returncode == 0 or "unsatisfiable" not in r.stdout:
            raise Fail(f"UNSAT input, but unigen exited with {r.returncode}\n{r.stdout[-3000:]}")
        return d, "unsat"
    if sols is None:
        nsamples = 200
    else:
        nsamples = min(max(300, t.opts.persol * len(sols)), max_samples(t, opts))
    r, samples = run_unigen(t, fname, opts, nsamples)
    if r.returncode != 0:
        if sols is None and "unsatisfiable" in r.stdout and not cms_sat(t, fname, []):
            return d, "unsat"
        raise Fail(f"unigen exited with {r.returncode}\n{r.stdout[-3000:]}")
    proj = check_samples(nvars, cls, show, opts, samples, nsamples)
    if sols is not None:
        for s in proj:
            if s not in sols:
                raise Fail(f"sample {s} is not a projected solution")
        msg = f"sols {len(sols)}" + chi2_check(t, "input", proj, sols, opts)
    else:
        msg = "sols >" + str(t.opts.maxsols)
        if opt_val(opts, "--fullsol") != "1":
            for s in t.rnd.sample(sorted(set(proj)), min(10, len(set(proj)))):
                if not cms_sat(t, fname, list(s)):
                    raise Fail(f"sample {s} is not a projected solution")
    if show and t.rnd.random() < t.opts.kernel:
        msg += kernel_test(t, fname, nvars, cls, show, opts, proj[0], sols)
    return d, msg


def worker(args):
    opts, seed = args
    t = Test(opts, seed)
    try:
        _, msg = one_test(t)
        status, out = "ok", f"{GREEN}OK{NC}   seed {seed}: {msg}"
    except Skip as e:
        status, out = "skip", f"{YELLOW}SKIP{NC} seed {seed}: {e}"
    except Fail as e:
        status = "fail"
        out = f"{RED}FAIL{NC} seed {seed}: {e}\n" + "\n".join("  " + c for c in t.cmds)
        out += f"\n  files: {t.dir}\n  repro: {sys.argv[0]} --seed {seed}"
    if opts.verbose:
        out += "\n" + "\n".join("  | " + s for s in t.lines)
    if status != "fail" and opts.seed is None:
        if t.dir: shutil.rmtree(t.dir, ignore_errors=True)
    return status, out


def main():
    opts = set_up_parser().parse_args()
    for exe in (opts.unigen, opts.cms):
        if not os.access(exe, os.X_OK):
            print(f"ERROR: not executable: {exe}")
            sys.exit(1)
    os.makedirs(opts.outdir, exist_ok=True)
    if opts.seed is not None:
        status, out = worker((opts, opts.seed))
        print(out)
        sys.exit(status == "fail")
    start = opts.start if opts.start is not None else random.randint(0, 10**9)
    seeds = itertools.count(start) if opts.num == 0 else range(start, start + opts.num)
    stats = {"ok": 0, "skip": 0, "fail": 0}
    with multiprocessing.Pool(opts.jobs) as pool:
        for status, out in pool.imap_unordered(worker, ((opts, s) for s in seeds)):
            stats[status] += 1
            print(out, flush=True)
    print(f"done: {stats}")
    sys.exit(stats["fail"] > 0)


if __name__ == "__main__":
    main()
