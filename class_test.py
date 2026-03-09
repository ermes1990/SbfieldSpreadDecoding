from tqdm import tqdm
from math import floor
import csv


class TestRun:
    def __init__(self, q, r, k, deletions, insertions, verbose=False):
        self.q = q
        self.r = r
        self.k = k
        self.n = r*k
        self.deletions = deletions
        self.insertions = insertions
        self.verbose = verbose

        # Fields
        self.Fq = GF(q, 'z')
        self.Fqn = GF(self.q^self.n, 'y')
        self.Fqk, self.embedding = self.Fqn.subfield(self.k, map=True)
        self.alpha = self.Fqk.gen()

        # Ambient vector space
        self.V = VectorSpace(self.Fq, self.n)

    # --- helper: multiply subfield element into big field ---
    def mul_subspace(self, R, z):
        R_basis = [self.Fqn(x.list()) for x in R.basis()]
        R_z_basis = [self.embedding(z) * x for x in R_basis]
        return self.V.subspace(R_z_basis)

    def expand_subspace(self, R_exp, R):
        if self.verbose:
            print(f"Expand: dim before = {R_exp.dimension()}")
        # pick random z not in prime subfield
        z = self.Fqk.random_element()
        while z in self.Fqk.prime_subfield():
            z = self.Fqk.random_element()
        R_z = self.mul_subspace(R, z)
        R_exp = R_exp + R_z
        if self.verbose:
            print(f"Dimension after expansion = {R_exp.dimension()}")
        return R_exp

    def contract_subspace(self, R_cont, R):
        if self.verbose:
            print(f"Contract: dim before = {R_cont.dimension()}")
        z = self.Fqk.random_element()
        while z in self.Fqk.prime_subfield():
            z = self.Fqk.random_element()
        R_z = self.mul_subspace(R, z)
        R_cont = R_cont.intersection(R_z)
        if self.verbose:
            print(f"Dimension after contraction = {R_cont.dimension()}")
        return R_cont

    def corrupt_encode_random(self):
        Vk = self.V.subspace([self.embedding(self.alpha**i) for i in range(self.k)])
        Fqk_basis = [self.Fqn(x.list()) for x in Vk.basis()]
        a = self.Fqn.random_element()
        while a == 0:
            a = self.Fqn.random_element()
        VA = self.V.subspace([a*x for x in Fqk_basis])
        VB = self.V.subspace([self.Fqn.random_element() for _ in range(self.insertions)])
        VA_red = self.V.subspace([VA.random_element() for _ in range(self.k - self.deletions)])
        while VA_red.dimension() != self.k - self.deletions:
            VA_red = self.V.subspace([VA.random_element() for _ in range(self.k - self.deletions)])
        R = VB + VA_red
        return R, VA

    def decode(self, R):
        max_exp = min(floor((self.V.dimension() - self.Fqk.degree()) / R.dimension()) - 1, self.k - 1)
        R_exp = self.V.subspace([self.Fqn(x.list()) for x in R.basis()])
        if max_exp > 0:
            for _ in range(max_exp - 1): # diemsnion starts from r with 0 expansion, doing (e-1) expansions will result in dimension ~ er
                R_exp = self.expand_subspace(R_exp, R)
        res = R_exp
        lim = self.Fqk.degree()
        res_old = res
        while res.dimension() > 0 and lim >= 0:
            res_old = res
            res = self.contract_subspace(res, R_exp)
            if res.dimension() == self.k:
                return res
            lim -= 1
        return res

    def decoding_advance(self, R):
        R_exp = self.V.subspace([self.Fqn(x.list()) for x in R.basis()])
        max_exp = min(floor((self.V.dimension()) / R.dimension()) - 1, self.k - 1)
        if max_exp > 0:
            for _ in range(max_exp - 1): # diemsnion starts from r with 0 expansion, doing (e-1) expansions will result in dimension ~ er
                R_exp = self.expand_subspace(R_exp, R)
        res = R_exp
        res_old = res
        lim = self.Fqk.degree()
        while res.dimension() > 0 and lim >= 0:
            res_old = res
            res = self.contract_subspace(res, R_exp)
            lim -= 1
        res = res_old
        while res.dimension() < self.Fqk.degree() and res.dimension() > 0:
            res = self.expand_subspace(res, res)
        return self.expand_subspace(res, res)

    def extract(self, R, dim_extract= 5):
        C = self.V.subspace([])
        n = self.V.dimension()
        r = R.dimension()
        max_exp = min(floor((self.V.dimension() - self.Fqk.degree() - 4) /( R.dimension() - 2)) - 1, self.k - 1)
        #print(f"ext exp:{max_exp}")
        while C.dimension() < dim_extract:
            y = self.Fqk.random_element()
            z = self.Fqk.random_element()
            yR = self.mul_subspace(R, y)
            zR = self.mul_subspace(R, z)
            for _ in range(max_exp - 1):
                zR = self.expand_subspace(zR, R)
            Extract = yR.intersection(zR)
            if Extract.dimension() < R.dimension():
                C = C + Extract
        return C

    def run(self, trials=100, mask=7):
        print(f"Test: ({self.deletions},{self.insertions})")
        counter_alg1 = 0
        counter_alg2 = 0
        counter_alg3 = 0
        for _ in tqdm(range(trials)):
            R, VA = self.corrupt_encode_random()
            if mask & 1:
                for _ in range(5):
                    dec_A = self.decode(R)
                    if dec_A.dimension() == self.k:
                        if dec_A == VA:
                            counter_alg1 += 1
                        break
            if mask & 2:
                for _ in range(5):
                    dec_A = self.decoding_advance(R)
                    if dec_A.dimension() == self.k:
                        if dec_A == VA:
                            counter_alg2 += 1
                        break
            if mask & 4:
                C = self.extract(R, dim_extract=5)
                for _ in range(5):
                    dec_A = self.decoding_advance(C)
                    if dec_A.dimension() == self.k:
                        if dec_A == VA:
                            counter_alg3 += 1
                        break
        return counter_alg1, counter_alg2, counter_alg3


# Set seed to reproduce same outcome
set_random_seed(67)

# Run multiple test, save them in csv file ready to plot in latex
results = []

for ins in range(1,16):
    tester = TestRun(2, 6, 5, 3, ins, verbose=False)
    alg1, alg2, alg3 = tester.run(1000)
    results.append([ins, alg1, alg2, alg3])

with open("results_6_5.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["insertions", "alg1", "alg2", "alg3"])
    writer.writerows(results)

# Run multiple test, save them in csv file ready to plot in latex
results = []

for ins in range(1,32):
    tester = TestRun(2, 8, 8, 6, ins, verbose=False)
    alg1, alg2, alg3 = tester.run(100)
    results.append([ins, alg1, alg2, alg3])

with open("results_8_8.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["insertions", "alg1", "alg2", "alg3"])
    writer.writerows(results)





