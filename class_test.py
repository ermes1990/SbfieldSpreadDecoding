from tqdm import tqdm
from math import floor
import csv


class TestRun:
    def __init__(self, q, r, k, deletions, insertions):
        self.q = q
        self.r = r
        self.k = k
        self.n = r*k
        self.deletions = deletions
        self.insertions = insertions

        # Fields
        self.Fq = GF(q, 'z')
        self.Fqn = GF(self.q^self.n, 'y')
        self.Fqk, self.embedding = self.Fqn.subfield(self.k, map=True)
        self.alpha = self.Fqk.gen()

        # Ambient vector space
        self.V = VectorSpace(self.Fq, self.n)

    # --- helper: multiply subfield element into big field ---
    # Input: space R,  z in Fqn
    # Output: zR = {zr | r in R}.
    def mul_subspace(self, R, z):
        R_basis = [self.Fqn(x.list()) for x in R.basis()]
        R_z_basis = [self.embedding(z) * x for x in R_basis]
        return self.V.subspace(R_z_basis)

    # Input: R_exp, R (subspaces)
    # Output R_exp + zR
    def expand_subspace(self, R_exp, R):
        # pick random z not in prime subfield
        z = self.Fqk.random_element()
        while z in self.Fqk.prime_subfield():
            z = self.Fqk.random_element()
        R_z = self.mul_subspace(R, z)
        R_exp = R_exp + R_z
        return R_exp

    #Input: R_cont, R (subspaces)
    #Output: R_cont meet zR 
    def reduce_subspace(self, R_cont, R):
        z = self.Fqk.random_element()
        while z in self.Fqk.prime_subfield():
            z = self.Fqk.random_element()
        R_z = self.mul_subspace(R, z)
        R_cont = R_cont.intersection(R_z)
        return R_cont

    
    # Output:  R=U+B, aFqk
    # U is a subspace of aFqk of dimension k - deletion
    # B is a random subspace of dimension insertion 
    # aFqk is the original message we want to recover
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

    
    # Expand Reduce algorithm
    def decode_ER(self, R):
        # Set first expansion size.
        #max_exp = min(ceil((self.V.dimension()) / R.dimension()) - 1, self.k - 1)
        max_exp = min(ceil((self.V.dimension() - self.k -self.r + 1) / (R.dimension()-2)) - 1, self.k - 1)
        #max_exp = min(ceil((self.V.dimension()-self.r) / R.dimension()) - 1, self.k - 1)
        R_exp = self.V.subspace([self.Fqn(x.list()) for x in R.basis()])
        #print(f"max_exp: {max_exp}")
        # Expand
        if max_exp > 0:
            for _ in range(max_exp - 1): # diemsnion starts from r with 0 expansion, doing (e-1) expansions will result in dimension ~ er
                R_exp = self.expand_subspace(R_exp, R)
        res = R_exp
        #print(f"R_expanded: {R_exp.dimension()}")
        # Reduce and keep the last non-zero value
        lim = self.Fqk.degree() + 5
        res_old = res
        while res.dimension() > 0 and lim >= 0:
            res_old = res
            #print(f"{res_old.dimension()}")
            res = self.reduce_subspace(res, R_exp)
            if res.dimension() == self.k:
                return res
            lim -= 1
        return res
    
    # Expand Reduce Expand algorithm
    def decoding_ERE(self, R):
        # Set first expansion size.
        R_exp = self.V.subspace([self.Fqn(x.list()) for x in R.basis()])
        #print(max_exp)
        # Expand
        max_exp = min(ceil((self.V.dimension()-self.r + 1) / R.dimension()) - 2, self.k - 1)
        if max_exp > 0:
            for _ in range(max_exp - 1): # diemsnion starts from r with 0 expansion, doing (e-1) expansions will result in dimension ~ er
                R_exp = self.expand_subspace(R_exp, R)
        res = R_exp
        # Reduce and keep last non-zero space
        res_old = res
        lim = self.Fqk.degree() + 5
        while res.dimension() > 0 and lim >= 0:
            res_old = res
            #print(res_old.dimension())
            res = self.reduce_subspace(res, R_exp)
            lim -= 1
        res = res_old
        # Expand until raching dimension k.
        while res.dimension() < self.Fqk.degree() and res.dimension() > 0:
            res = self.expand_subspace(res, res)
        return self.expand_subspace(res, res)

    # Filter algorith to use in filtered ERE
    def filter(self, R, dim_extract= 5):
        # Set expansion size
        max_exp1 = min(ceil((self.V.dimension() -self.r + 1) / R.dimension()) - 2, self.k - 1)
        max_exp2 = min(floor((self.V.dimension() - self.k)/(R.dimension() - 2)) - 1, self.k -1)
        max_exp = min(max_exp1, max_exp2)
        #max_exp = min(ceil((self.V.dimension()-self.r) / R.dimension()) - 1, self.k - 1)
        #print(f"ext exp:{max_exp}")
        C = self.V.subspace([]) # Initialize filter space to 0 dimensional subspace
        # Start collecting intersections in  filter space until reach threshold dimension dim_extract.
        while C.dimension() < dim_extract:
            # choose two Fq- l.i. elements in Fqn  
            y = self.Fqk.random_element()
            z = self.Fqk.random_element()
            while z==0 or y/z in self.Fq: # Check l.i.
                y = self.Fqk.random_element()
                z = self.Fqk.random_element()
            # multiply received space R -> yR, zR
            yR = self.mul_subspace(R, y)
            zR = self.mul_subspace(R, z)
            # Expand zR
            for _ in range(max_exp - 1):
                zR = self.expand_subspace(zR, R)
            # Intersect yR, Expanded(zR)
            Extract = yR.intersection(zR)
            C = C + Extract # add extraction to filter
        return C
    
    
    # Run test on the setting of self
    # keep track of how many succesful recoveries for each algorithm.
    # alg1 = ER, alg2 = ERE, alg3 = Filtered ERE
    # Params:
    # trials: numer of runs
    # mask: integer to binary chose which test to run (e.g. mask=6 -> 110: run alg 3 and 2; mask = 1 -> 001 run only alg1
    # extract: Set size of filter before returning
    # max_ERE: Max number of to times repeat ERE when the answer has dimension != k (i.e. we know the algorithm failed so we can repeat)
    def run(self, trials=100, mask=7, extract=5, max_ERE=5):
        print(f"Test: ({self.deletions},{self.insertions})")
        counter_alg1 = 0
        counter_alg2 = 0
        counter_alg3 = 0
        for _ in tqdm(range(trials)):
            R, VA = self.corrupt_encode_random()
            if mask & 1:
                for _ in range(5):
                    dec_A = self.decode_ER(R)
                    if dec_A.dimension() == self.k:
                        if dec_A == VA:
                            counter_alg1 += 1
                        break
            if mask & 2:
                for _ in range(max_ERE):
                    dec_A = self.decoding_ERE(R)
                    if dec_A.dimension() == self.k:
                        if dec_A == VA:
                            counter_alg2 += 1
                        break
            if mask & 4:
                C = self.filter(R, dim_extract=extract)
                for _ in range(5):
                    dec_A = self.decoding_ERE(C)
                    #print(dec_A.dimension())
                    #print(VA.is_subspace(dec_A))
                    if dec_A.dimension() == self.k:
                        if dec_A == VA:
                            counter_alg3 += 1
                        break
        print(f"({counter_alg1}, {counter_alg2}, {counter_alg3}) \n")
        return counter_alg1, counter_alg2, counter_alg3

    # Try ERE different number of times.
    # Save when first success
    # Report how many times it was successful at different iteration
    # Input: 
    #       max_ERE := nomber of different max number of rtepetition allowed.  example [1,5,10]
    #       trials: number of test to run.
    # Output: successful count after iteration idicated by max_ERE.
    def run_ERE(self, max_ERE, trials=100):
        print(f"Test_ERE: ({self.deletions},{self.insertions})")
        counter = [0]*len(max_ERE)
        for _ in tqdm(range(trials)):
            # Create corrupted codeword to attempt to decode
            R, VA = self.corrupt_encode_random()
            success_at = None
            # Start decoding
            for attempt in range(max(max_ERE)):
                dec_A = self.decoding_ERE(R)
                # When dimension is not k, we can do another attempt
                if dec_A.dimension() != self.k:
                    continue
                # When k is either the correct or we got an error, no more repetition.    
                if dec_A.dimension() == self.k:
                    # If correct we save the attempt number at which we got success
                    if dec_A == VA:
                       success_at = attempt
                    break # we end this attempt either it was successful or we got a mistake.
            
            # Increment counters corresponding below attempt number
            if success_at != None:    
                for j in range(len(counter)):
                    if success_at < max_ERE[j]:
                        counter[j] += 1 
        print(f"{counter}")
        return counter

# Set seed to reproduce same outcome
set_random_seed(67)

# Run multiple test, save them in csv file ready to plot in latex
results = []

for ins in range(1,16):
    tester = TestRun(2, 6, 5, 3, ins)
    alg1, alg2, alg3 = tester.run(1000,mask=7, extract=3, max_ERE=5)
    results.append([ins, alg1, alg2, alg3])
    

with open("results_6_5.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["insertions", "alg1", "alg2", "alg3"])
    writer.writerows(results)


# Set seed to reproduce same outcome
set_random_seed(67)

# Run multiple test, save them in csv file ready to plot in latex
results = []

for ins in range(16,19):
    tester = TestRun(2, 8, 8, 6, ins)
    alg1, alg2, alg3 = tester.run(100,extract=3, mask=7, max_ERE=5)
    results.append([ins, alg1, alg2, alg3])

with open("results_8_8.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["insertions", "alg1", "alg2", "alg3"])
    writer.writerows(results)


# Set seed to reproduce same outcome
set_random_seed(67)

results = []
# Set max iteration to compare
max_ERE = [1,5,10,25,100]
for ins in range(1,32):
    tester = TestRun(2, 8, 8, 6, ins)
    results.append([ins] + tester.run_ERE(max_ERE, trials=100))

# Write file with results
with open("results_ERE.csv", "w", newline="") as f:
    labels = ["insertions"] + [f"ERE_{m}" for m in max_ERE]
    writer = csv.writer(f)
    writer.writerow(labels)
    writer.writerows(results)




