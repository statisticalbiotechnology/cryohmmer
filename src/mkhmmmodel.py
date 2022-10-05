import math

inputFile = "../data/output_aa_probabilities.aap"
delta = 0.05


def printModel(model, conf, id):
    ingress = f'''HMMER3/f [3.2 | April 2018]
NAME  {id}
LENG  {len(model)}
ALPH  amino
RF    no
MM    no
CONS  yes
CS    no
MAP   yes
STATS LOCAL MSV       -9.9014  0.70957
STATS LOCAL VITERBI  -10.7224  0.70957
STATS LOCAL FORWARD   -4.1637  0.70957
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.36553  4.52577  2.96709  2.70473  3.20818  3.02239  3.41069  2.90041  2.55332  2.35210  3.67329  3.19812  3.45595  3.16091  3.07934  2.66722  2.85475  2.56965  4.55393  3.62921
          2.68640  4.42247  2.77497  2.73145  3.46376  2.40504  3.72516  3.29302  2.67763  2.69377  4.24712  2.90369  2.73719  3.18168  2.89823  2.37879  2.77497  2.98431  4.58499  3.61525
          0.57544  1.78073  1.31293  1.75577  0.18968  0.00000        *
'''
    with open(f"../hmm/{id}.hmm","w") as of:
        print(ingress, end='', file=of)
        for ix in range(len(model)):
            ## Match
            print(f'{ix+1: >7} ', end='', file=of)
            low, low_aa = 1.e10, ''
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                if model[ix][aa]<low:
                    low, low_aa = model[ix][aa], aa
                print(f'  {model[ix][aa]: >.5f}', end='', file=of)
            print(f'{ix+1: >7} {low_aa} - - -', file=of)
            ## Insert
            print('        ', end='', file=of)
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                print(f'  {model[ix][aa]: >.5f}', end='', file=of)
            print('', file=of)
            ## Transitions
            # m->m     m->i     m->d     i->m     i->i     d->m     d->d
            mm = max(conf[ix]-delta, 0.5)
            print('        ', end='', file=of)
            # m->m
            print(f'  {-math.log(mm): >.5f}', end='', file=of)
            # m->i
            print(f'  {-math.log((1.-mm)/2.): >.5f}', end='', file=of)
            # m->d
            print(f'  {-math.log((1.-mm)/2.): >.5f}', end='', file=of)
            # i->m
            print(f'  {-math.log(1.-delta): >.5f}', end='', file=of)
            # i->i
            print(f'  {-math.log(delta): >.5f}', end='', file=of)
            # d->m
            print(f'  {-math.log(1-delta): >.5f}', end='', file=of)
            # d->d
            print(f'  {-math.log(delta): >.5f}', end='', file=of)
            print('', file=of)
        print('//', end='', file=of)


with open(inputFile) as file:
    line = file.readline()
    while line.startswith("======="):
        id_line = file.readline().rstrip()
        print(id_line)
        chain = int(id_line.split(":")[1])
        id=f"chain{chain}"
        file.readline() # "Not pruned"
        len_line = file.readline().rstrip().split(": ")[1]
        model_len = int(len_line)
        conf_line = file.readline().rstrip().split(":")[1]
        conf = [ float(flt)/100. for flt in conf_line.split(",") ]
        dummy = file.readline() # "Amino acid probability per residue:"
        model = [{} for _ in range(model_len)]
        try:
            for line in file:
                residue, prob_line = line.rstrip().split(':')
                logprob = [ -math.log(float(flt)) for flt in prob_line.split(",") ]
                for ix in range(model_len):
                    model[ix][residue] = logprob[ix]
        except ValueError:
            pass
        printModel(model, conf, id)
