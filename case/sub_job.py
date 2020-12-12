import subprocess
import numpy as np


def gen_seeds(path):
    nseeds = 5
    seeds = np.random.randint(999, 99999, nseeds)
    with open('{}/seeds.txt'.format(path), 'w') as f:
        for i in range(nseeds):
            f.write('variable seed{} equal {}\n'.format(i+1, seeds[i]))

def create_folder(rbc_tp, n=5):
    for i in range(n):
        subprocess.call('mkdir {}{}'.format(rbc_tp, i), shell=True)
        subprocess.call('cp {}_temp/* {}{}'.format(rbc_tp, rbc_tp, i), shell=True)
        gen_seeds('{}{}'.format(rbc_tp, i))
        subprocess.call('cd {}{} && mkdir data data/traj && cp restart* data/ && sbatch adr.sh'.format(\
                        rbc_tp, i), shell=True)

def main():
    create_folder('nrbc', 5)
    create_folder('drbc', 5)


if __name__ == "__main__":
    main()
