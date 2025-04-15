from atmons import *

config, grid =  loader()

# config['spec_key'] = 'be'
config['spec_key'] = 'wfc'
grid['n_nu'] = 300
grid['n_theta'] = 30
grid['n_phi'] = 30

dumper('paper', config, grid)
config, grid = loader("paper")


NS = NeurtonStar(config, grid)

burst = NS.burst()
counter = 1
Lum = []
import matplotlib.pyplot as plt

for shot in burst:
    print(f"{counter} | {shot.lum}")
    plt.plot(log(shot.E_null), log(shot.B_real) - 36,)
    counter += 1
    if counter>1:
        break

plt.show()