from atmons import *

class BaseModel(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.const = Const()

    def do_experiment(self, experiment='test'):
        config = {
            'v_rot': 600,
            'i_ang': 45,
            'chem': 's1',
            'spec_key': 'wfc',
        }
        ns = build_Neutron_Star(**config)
        shot = ns.atmosphere
        sp_layer = ns.sp_layer

        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        print(f"{'wfc':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")

        return 0


model = BaseModel(name = 'discussion')

model()