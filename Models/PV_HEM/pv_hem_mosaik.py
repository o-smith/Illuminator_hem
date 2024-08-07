import itertools
import mosaik_api
import Models.PV_HEM.pv_hem_model as PV_hem_model
"""try:
    import Models.PV_HEM.pv_model as PV_model
except ModuleNotFoundError:
    import pv_model as PV_model
else:
    import Models.PV.pv_model as PV_model"""
import pandas as pd
import itertools
from core.simulation_time import SimulationTime
from core.external_conditions import ExternalConditions


meta = {
    'type': 'event-based', #if reading from a csv file then it is time based
    'models': {
        'PVset': {
            'public': True,
            'params': ['peak_power', 'ventilation_strategy', 'pitch', 'orientation',
                       'base_height', 'height', 'width', 'shading', 'inverter_peak_power',
                        'inverter_is_inside'],
            # and are attrs the specific outputs we want from the code? to connect with other models
            'attrs': ["air_temperatures",
                        "wind_speeds",
                        "wind_directions",
                        "diffuse_horizontal_radiation",
                        "direct_beam_radiation",
                        "solar_reflectivity_of_ground",
                        'pv_gen',
                        'total_irr'
                    ]
        },
    },
}

class PvHemAdapter(mosaik_api.Simulator):
    def __init__(self):
        super(PvHemAdapter, self).__init__(meta)
        self.eid_prefix='pv_'
        self.entities = {}  # every entity that we create of PV gets stored in this dictionary as a list
        self.mods = {}
        self._cache = {}  #we store the final outputs after calling the python model (#PV1) here.

    def init(self, _, time_resolution = .15):
        self.time_resolution = time_resolution
        return self.meta

    def create(self, num, model, sim_start=0, **model_params):
        
        self.start = pd.to_datetime(sim_start)
        entities = []
        for i in range (num):
            eid = '%s%d' % (self.eid_prefix, i)

            # we are creating an instance for PV and call the python file for that. **model_params refers to the
            # parameters we have mentioned above in the META. New instance will have those parameters.
            model_instance = PV_hem_model.PhotovoltaicSystem(**model_params)
            
            self.entities[eid] = model_instance
            entities.append({'eid': eid, 'type': model})
        return entities
    

    def step(self, time, inputs, max_advance):
        # in this method, we call the python file at every data interval and perform the calculations.
        current_time = (self.start + pd.Timedelta(time * self.time_resolution,
                                                  unit='seconds'))  # timedelta represents a duration of time
        print('from pv %%%%%%%%%', current_time)
        for eid, attrs in inputs.items():
            v = []  # we create this empty list to hold all the input values we want to give since we have more than 2
            for _, vals in attrs.items():
                u = list(vals.values())
                v.append(u)

            v_merged = list(itertools.chain(*v))
            self._cache[eid] = self.entities[eid].connect(v_merged[0], v_merged[1], v_merged[2], v_merged[3],
                                                          v_merged[4], v_merged[5])

    def get_data(self, _):
        data = {}
        for eid in self.entities.keys():
            data[eid] = {}
            # print("HERE IT IS:")
            print(self._cache[eid])
            data[eid]["pv_gen"] = self._cache[eid]["pv_gen"]
            print(f'PV GENERATED = {self._cache[eid]["pv_gen"]}')
            data[eid]["total_irr"] = 10.0
        return data
