import mosaik_api_v3
import numpy as np
from core.external_conditions import ExternalConditions
from core.simulation_time import SimulationTime
from core.energy_supply.energy_supply import EnergySupply, EnergySupplyConnection
from  core.energy_supply.pv import PhotovoltaicSystem

META = {
    'models': {
        'PVSystem': {
            'public': True,
            'params': [
                'peak_power', 'ventilation_strategy', 'pitch', 'orientation',
                'base_height', 'height', 'width', 'inverter_peak_power', 'inverter_is_inside'
            ],
            'attrs': ['energy_produced', 'energy_lost'],
        },
    }
}



class PVSimulator(mosaik_api_v3.Simulator):
    def __init__(self):
        super().__init__(META)
        self.entities = {}
        self.next_eid = 1
        self.simulation_time = SimulationTime(0, 8, 1)
        self.external_conditions = None  # Will be set later
        self.energy_supply = EnergySupply("electricity", self.simulation_time)
        self.energy_supply_conn = self.energy_supply.connection("pv generation without shading")

    def init(self, sid, air_temperatures, wind_speeds, diffuse_horizontal_radiation, direct_beam_radiation,
             solar_reflectivity_of_ground, latitude, longitude, timezone, start_day, end_day, time_series_step,
             january_first, daylight_savings, leap_day_included, direct_beam_conversion_needed, shading_segments):
        self.external_conditions = ExternalConditions(
            self.simulation_time,
            air_temperatures,
            wind_speeds,
            diffuse_horizontal_radiation,
            direct_beam_radiation,
            solar_reflectivity_of_ground,
            latitude,
            longitude,
            timezone,
            start_day,
            end_day,
            time_series_step,
            january_first,
            daylight_savings,
            leap_day_included,
            direct_beam_conversion_needed,
            shading_segments
        )
        return self.meta

    def create(self, num, model, **params):
        entities = []
        for _ in range(num):
            eid = self.next_eid
            self.next_eid += 1
            pv_system = PhotovoltaicSystem(
                params['peak_power'],
                params['ventilation_strategy'],
                params['pitch'],
                params['orientation'],
                params['base_height'],
                params['height'],
                params['width'],
                self.external_conditions,
                self.energy_supply_conn,
                self.simulation_time,
                [],  # Shading, can be parameterized as needed
                params['inverter_peak_power'],
                params['inverter_is_inside']
            )
            self.entities[eid] = pv_system
            entities.append({'eid': eid, 'type': model})
        return entities

    def step(self, time, inputs): # inputs not used since it relies on its external conditions?
        for eid, entity in self.entities.items():
            energy_produced, energy_lost = entity.produce_energy()
            self.data[eid] = {
                'energy_produced': energy_produced,
                'energy_lost': energy_lost,
            }
        return time + 1

    def get_data(self, outputs):
        data = {}
        for eid, attrs in outputs.items():
            data[eid] = {attr: self.data[eid][attr] for attr in attrs}
        return data

if __name__ == '__main__':
    mosaik_api_v3.start_simulation(PVSimulator())