# -*- coding: utf-8 -*-
"""
File containing Poles And Zeros (PAZ) and gain for common seismometers.
The instruments must be corrected before to velocity, which is the RESP/SEED
standard.

The seismometer is represented as a dictionary containing the fields:

:type poles: List of Complex Numbers
:ivar poles: Poles of the seismometer to simulate
:type zeros: List of Complex Numbers
:ivar zeros: Zeros of the seismometer to simulate
:type gain: Float
:ivar gain: Gain factor of seismometer to simulate

Currently contained seismometers::

    PAZ_WOOD_ANDERSON

Note, there is only one zero and two poles. That is when simulating the
Wood Anderson, the signal is automatically integrated (most probably to
meter)
"""

# Im Wood-Anderson
# findet aber keine Umrechnung statt, es gibt Länge (displacement) als
# Länge (Ausschlag der "Nadel" (Lichtstrahl auf Photoplatte)) wieder aus.
# D.h. diese 2800 sind einfach ein Vergrößerungsfaktor, deine Einheit
# kommt aus dem, worauf du vorher dein anderes Seismometer runtergerechnet
# hast. (thanks to Christian Sippl)
PAZ_WOOD_ANDERSON = {
    'poles': [-6.2832 - 4.7124j,
              - 6.2832 + 4.7124j],
    'zeros': [0.0 + 0.0j] * 1,
    'gain': 2800
}

#CX stations
PAZ_STS2 = {'normalization_factor': 60077000.0,
            'name': 'GFZ:STS-2/N/g=1500',
            'sensitivity': 629145000.0,
            'normalization_frequency': 1.0,
            'sensor_manufacturer': 'Streckeisen',
            'sensitivity_unit': 'M/S',
            'sensitivity_frequency': 0.02,
            'poles': [(-0.037004 + 0.037016j), (-0.037004 - 0.037016j),
                      (-251.33 + 0j), (-131.04 - 467.29j), (-131.04 + 467.29j)],
            'gain': 60077000.0,
            'zeros': [0j, 0j],
            'sensor_model': 'STS-2/N'}
# sensitivity at some stations 'sensitivity': 2516580000.0


#LVC
PAZ_GURALP = {'normalization_factor': 71367000.0,
              'name': 'GFZ:Guralp_CMG-3T/120F/g=1500',
              'sensitivity': 33500713000.0,
              'normalization_frequency': 0.1,
              'sensor_manufacturer': 'Guralp',
              'sensitivity_unit': 'M/S',
              'sensitivity_frequency': 0.04,
              'poles': [(-0.03701 + 0.03701j), (-0.03701 - 0.03701j),
                        (-197.9 + 197.9j), (-197.9 - 197.9j), (-911.1 + 0j)],
              'gain': 71367000.0, 'zeros': [0j, 0j],
              'sensor_model': 'Guralp_CMG-3T/120F'}
