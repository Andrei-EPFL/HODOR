#!/usr/bin/env python3
import numpy as np
from scipy.special import erf
import warnings
from halotools.empirical_models import occupation_model_template, model_defaults
from halotools.custom_exceptions import HalotoolsError

# The models in this file are modulates of the Zheng07 model, see
# https://halotools.readthedocs.io/en/latest/_modules/halotools/empirical_models/occupation_models/zheng07_components.html

class MWCens(occupation_model_template.OccupationComponent):
  '''
    The central occupation model used by Martin White:
    <Ncen> = 0.5 * [ 1 + Erf( x / sqrt{2} / sigma ) ]
    where x = ln( Mhalo / Mcut )
  '''
  def __init__(self, threshold=model_defaults.default_luminosity_threshold, \
      prim_haloprop_key=model_defaults.prim_haloprop_key, redshift=0, \
      **kwargs):
    '''
      Examples
      --------
      cen_model = MWCens()
      cen_model = MWCens(threshold=-19.5)
      cen_model = MWCens(prim_haloprop_key='halo_m200b')
    '''
    upper_occupation_bound = 1.0
    super(MWCens, self).__init__(gal_type='centrals', \
        threshold=threshold, upper_occupation_bound=upper_occupation_bound, \
        prim_haloprop_key=prim_haloprop_key, **kwargs)
    self.redshift = redshift
    self.param_dict = self.get_default_parameters(self.threshold)

  def mean_occupation(self, **kwargs):
    if 'table' in list(kwargs.keys()):
      mass = kwargs['table'][self.prim_haloprop_key]
    elif 'prim_haloprop' in list(kwargs.keys()):
      mass = np.atleast_1d(kwargs['prim_haloprop'])
    else:
      msg = ('\nYou must pass either a `table` or `prim_haloprop` argument \n'
          'to the `mean_occupation` function of the `MWCens` class.\n')
      raise HalotoolsError(msg)

    logM = np.log(mass)
    inv_sqrt2 = 0.7071067811865475244
    ln10 = 2.30258509299404568402
    mean_ncen = 0.5 * (1.0 + erf((logM - ln10 * self.param_dict['logMcut']) * \
        inv_sqrt2 / self.param_dict['sigma_logM']))
    return mean_ncen

  def get_default_parameters(self, threshold):
    '''
      Best-fit HOD parameters from Martin White
    '''
    param_dict = (
        {'logMcut': 11.5,
        'sigma_logM': 1.0}
        )
    return param_dict


class MWSats(occupation_model_template.OccupationComponent):
  '''
    The satellite occupation model used by Martin White:
    <Nsat> = [ ( Mhalo - k * Mcut ) / M1 ]^alpha
  '''
  def __init__(self, threshold=model_defaults.default_luminosity_threshold, \
      prim_haloprop_key=model_defaults.prim_haloprop_key, redshift=0, \
      modulate_with_cenocc=False, cenocc_model=None, **kwargs):
    upper_occupation_bound = float("inf")
    super(MWSats, self).__init__(gal_type='satellites', \
        threshold=threshold, upper_occupation_bound=upper_occupation_bound, \
        prim_haloprop_key=prim_haloprop_key, **kwargs)
    self.redshift = redshift
    self.param_dict = self.get_default_parameters(self.threshold)

    if cenocc_model is None:
      cenocc_model = MWCens(prim_haloprop_key=prim_haloprop_key, \
          threshold=threshold)
    else:
      if modulate_with_cenocc is False:
        msg = ("You chose to input a `cenocc_model`, but you set the \n"
            "`modulate_with_cenocc` keyword to False, so your "
            "`cenocc_model` will have no impact on the model's behavior.\n"
            "Be sure this is what you intend before proceeding.\n")
        warnings.warn(msg)

    self.modulate_with_cenocc = modulate_with_cenocc
    if self.modulate_with_cenocc:
      try:
        assert isinstance(cenocc_model, \
            occupation_model_templateOccupationComponent)
      except AssertionError:
        msg = ('The input `cenocc_model` must be an instance of \n'
            '`OccupationComponent` or one of its sub-classes.\n')
        raise HalotoolsError(msg)

      self.central_occupation_model = cenocc_model
      self.param_dict.update(self.central_occupation_model.param_dict)

  def mean_occupation(self, **kwargs):
    if self.modulate_with_cenocc:
      for key, value in self.param_dict.items():
        if key in self.central_occupation_model.param_dict:
          self.central_occupation_model.param_dict[key] = value

    if 'table' in list(kwargs.keys()):
      mass = kwargs['table'][self.prim_haloprop_key]
    elif 'prim_haloprop' in list(kwargs.keys()):
      mass = np.atleast_1d(kwargs['prim_haloprop'])
    else:
      msg = ('\nYou must pass either a `table` or `prim_haloprop` argument \n'
          'to the `mean_occupation` function of the `MWSats` class.\n')
      raise HalotoolsError(msg)

    Mcut = 10.**self.param_dict['logMcut']
    M1 = 10.**self.param_dict['logM1']

    mean_nsat = np.zeros_like(mass)
    idx_nonzero = np.where(mass - self.param_dict['k'] * Mcut > 0)[0]
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", RuntimeWarning)
      mean_nsat[idx_nonzero] = \
          ((mass[idx_nonzero] - self.param_dict['k'] * Mcut) / \
          M1)**self.param_dict['alpha']
#    mean_nsat = ((mass - self.param_dict['k'] * Mcut) / \
#        M1)**self.param_dict['alpha']

    if self.modulate_with_cenocc:
      mean_ncen = self.central_occupation_model.mean_occupation(**kwargs)
      mean_nsat *= mean_ncen

    return mean_nsat

  def get_default_parameters(self, threshold):
    '''
      Best-fit HOD parameters from Martin White
    '''
    param_dict = (
        {'logM1': 13.75,
        'logMcut': 11.5,
        'k': 0.5,
        'alpha': 0.5}
        )
    return param_dict


class ContrerasCens(occupation_model_template.OccupationComponent):
  '''
    The central occupation model from Contreras (1301.3497)
  '''
  def __init__(self, threshold=model_defaults.default_luminosity_threshold, \
      prim_haloprop_key=model_defaults.prim_haloprop_key, redshift=0, \
      **kwargs):
    '''
      Examples
      --------
      cen_model = ContrerasCens()
      cen_model = ContrerasCens(threshold=-19.5)
      cen_model = ContrerasCens(prim_haloprop_key='halo_m200b')
    '''
    upper_occupation_bound = 1.0
    super(ContrerasCens, self).__init__(gal_type='centrals', \
        threshold=threshold, upper_occupation_bound=upper_occupation_bound, \
        prim_haloprop_key=prim_haloprop_key, **kwargs)
    self.redshift = redshift
    self.param_dict = self.get_default_parameters(self.threshold)

  def mean_occupation(self, **kwargs):
    if 'table' in list(kwargs.keys()):
      mass = kwargs['table'][self.prim_haloprop_key]
    elif 'prim_haloprop' in list(kwargs.keys()):
      mass = np.atleast_1d(kwargs['prim_haloprop'])
    else:
      msg = ('\nYou must pass either a `table` or `prim_haloprop` argument \n'
          'to the `mean_occupation` function of the `ContrerasCens` class.\n')
      raise HalotoolsError(msg)

    logM = np.log10(mass)
    fac = np.zeros_like(mass)
    idx = (logM < self.param_dict['logMc'])

    fac[idx] = (logM[idx] - self.param_dict['logMc'])/self.param_dict['siga']
    fac[~idx] = (logM[~idx] - self.param_dict['logMc'])/self.param_dict['sigb']

    mean_ncen = self.param_dict['Fb'] * (1 - self.param_dict['Fa']) * \
        np.exp(-0.5 * fac**2) + self.param_dict['Fa'] * (1 + erf(fac))
    return mean_ncen

  def get_default_parameters(self, threshold):
    '''
      Best-fit HOD parameters
    '''
    param_dict = (
        {'Fa': 4.502e-3,
        'Fb': 0.255,
        'logMc': 11.332,
        'siga': 0.1,
        'sigb': 0.323}
        )
    return param_dict


class ContrerasSats(occupation_model_template.OccupationComponent):
  '''
    The satellite occupation model from Contreras (1301.3497)
  '''
  def __init__(self, threshold=model_defaults.default_luminosity_threshold, \
      prim_haloprop_key=model_defaults.prim_haloprop_key, redshift=0, \
      modulate_with_cenocc=False, cenocc_model=None, **kwargs):
    upper_occupation_bound = float("inf")
    super(ContrerasSats, self).__init__(gal_type='satellites', \
        threshold=threshold, upper_occupation_bound=upper_occupation_bound, \
        prim_haloprop_key=prim_haloprop_key, **kwargs)
    self.redshift = redshift
    self.param_dict = self.get_default_parameters(self.threshold)

    if cenocc_model is None:
      cenocc_model = ContrerasCens(prim_haloprop_key=prim_haloprop_key, \
          threshold=threshold)
    else:
      if modulate_with_cenocc is False:
        msg = ("You chose to input a `cenocc_model`, but you set the \n"
            "`modulate_with_cenocc` keyword to False, so your "
            "`cenocc_model` will have no impact on the model's behavior.\n"
            "Be sure this is what you intend before proceeding.\n")
        warnings.warn(msg)

    self.modulate_with_cenocc = modulate_with_cenocc
    if self.modulate_with_cenocc:
      try:
        assert isinstance(cenocc_model, \
            occupation_model_templateOccupationComponent)
      except AssertionError:
        msg = ('The input `cenocc_model` must be an instance of \n'
            '`OccupationComponent` or one of its sub-classes.\n')
        raise HalotoolsError(msg)

      self.central_occupation_model = cenocc_model
      self.param_dict.update(self.central_occupation_model.param_dict)

  def mean_occupation(self, **kwargs):
    if self.modulate_with_cenocc:
      for key, value in self.param_dict.items():
        if key in self.central_occupation_model.param_dict:
          self.central_occupation_model.param_dict[key] = value

    if 'table' in list(kwargs.keys()):
      mass = kwargs['table'][self.prim_haloprop_key]
    elif 'prim_haloprop' in list(kwargs.keys()):
      mass = np.atleast_1d(kwargs['prim_haloprop'])
    else:
      msg = ('\nYou must pass either a `table` or `prim_haloprop` argument \n'
          'to the `mean_occupation` function of the `ContrerasSats` class.\n')
      raise HalotoolsError(msg)

    Mmin = 10.**self.param_dict['logMmin']

    mean_nsat = np.zeros_like(mass)
    mean_nsat = self.param_dict['Fs'] * (1 + erf((np.log10(mass) - \
        self.param_dict['logMmin']) / self.param_dict['deltaM'])) * \
        (mass / Mmin)**self.param_dict['alpha']

    if self.modulate_with_cenocc:
      mean_ncen = self.central_occupation_model.mean_occupation(**kwargs)
      mean_nsat *= mean_ncen

    return mean_nsat

  def get_default_parameters(self, threshold):
    '''
      Best-fit HOD parameters
    '''
    param_dict = (
        {'Fs': 3.502e-3,
        'logMmin': 11.511,
        'deltaM': 0.156,
        'alpha': 0.73}
        )
    return param_dict
