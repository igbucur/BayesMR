## Section 4 - Bayesian model averaging and inference
## 4.3 - Example: Instrumental Variable Setting

source('R/write_configuration_file.R')
source('R/PolyChord_interface.R')

write_configuration_file(
  config_file_name = 'ini/smmr_4_4_dir.ini',
  SS_filename = 'data/near_LCD_SSc.txt',
  sigmaG_filename = 'data/near_LCD_SiG.txt',
  num_instruments = 1,
  num_observations = 10000,
  PolyChord_control = list(num_live_points = 1000)
)

write_configuration_file(
  config_file_name = 'ini/smmr_4_4_rev.ini',
  SS_filename = 'data/near_LCD_SSr.txt',
  sigmaG_filename = 'data/near_LCD_SiG.txt',
  num_instruments = 1,
  num_observations = 10000,
  PolyChord_control = list(num_live_points = 1000)
)


system('bin/polychord_MR ini/smmr_4_4_dir.ini')
system('bin/polychord_MR ini/smmr_4_4_rev.ini')


print(derive_polychord_Bayes_factor('chains/smmr_4_4_dir.stats', 'chains/smmr_4_4_rev.stats'))


