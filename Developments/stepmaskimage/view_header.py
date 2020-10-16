from darepype.drp import DataFits
from astropy.io import fits

config = '/Users/josh/pipeline/pipeline/Developments/stepwebastrometry/pipeconf_stonedge_auto.txt'
fp = '/Users/josh/Desktop/pipeline_test/data/M5_r-band_60s_bin2_200711_053415_itzamna_seo_0_RAW_TABLE.fits'

fts = DataFits(config=config)
fts.load(fp)
# print(repr(fits.HDUList(file=fp)))
# fts.header['RA'] = 0
# fts.header['Dec'] = 0
print(repr(fts.header))
print(repr(fts.image))
# print(repr(fts.table))