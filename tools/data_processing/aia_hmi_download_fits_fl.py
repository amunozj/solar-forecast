#!/usr/bin/env python

import sys

sys.dont_write_bytecode = True

import datetime
import urllib
from astropy.io import fits
import numpy as np
import json
import pickle
import time
import urllib
import glob
import pandas as pd
# import feather
import os.path
from sunpy.net import jsoc, attrs, vso
import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
import shutil

HMI_tmp_path = '/d1/flares/tmp'

Coeff_file_path = 'aia_degrade.txt'
Coeff_data = pd.read_csv('aia_degrade.txt', header=0)
Coeff_init_date = datetime.datetime.strptime(Coeff_data['Date'].values[0][:-4], '%Y-%m-%dT%H:%M:%S')


def get_Coeff_index(date_s):
    date = datetime.datetime.strptime(date_s, '%Y%m%d')
    delta = date - Coeff_init_date
    delta = delta.days + 1
    return int(delta)


flaresCSV = pd.read_csv('flares_CMX_Flux_2010_2017.csv', header=0)
flares = flaresCSV.values[flaresCSV.values[:, 11] >= 1e-6, :]

aia_ls = [94, 131, 171, 193, 211, 304, 335, 1600]
columns = aia_ls.copy()

columns.append('HMI')

count = 0;
it = np.nditer(flares[:, 11], flags=['f_index'])
while not it.finished:
    date_s = '%4d%02d%02d_%02d%02d' % (
    flares[it.index, 0], flares[it.index, 1], flares[it.index, 2], flares[it.index, 9], flares[it.index, 10])
    date = datetime.datetime.strptime(date_s, '%Y%m%d_%H%M')
    # print date
    date -= datetime.timedelta(seconds=3600 * 6)
    # print date
    date_s = datetime.datetime.strftime(date, '%Y%m%d_%H%M')
    date_e = date + datetime.timedelta(seconds=2*720)

    print('{0}: From {1} to {2}'.format(count,date_s, date_e))

    client = jsoc.JSOCClient()
    response = client.search(jsoc.attrs.Time(date.isoformat(), date_e.isoformat()), jsoc.attrs.Series('hmi.M_720s'),
                             jsoc.attrs.Notify('amunozj@boulder.swri.edu'))

    if not os.path.exists(HMI_tmp_path + '/%05d' % count):
        os.makedirs(HMI_tmp_path + '/%05d' % count)

        print(response)

        try:
            if response.table.columns['INSTRUME'][0] != 'HMI_SIDE1':
                continue
        except KeyError:
            if len(response.table) == 0:
                continue
        try:
            res = client.fetch(response, path=HMI_tmp_path + '/%05d' % count)
            res.wait(progress=True)
        except BaseException:
            continue

        # Recasting date to match HMI's availability

        hmi_files = os.listdir(HMI_tmp_path + '/%05d' % count)
        date_s = hmi_files[0][11:24]
        date = datetime.datetime.strptime(date_s, '%Y%m%d_%H%M')

        print(date_s, date)

        # print date_s[0:8]
        inx = get_Coeff_index(date_s[0:8])

        my_file = '/d1/flares/FlaresHMI/XrF%05.2f_AIA8ch_HMI_' % (np.log10(flares[it.index, 11]*1e10)) + date_s + '_1k_006h.fits'
        if not os.path.isfile(my_file):
            print(my_file)

            # data = []
            is_complete = 1

            new_hdul = fits.HDUList()
            for aia_l in aia_ls:

                print(date.year, date.month, date.day, date.hour, date_s, aia_l)
                url = "http://jsoc.stanford.edu/data/aia/synoptic/%d/%02d/%02d/H%02d00/AIA%s_%04d.fits" % (
                date.year, date.month, date.day, date.hour, date_s, aia_l)
                print(url)

                try:
                    response = urllib.request.urlopen(url)
                except:
                    is_complete = 0
                    break

                if response.getcode() == 200:

                    while True:
                        print('while loop')
                        try:
                            chromosphere_image = fits.open(url, cache=False)
                            chromosphere_image.verify("fix")
                            break
                        except:
                            print('Timed out, trying again')

                else:
                    is_complete = 0
                    break

                exptime = chromosphere_image[1].header['EXPTIME']
                if exptime == 0:
                    chromosphere_image.close()
                    is_complete = 0
                    break

                if chromosphere_image[1].header['QUALITY'] != 0:
                    chromosphere_image.close()
                    is_complete = 0
                    break

                image_hud = fits.ImageHDU(
                    data=np.array(chromosphere_image[1].data / exptime / Coeff_data[str(aia_l)].values[inx], dtype=np.float32),
                    header=chromosphere_image[1].header)

                new_hdul.append(image_hud)
                chromosphere_image.close()

            if is_complete:



                print(HMI_tmp_path + '/' + hmi_files[0])
                hmiFile = fits.open(HMI_tmp_path + '/%05d' % count + '/' + hmi_files[0], cache=False)
                hmiFile.verify("fix")

                #Create AIA map to capture the frame of reference
                aiaMap = sunpy.map.Map(chromosphere_image[1].data, chromosphere_image[1].header)

                #Create HMI map to interpolate to AIA field of view
                hmiMap = sunpy.map.Map(hmiFile[1].data, hmiFile[1].header)

                #Set values outside the solar radius to 0
                x, y = np.meshgrid(*[np.arange(v.value) for v in hmiMap.dimensions]) * u.pixel
                hpcCoords = hmiMap.pixel_to_world(x, y)
                r = np.sqrt(hpcCoords.Tx ** 2 + hpcCoords.Ty ** 2) / hmiMap.rsun_obs
                hmiMap.data[r>1] = 0

                #Rotate to account for P angle and re-scale to match AIA scale
                hmiMap = hmiMap.rotate(rmatrix=hmiMap.rotation_matrix, scale=hmiMap.scale[0] / aiaMap.scale[0],
                                       recenter=True)

                #Find the pixel coordinates of the AIA bottom-left in HMI
                aiaBottomLeft = SkyCoord(aiaMap.bottom_left_coord.Tx, aiaMap.bottom_left_coord.Ty, frame=hmiMap.coordinate_frame)
                pxBottomLeft = hmiMap.world_to_pixel(aiaBottomLeft)
                pxBottomLeft = np.array([int(np.round(pxBottomLeft[0]/u.pixel)), int(np.round(pxBottomLeft[1]/u.pixel))])
                pxTopRight = pxBottomLeft + 1024

                #Get only the pixels that match the AIA field of view
                hmiMap = hmiMap.submap(pxBottomLeft*u.pixel, top_right=pxTopRight*u.pixel)

                # fig = plt.figure()
                # aiaMap.plot(cmap=plt.cm.Greys_r)
                # plt.show()
                #
                # fig = plt.figure()
                # hmiMap.plot(cmap=plt.cm.Greys_r)
                # plt.show()


                hmiMap.meta['history'] = 'Remapped to AIA FOV'
                hmiMap.meta['comment'] = 'Removed history and comments to avoid astropy errors'
                hmiMap.meta[''] = 'NaN'

                # Append HMI magnetogram
                image_hud = fits.ImageHDU(
                    data=np.array(hmiMap.data, dtype=np.float32),
                    header=fits.Header(hmiMap.meta))


                new_hdul.append(image_hud)
                # data.append(hmiMap.data.flatten())

                # data = np.array(data, dtype=np.float32).transpose()
                # # print data.shape



                # df = pd.DataFrame(data, columns=columns)
                # print df

                # np.save(file('AIA_data_2014/%05d_%s_AIA_%02d_1024_1024.dat'%(count,date_s,len(aia_ls)),'w'),data)

                # np.save(file('/home/solardynamo/AIA_data_Flares/%05d_AIA'%(count) + date_s + '_8chnls_1024_012h.hdr','w'),chromosphere_image[1].header)
                # np.save(file('/home/solardynamo/AIA_data_Flares/%05d_12m_AIA'%(count) + date_s + '_8chnls_1024.dat','w'),data)
                # df = feather.write_dataframe(df, '/d1/flares/FlaresHMI/%05d_AIA' % (
                # count) + date_s + '_8chnls_1024_004h.fthr')

                new_hdul.writeto('/d1/flares/FlaresHMI/XrF%05.2f_AIA8ch_HMI_' % (np.log10(flares[it.index, 11]*1e10)) + date_s + '_1k_006h.fits', output_verify='ignore')

                hmiFile.close()



    count += 1

    it.iternext()

