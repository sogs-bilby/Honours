import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import sys
# weird correcting things

def convert(filename, instr_index, type_time="jd"):
    # load data
    data = np.loadtxt(f"data/{filename}")
    raw_times = data[:, 0]

    if data.shape[1] < 4:
        new_col = np.repeat([instr_index], data.shape[0]).reshape(-1, 1)
        data = np.append(data, new_col, axis=1)

    print(data)

    if type_time != "bjd":
        # define coordiantes of star in sky
        star = SkyCoord(ra='08h42m25.12195s', dec='+04d34m41.1457s', frame='icrs')

        # location of observatory. using the general geocentre doesn't seem to make any real difference
        # over defining the location as the actual instrument's location

        observatory = EarthLocation.of_site("Observatoire de Haute Provence")
        #observatory = EarthLocation.from_geocentric(0, 0, 0, unit="m")

        # astropy time object
        t = Time(raw_times, format='jd', scale='utc', location=observatory)

        if type_time == "hjd":
            # convert HJD to JD
            t -= t.light_travel_time(star, kind="heliocentric")

        # convert to BJD in Barycentric Dynamical Time (TDB) to account for relativity
        t_tdb = t.tdb

        # get light travel time correction to barycentre
        ltt_bary = t.light_travel_time(star, kind="barycentric")

        # do correction to get final value for bjd in tdb
        bjd_tdb = t_tdb + ltt_bary

        # save
        data[:, 0] = bjd_tdb.value

    np.savetxt(f"data/{filename}_corrected.txt", data, fmt="%.6f")

    print(f"Successfully converted {filename} to BJD")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 clean_data.py <filename> <instrument index> <time type>")
        sys.exit(1)

    convert(sys.argv[1], int(sys.argv[2]), sys.argv[3])