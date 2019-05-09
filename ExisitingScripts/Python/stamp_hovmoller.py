import iris
import iris.analysis.calculus
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from plotting_tools import *
from storm_tools import *
from cube_tools import *
import matplotlib.cm as mpl_cm
import iris.plot as iplt
from vorticity import calc_vrt_spherical


def load_ens_members(em, fpath, x0, y0):
    md = 1203
    TT = 12
    stash_code1 = 'm01s15i201'
    stash_code2 = 'm01s15i202'
    data_constraint1 = iris.AttributeConstraint(STASH=stash_code1)
    data_constraint2 = iris.AttributeConstraint(STASH=stash_code2)
    p_constraint = iris.Constraint(pressure=850)
    phi_interval = np.pi / 8
    ranges = np.arange(0, 500, 5)
    phis = np.arange(0, 2 * np.pi, phi_interval)
    df = fpath + '{0:02d}.pp'.format(em)
    u = iris.load_cube(df, data_constraint1).extract(p_constraint)
    v = iris.load_cube(df, data_constraint2).extract(p_constraint)

    i = 0

    for u_slc, v_slc in zip(u.slices(['latitude', 'longitude']), v.slices(['latitude', 'longitude'])):
        if i < 2:
            minlat = y0[0] - 5
            maxlat = y0[0] + 5
            minlon = x0[0] - 5
            maxlon = x0[0] + 5
        else:
            minlat = y0[i / 2 - 1] - 5
            maxlat = y0[i / 2 - 1] + 5
            minlon = x0[i / 2 - 1] - 5
            maxlon = x0[i / 2 - 1] + 5

        b_constraint = box_constraint(minlat, maxlat, minlon, maxlon)
        u_box = u_slc.extract(b_constraint)
        v_box = u_slc.extract(b_constraint)
        vrt = calc_vrt_spherical(u_box, v_box)
        cenlat, cenlon, max_val = max_vals(vrt)
        for r in ranges:
            #print r
            for phi in phis:
                #print phi
                xpoi = cenlon + 0.009 * r * np.cos(phi)  # Check the 0.009
                ypoi = cenlat + 0.009 * r * np.sin(phi)
                new_point = [('latitude', ypoi), ('longitude', xpoi)]
                newu = u_slc.interpolate(
                    new_point, iris.analysis.Linear()).data
                newv = v_slc.interpolate(
                    new_point, iris.analysis.Linear()).data
                aziv = -newu * np.sin(phi) + newv * np.cos(phi)
                radu = newu * np.cos(phi) + newv * np.sin(phi)
                if phi == 0:
                    vazi = aziv
                    urad = radu
                else:
                    vazi = np.append(vazi, aziv)
                    urad = np.append(urad, radu)

            if r == 0:
                Vazi = np.mean(vazi)
                Urad = np.mean(urad)
            else:
                Vazi = np.append(Vazi, np.mean(vazi))
                Urad = np.append(Urad, np.mean(urad))
        #print i
        if i == 0:
            v_azi_all = Vazi
            u_rad_all = Urad
        else:
            v_azi_all = np.dstack((v_azi_all, Vazi))
            u_rad_all = np.dstack((u_rad_all, Urad))

        i = i + 1

    return v_azi_all


def max_vals(cube):
    # Find latitude/longitude coordinates of maximum value of data in cube.
    index = np.argmax(cube.data)
    indices = np.unravel_index(index, cube.data.shape)
    # Find the coordinates for this central position
    cenlat = cube.coord('latitude').points[indices[0]]
    cenlon = cube.coord('longitude').points[indices[1]]
    maxval = cube.data[indices]
    # print maxval
    return cenlat, cenlon, maxval


def plot_hovmoller(v_azi, outfile):
    ranges = np.arange(0, 500, 5)
    fig = plt.figure()

    data = v_azi[0]
    data = np.swapaxes(data, 0, 1)
    #print ranges.shape
    #print data.shape
    times = np.arange(41)
    fig = plt.figure(1)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Radius (km)', fontsize=18)
    ax.set_ylabel('Forecast time (Z)', fontsize=18)
    hovmol = ax.pcolormesh(ranges, times, data, cmap='jet')
    cbar = plt.colorbar(hovmol)
    cbar.set_label('Azimuthal velocity (ms$^{-1}$)', size=14)
    cbar.ax.tick_params(labelsize=14)
    # plt.show()
    plt.savefig(outfile, dpi=500)
    plt.close()


def main():
    monthday = [1203]
    times = [12]
    ensemble_members = np.arange(12)
    for md in monthday:
        for TT in times:
            fpath = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/4p4/{0:04d}_{1:02d}Z/wind_plev/uvplev_4p4_{0:04d}_{1:02d}Z_em'.format(
                md, TT)
            for em in ensemble_members:
                outfile = 'hovs/glm/{0:04d}_{1:02d}Z/ensemble_plots/hovmoller_4p4_{0:04d}_{1:02d}Z_em{2:02d}.png'.format(
                    md, TT, em)
                filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver_4p4_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(
                    md, TT, em)
                track_data = np.load(filepath)
                lats = track_data['lats']
                lons = track_data['lons']

                [y0, x0] = [lats, lons]  # centre of storm
                df = fpath + '{0:02d}.pp'.format(em)
                v_azi_list = load_ens_members(em, fpath, x0, y0)
                plot_hovmoller(v_azi_list, outfile)


if __name__ == '__main__':
    main()
