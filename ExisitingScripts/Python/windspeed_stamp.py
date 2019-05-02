'''
Produces a simple stamp plot of windspeed for an ensemble of size N.
John Ashcroft, February 2019.
'''

import numpy as np
import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib
from matplotlib.colors import from_levels_and_colors
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.cm as mpl_cm
# University System python may be broken
# If some one insists on using it...
backend = mpl.get_backend()
if backend == 'Qt4Agg' and sys.version_info[0] == 2:
    # Fix the backend
    print('swapping to Agg Backend')
    print('Please consider using anaconda')
    mpl.use('Agg')
# DO NOT MOVE ABOVE BACKEND FIX
import matplotlib.pyplot as plt  # KEEP ME HERE!!!
################################


def main():
    # Define all constraints and locate data.
    data_loc = '/nfs/a37/scjea/model_runs/Hagupit/u-ap087/data/4p4/1203_12Z/wind_plev/'
    # Should include everything but 'NN.pp' where NN is the ensemble member.
    data_name = 'uvplev_4p4_1203_12Z_em'
    outfile_loc = './windspeed_plots/'
    # time data to be added to the outfile
    outfile_name = 'wspeed_1203_12Z_4p4_850hPa_'

    ens_members = np.arange(12)  # Ensemble members to use (i.e. 12)
    plev = 850  # Must be alevel that is in the data.

    md = 1203
    TT = 12  # Forecast initial times
    model = '4p4'

    v_times = [0, 3, 6, 9, 12, 15, 18, 21]  # Times during day to plot
    init_day = 3
    init_time = 12  # The first plot will be of the time after this time, i.e. in this case the 4th at 15Z
    final_day = 8
    final_time = 9  # Last time to plot
    v_days = np.arange(init_day, init_day + 6, 1)
    yr = 2014
    mth = 12

    domain_rad = 4.  # Domain is 2*domain_rad x 2*domain_rad

    u_stash = 'm01s15i201'
    v_stash = 'm01s15i202'
    slp_stash = 'm01s16i222'
    u_constraint = iris.AttributeConstraint(STASH=u_stash)
    v_constraint = iris.AttributeConstraint(STASH=v_stash)
    slp_constraint = iris.AttributeConstraint(STASH=slp_stash)

    # Determine the dimensions of the stamp plot according to the number of members
    n_ems = len(ens_members)
    nrows, ncols = find_subplot_dims(n_ems)

    p_constraint = iris.Constraint(pressure=plev)

    for dd in v_days:
        for hr in v_times:
            if dd < init_day or (dd == init_day and hr < init_time + 1):
                continue
            if dd > final_day or (dd > final_day - 1 and hr > final_time - 1):
                continue
            lead_time = (dd - init_day) * 24 + (hr - init_time)
            print 'Plotting data for dd = {0}, hr = {1}'.format(dd, hr)
            lead_time = (dd - init_day) * 24 + (hr - init_time)
            time_constraint = iris.Constraint(
                time=iris.time.PartialDateTime(year=yr, month=mth, day=dd, hour=hr))
            outfile = outfile_loc + outfile_name + \
                '{0:02d}_{1:02d}Z.png'.format(dd, hr)

            # Fix domain according to position of ensemble member 0, this assumes the position is similar in ensemble members
            [cenlat, cenlon] = find_centre_3h(md, TT, 0, dd, hr, model)
            minlat = cenlat - domain_rad
            maxlat = cenlat + domain_rad
            minlon = cenlon - domain_rad
            maxlon = cenlon + domain_rad
            b_constraint = box_constraint(minlat, maxlat, minlon, maxlon)

            # Create figure
            fig, axs = plt.subplots(nrows, ncols, dpi=100, subplot_kw={
                                    'projection': ccrs.PlateCarree()})
            for i, em in enumerate(ens_members):

                # Determine figure coordinate
                ic = i % ncols
                ir = i / ncols
                if len(axs.shape) > 1:
                    ax = axs[ir, ic]
                else:
                    ax = axs[ic]

                # Load the data for this ensemble member at this time
                df = data_loc + data_name + '{0:02d}.pp'.format(em)
                with iris.FUTURE.context(cell_datetime_objects=True):
                    u = iris.load_cube(df, u_constraint).extract(
                        time_constraint & p_constraint & b_constraint)
                    v = iris.load_cube(df, v_constraint).extract(
                        time_constraint & p_constraint & b_constraint)
                ws = (u**2 + v**2)**0.5

                coord = iris.coords.AuxCoord(em, long_name='ensemble_member')
                u.add_aux_coord(coord)
                v.add_aux_coord(coord)
                ws.add_aux_coord(coord)

                wspeed_plt = plot_wspeed(ax, ws)
                plot_winds(ax, u, v, model)

                # Make sure only the outside axis are labelled
                if i % ncols != 0:
                    left_label = False
                else:
                    left_label = True
                if i < ncols * (nrows - 1):
                    bottom_label = False
                else:
                    bottom_label = True

                # Add grid lines and ticks etc.
                map_formatter(ax, bottom_label=bottom_label, left_label=left_label,
                              labelsize=20, tick_base_x=2, tick_base_y=2)

                # Annotate to add the ensemble member onto the plot
                ax.annotate('{0:02d}'.format(em), xy=(0.97, 0.03), xycoords='axes fraction', horizontalalignment='right',
                            verticalalignment='bottom', color='k', backgroundcolor='white', fontsize=20)

            # Reduce white space and then make whole figure bigger, keeping the aspect ratio constant.
            plt.gcf().subplots_adjust(hspace=0.025, wspace=0.025,
                                      bottom=0.05, top=0.95, left=0.075, right=0.925)
            w, h = plt.gcf().get_size_inches()
            plt.gcf().set_size_inches(w * 3, h * 3)
            ax = fig.get_axes()

            # Format colourbar
            cbar = plt.colorbar(wspeed_plt, ax=fig.get_axes(
            ), orientation='horizontal', extend='both', fraction=0.046, pad=0.09)
            cbar.ax.tick_params(labelsize=18)
            cbar.set_label('ms$^{-1}$', size=18)
            axs = fig.gca()

            # Add initial time, valid time etc.
            axs.annotate('Initial time: {0}/{1}/{2}, {3:02d}Z'.format(str(md)[-2:], str(md)[:2], yr, TT), xy=(
                0.95, 0.95), xycoords='figure fraction', horizontalalignment='right', verticalalignment='top', color='k', fontsize=15)
            axs.annotate('Valid time: {0:02d}/{1:02d}/{2}, {3:02d}Z'.format(dd, mth, yr, hr), xy=(0.95, 0.925),
                         xycoords='figure fraction', horizontalalignment='right', verticalalignment='top', color='k', fontsize=15)
            axs.annotate('T+{0}Z'.format(lead_time), xy=(0.95, 0.9), xycoords='figure fraction',
                         horizontalalignment='right', verticalalignment='top', color='k', fontsize=15)

            # plt.show()
            plt.savefig(outfile)
            plt.close()


def find_subplot_dims(n):
    n_cols = np.ceil(n**0.5).astype(int)
    n_rows = np.ceil(1. * n / n_cols).astype(int)
    return n_rows, n_cols


def plot_wspeed(ax, ws):
    levels = np.arange(11) * 5
    cmap = mpl_cm.get_cmap('plasma')
    x = ws.coord('longitude').points
    y = ws.coord('latitude').points
    X, Y = np.meshgrid(x, y)

    wspeed = ax.contourf(X, Y, ws.data, cmap=cmap,
                         levels=levels, extend='both')
    return wspeed


def plot_winds(ax, u, v, mod):
    if mod == '4p4':
        scale = 300.
        step = 12
    else:
        scale = 400.
        step = 8  # Need to trial and error these.

    x = u.coord('longitude').points
    y = u.coord('latitude').points
    X, Y = np.meshgrid(x, y)

    arrows = ax.quiver(X[::step, ::step], Y[::step, ::step], u.data[::step, ::step],
                       v.data[::step, ::step], angles='xy', scale=scale, transform=ccrs.PlateCarree())
    if u.coord('ensemble_member').points[0] == 0:
        ax.quiverkey(arrows, 0.85, 0.1, 20,
                     '20ms$^{-1}$', coordinates='figure', fontproperties={'size': 15})


def map_formatter(ax, tick_base_x=15.0, tick_base_y=15.0, labelsize=20, top_label=False, bottom_label=True, right_label=False, left_label=True, central_longitude=0.0, res='\
10m'):
    '''
    Adds gridlines, countries and lat/lon labels to the plot.
    '''
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.75, color='k', linestyle=':')
    gl.xlabels_top = top_label
    gl.ylabels_right = right_label
    gl.ylabels_left = left_label
    gl.xlabels_bottom = bottom_label
    gl.xlocator = mticker.MultipleLocator(base=tick_base_x)
    gl.ylocator = mticker.MultipleLocator(base=tick_base_y)
    ax.coastlines(resolution=res, color='k', linewidth=1)
    gl.xlabel_style = {'size': labelsize}
    gl.ylabel_style = {'size': labelsize}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER


def find_centre_3h(md, TT, em, v_day, v_time, mod):
    '''
    Loads the track data for a specific time. This function is tailored to my own data, will have to be rewritten depending on how the users track data is stored.
    If you want just a single domain simply create a function that returns the centre of this domain.
    '''
    if mod == '4p4':
        filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver3hr_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(
            md, TT, em, mod)
    elif mod == 'glm':
        filepath = "/nfs/a37/scjea/plot_scripts/tc_analysis/intensity_track/track_data/tcver3hr_newinterp_{3}_{0:04d}_{1:02d}Z_en{2:02d}.npz".format(
            md, TT, em, mod)
    track_data = np.load(filepath)
    lats = track_data['lats']
    lons = track_data['lons']
    days = track_data['vt_days']
    times = track_data['vt_times']
    index1 = np.where(days == v_day)
    lats = lats[index1]
    lons = lons[index1]
    times = times[index1]
    index2 = np.where(times == v_time)
    cenlat = lats[index2]
    cenlon = lons[index2]
    return cenlat[0], cenlon[0]


def box_constraint(minlat, maxlat, minlon, maxlon):
    # Create constraint to extract data from cube over a certain region
    longitude_constraint1 = iris.Constraint(
        longitude=lambda cell: cell > minlon)
    longitude_constraint2 = iris.Constraint(
        longitude=lambda cell: cell < maxlon)
    latitude_constraint1 = iris.Constraint(latitude=lambda cell: cell > minlat)
    latitude_constraint2 = iris.Constraint(latitude=lambda cell: cell < maxlat)
    box_constraint = longitude_constraint1 & longitude_constraint2 & latitude_constraint1 & latitude_constraint2
    return box_constraint


if __name__ == '__main__':
    main()
