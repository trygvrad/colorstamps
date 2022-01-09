import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from .stamps import get_cmap, parse_name_postfix, colorspace
import colorspacious



def apply_stamp(data_0, data_1, cmap, vmin=None, vmax = None, vmin_0 = None, vmax_0 = None, vmin_1 = None, vmax_1 = None, clip = None,
        l = None, rot = None, J = None, sat = None, limit_sat = None, a= None, b = None):
    '''
        Applies a cmap to data_0 and data_1 which formas axis 0 and 1 on the cmap
        if vmin/vmax is not specified the limits are set by set to np.min() and np.max() for data_0 and data_1

        Args:
            data_0: np array of shape(n,m) to apply to cmap axis 0
            data_1: np array of shape(n,m) to apply to cmap axis 1
            cmap: np array of shapek (k,k,3) with rgb values or a string descrbing a named colormap
            vmin[max]: optional float, lower[upper] bound for both data_0 and data_1
            vmin[max]_0[1]: optional float, lower[upper] bound for data_0[1], overrides vmin[max]
            clip: (string),
                'circle'  clips the data to fit a circular colormap using radial projection,
                'square', clips the data to fit a square colormap using radial projection
                'box',    clips the data to fit a square colormap using by limiting the first and second axis independently
                'none'    data out of bounds will have r,g,b = np.nan
                default behaviour is 'circle' or 'square' square if a radial-type cmap is called by name, otherwise 'box'

            l, rot, J, sat, limit_sat, a, b: arguments passed to stamps.get_cmap() when a colormap is called by name, see stamps.get_cmap() for explanation
            a and b: if clip = 'circle' these also modify the clip

        Returns:
            tuple (rgb, stamp)
            where rgb is a numpy array of (n,m,3) with values 0 -> 1 defining sRGB color in image
            and stamp is a colorstamps.helpers.Stamp object.

    '''

    if type(cmap) == type(''):
        cmap, a, b = parse_name_postfix(cmap, a, b) # if there is a _ in the name, set a and b according to the postfix, and remove the postfix

    if clip is None:
        if type(cmap) == type('') and cmap in ['disk','cone','funnel']:
            clip = 'circle'
        elif type(cmap) == type('') and cmap in ['flat','peak','abyss','hsv']:
            clip = 'square'
        else:
            clip = 'box'
    stamp = Stamp(cmap, l=l, rot = rot, J = J, sat = sat, limit_sat = limit_sat, a = a, b = b)
    cmap = stamp.cmap
    if vmin_0 is None:
        vmin_0 = vmin
    if vmin_1 is None:
        vmin_1 = vmin
    if vmax_0 is None:
        vmax_0 = vmax
    if vmax_1 is None:
        vmax_1 = vmax

    if vmin_0 is None:
        stamp.vmin_0 = np.min(data_0)
    else:
        stamp.vmin_0 = vmin_0
    if vmax_0 is None:
        stamp.vmax_0 = np.max(data_0)
    else:
        stamp.vmax_0 = vmax_0

    if vmin_1 is None:
        stamp.vmin_1 = np.min(data_1)
    else:
        stamp.vmin_1 = vmin_1
    if vmax_1 is None:
        stamp.vmax_1 = np.max(data_1)
    else:
        stamp.vmax_1 = vmax_1

    l = cmap.shape[0]
    if clip == 'box':
        colorax_0 = np.array(l*(data_0-stamp.vmin_0)/(stamp.vmax_0-stamp.vmin_0), dtype = int)
        colorax_1 = np.array(l*(data_1-stamp.vmin_1)/(stamp.vmax_1-stamp.vmin_1), dtype = int)

        colorax_0[colorax_0<0] = 0
        colorax_1[colorax_1<0] = 0
        colorax_0[colorax_0>l-1] = l-1
        colorax_1[colorax_1>l-1] = l-1

        rgb = stamp.cmap[colorax_0,colorax_1]

    if clip == 'square':

        # get a and b values of data (y and x coordinates in the cmap)
        # also calculate the phi and rin the cmap space
        a = 2*(data_0-stamp.vmin_0)/(stamp.vmax_0-stamp.vmin_0)-1 # -1 -> 1 for vmin to vmax
        b = 2*(data_1-stamp.vmin_1)/(stamp.vmax_1-stamp.vmin_1)-1 # -1 -> 1 for vmin to vmax
        phi = np.arctan2(a,b)
        r = np.sqrt(a**2+b**2) + 3/l

        colorax_0 = np.array(l*0.5*(a+1), dtype = int)
        colorax_1 = np.array(l*0.5*(b+1), dtype = int)


        mask = np.zeros(colorax_0.shape, dtype = bool)
        mask[colorax_0<0] = 1
        mask[colorax_1<0] = 1
        mask[colorax_0>l-1] = 1
        mask[colorax_1>l-1] = 1

        colorax_0[colorax_0<0] = 0
        colorax_1[colorax_1<0] = 0
        colorax_0[colorax_0>l-1] = l-1
        colorax_1[colorax_1>l-1] = l-1

        rgb = stamp.cmap[colorax_0,colorax_1]

        # make a mask that describes the square at the edge of valid colors:

        cmap_mask = np.ones((l,l))
        cmap_mask[1:-1,1:-1] = 0

        lspace = np.linspace(-1,1,l)
        cmap_atan = np.arctan2(lspace[:,np.newaxis]*np.ones((l,l)),lspace[np.newaxis,:]*np.ones((l,l)))

        # get the coordinates and angle of the edge of the cmap mask sorted according to the angle
        where = np.argwhere(cmap_mask)
        a_vec = where[:,0]
        b_vec = where[:,1]
        angle = cmap_atan[a_vec, b_vec]
        order = angle.argsort()
        angle = angle[order]
        a_vec = a_vec[order]
        b_vec = b_vec[order]

        # for each pixel that does not fit in the area, find the pixel in [a_vec,b_vec] with the closest angle
        phis_to_fit = phi[mask]
        idx = np.searchsorted(angle, phis_to_fit, side="left")
        idx += 1
        angle = np.concatenate(([angle[-1]-2*np.pi],angle,[2*np.pi+angle[0]]))
        idx[ np.abs(phis_to_fit - angle[idx-1]) < np.abs(phis_to_fit - angle[idx])] -=1
        idx -= 1
        loc_a = a_vec[idx]
        loc_b = b_vec[idx]

        # for any rgb not correctly set, update rgb
        rgb[mask] = stamp.cmap[loc_a,loc_b]

    elif clip == 'circle':

        if a is None:
            a = (-1,1)
        if b is None:
            b = (-1,1)
        # get a and b values of data (y and x coordinates in the cmap)
        # also calculate the phi and r in the cmap space
        a_img = (data_0-stamp.vmin_0)/(stamp.vmax_0-stamp.vmin_0) # 0 -> 1 for vmin to vmax
        b_img = (data_1-stamp.vmin_1)/(stamp.vmax_1-stamp.vmin_1) # 0 -> 1 for vmin to vmax
        a_img = a_img*(a[1]-a[0]) + a[0]  # a[0] -> a[1] for vmin to vmax
        b_img = b_img*(b[1]-b[0]) + b[0]  # b[0] -> b[1] for vmin to vmax

        phi = np.arctan2(a_img, b_img)
        r = np.sqrt(a_img**2 + b_img**2)
        mask = r > 1
        colorax_0 = np.array(l/(a[1]-a[0])*(a_img-a[0]), dtype = int)
        colorax_1 = np.array(l/(b[1]-b[0])*(b_img-b[0]), dtype = int)

        colorax_0[colorax_0 < 0] = 0
        colorax_1[colorax_1 < 0] = 0
        colorax_0[colorax_0 > l-1] = l-1
        colorax_1[colorax_1 > l-1] = l-1

        rgb = stamp.cmap[colorax_0, colorax_1]


        mask += np.isnan(rgb[:,:,0])

        # make a mask that describes the circle at the edge of valid colors:
        '''
        ar = np.linspace(a[0],a[1],l)
        br = np.linspace(b[0],b[1],l)
        cmap_ab = np.sqrt(ar[:,np.newaxis]**2+br[np.newaxis,:]**2)
        cmap_atan = np.arctan2(ar[:,np.newaxis]*np.ones((l,l)),br[np.newaxis,:]*np.ones((l,l)))
        corr = np.sin(np.abs(2*(cmap_atan%(0.5*np.pi)-0.25*np.pi)))
        cmap_mask = (cmap_ab > 1.0-1.45/l-0.58/l*corr )& (cmap_ab < 1)
        # include edge if edge less than 1
        edge_mask = cmap_ab<1.0
        edge_mask[1:-1,1:-1] = 0
        cmap_mask += edge_mask'''

        ar = np.linspace(a[0], a[1],l)
        br = np.linspace(b[0], b[1],l)
        cmap_ab = np.sqrt(ar[:,np.newaxis]**2 + br[np.newaxis,:]**2)
        cmap_atan = np.arctan2(ar[:,np.newaxis]*np.ones((l,l)), br[np.newaxis,:]*np.ones((l,l)))

        edge_mask = np.array(cmap_ab < 1.0, dtype = int)

        cmap_mask = np.zeros(edge_mask.shape, dtype = int)
        cmap_mask[:-1,:] +=  edge_mask[:-1,:] > edge_mask[1:,:]
        cmap_mask[ 1:,:] +=   edge_mask[1:,:] > edge_mask[:-1,:]
        cmap_mask[:,:-1] +=  edge_mask[:,:-1] > edge_mask[:,1:]
        cmap_mask[:, 1:] +=   edge_mask[:,1:] > edge_mask[:,:-1]

        cmap_mask[:,0] += edge_mask[:,0]
        cmap_mask[0,:] += edge_mask[0,:]
        cmap_mask[:,-1] += edge_mask[:,-1]
        cmap_mask[-1,:] += edge_mask[-1,:]

        cmap_mask = cmap_mask>0

        # get the coordinates and angle of the edge of the cmap mask sorted according to the angle
        where = np.argwhere(cmap_mask)
        a_vec = where[:,0]
        b_vec = where[:,1]
        angle = cmap_atan[a_vec, b_vec]
        order = angle.argsort()
        angle = angle[order]
        a_vec = a_vec[order]
        b_vec = b_vec[order]

        # for each pixel that does not fit in the area, find the pixel in [a_vec,b_vec] with the closest angle
        phis_to_fit = phi[mask]
        idx = np.searchsorted(angle, phis_to_fit, side="left")
        idx += 1
        angle = np.concatenate(([angle[-1]-2*np.pi],angle,[2*np.pi+angle[0]]))
        idx[ np.abs(phis_to_fit - angle[idx-1]) < np.abs(phis_to_fit - angle[idx])] -=1
        idx -= 1
        loc_a = a_vec[idx]
        loc_b = b_vec[idx]

        # for any rgb not correctly set, update rgb
        rgb[mask] = stamp.cmap[loc_a,loc_b]

    elif clip == 'none':
        # 0.9999999 to approximate inclusive vmax
        colorax_0 = np.array(0.999999999*l*(data_0-stamp.vmin_0)/(stamp.vmax_0-stamp.vmin_0), dtype = int)
        colorax_1 = np.array(0.999999999*l*(data_1-stamp.vmin_1)/(stamp.vmax_1-stamp.vmin_1), dtype = int)

        none_mask = np.zeros(colorax_0.shape, dtype = int)
        none_mask[colorax_0<0] += 1
        none_mask[colorax_1<0] += 1
        none_mask[colorax_0>l-1] += 1
        none_mask[colorax_1>l-1] += 1

        colorax_0[none_mask>0] = 0
        colorax_1[none_mask>0] = 0

        rgb = stamp.cmap[colorax_0,colorax_1]
        rgb[:,:,0][none_mask>0] = None
        rgb[:,:,1][none_mask>0] = None
        rgb[:,:,2][none_mask>0] = None

    return rgb, stamp


class Stamp:
    '''
    Class for holding 2d colormaps

    In normal use class objects of this class will be generated by apply_stamp(), and not generated.
    May be called directly to inspect colormaps using the member funciton Stamp.eval() without having any data, as needed for apply_stamp().
    Member variables self.vmin_0, self.vmin_1, self.vmax_0, self.vmax_1 are set when invoked by apply_stamp() and determine the limts.

    Args:
        cmap:
            numpy array of shape (l,l,3) with rgb values
            or string for a valid colormap available with stamps.get_cmap()

        l, rot, J, sat, limit_sat, a, b: passed to stamps.get_cmap() if cmap is a string, otherwise ignored

    '''

    def __init__(self, cmap, l = None, rot = None, J = None, sat = None, limit_sat = None, a = None, b = None):
        if type(cmap) == type(''):
            self.name = cmap
            cmap = get_cmap(cmap, l=l, rot = rot, J = J, sat = sat, limit_sat = limit_sat, a = a, b = b)
        self.cmap = cmap
        self.l = cmap.shape[0]
        self.vmin_0 = 0
        self.vmin_1 = 0
        self.vmax_0 = 1
        self.vmax_1 = 1

    def overlay_ax(self, ax, lower_left_corner = [0.66, 0.66], width = None, height = 0.25):
        '''
        overlays new axes over the existing axes ('ax') with the colormap

        Args:
            ax: axes to position new axes over
            lower_left_corner: [y,x] coordinates in fraction coorniates of ax descrbing position of new ax
            width: float (optional), width of new axes
            heght: float (optional) height of the new axes (if width not specified)
        returns:
            cmap_ax, matplotlib axes object with the colormap
        '''

        ax_pos = ax.get_position() # get the original position
        fig = ax.figure
        fig_width, fig_height = fig.get_size_inches()
        ax_height_in_inches = fig_height*ax_pos.height
        ax_width_in_inches = fig_width*ax_pos.width

        if not width is None:
            new_ax_width_in_inches = width*ax_width_in_inches
            new_ax_height_in_inches = new_ax_width_in_inches
        else:
            new_ax_height_in_inches = height*ax_height_in_inches
            new_ax_width_in_inches = new_ax_height_in_inches

        new_ax_width = new_ax_width_in_inches/fig_width
        new_ax_height = new_ax_height_in_inches/fig_height

        new_pos = [ax_pos.x0 + ax_pos.width*lower_left_corner[1], ax_pos.y0 + ax_pos.height*lower_left_corner[0],
                            new_ax_width, new_ax_height]
        cmap_ax = fig.add_axes(new_pos)
        extent = (self.vmin_1, self.vmax_1, self.vmin_0, self.vmax_0)
        cmap_ax.imshow(self.cmap, origin = 'lower', extent = extent, aspect='auto')
        return cmap_ax

    def show_in_ax(self, cmap_ax):
        '''
        plots the colormap as an image in ax with correct limits

        Args:
            cmap_ax: matplotlib axes object
        returns:
            None
        '''

        extent = (self.vmin_1, self.vmax_1, self.vmin_0, self.vmax_0)
        dxdy = (self.vmax_1-self.vmin_1)/(self.vmax_0-self.vmin_0)
        cmap_ax.imshow(self.cmap, origin = 'lower', extent = extent, aspect=dxdy)
        return


    def eval(self, axes = None):
        '''
        plots the colormap overlaid with contour plots of J, hue and sat and simulated partial colorblindness
        (partial deuteranomaly appears to be the most common kind https://en.wikipedia.org/wiki/Color_blindness)

        Args:
            axes: default None, alternatively a list of 5 matplotlib axes
                if None, a fig with 5 axes will be created
        returns:
            fig, axes if axes are not declared as input, otherwise None
        '''

        if axes is None:
            fig, axes = plt.subplots(1,5, figsize = (10,2))
        axes = axes.ravel()
        for ax in axes:
                ax.set_xticks([])
                ax.set_yticks([])

        axes[0].imshow(self.cmap, origin = 'lower')
        if hasattr(self, 'name'):
            axes[0].set_title(self.name)


        rgb = np.nan_to_num(self.cmap)
        Jab = colorspacious.cspace_convert(rgb, "sRGB1", colorspace)
        phi = np.arctan2(Jab[:,:,1],Jab[:,:,2])
        sat = np.sqrt(Jab[:,:,1]**2 + Jab[:,:,2]**2)


        for ax in axes[1:4]:
            ax.imshow(self.cmap, origin = 'lower', alpha = 0.5)
        lw = 0.5
        axes[1].contour(Jab[:,:,0], levels = 21, colors = [[0,0,0]], linewidths = lw)
        axes[1].set_title('J')
        #axes[4].contour(Jab[:,:,1], levels = 21, colors = [[0,0,0]], linewidths = lw)
        #axes[4].set_title('a')
        #axes[5].contour(Jab[:,:,2], levels = 21, colors = [[0,0,0]], linewidths = lw)
        #axes[5].set_title('b')

        axes[2].contour(phi, levels = 21, colors = [[0,0,0]], linewidths = lw)
        axes[2].set_title('hue')
        axes[3].contour(sat, levels = 21, colors = [[0,0,0]], linewidths = lw)
        axes[3].set_title('saturation')

        cvd_space = {"name": "sRGB1+CVD",
                    "cvd_type": "deuteranomaly",
                    "severity": 50}
        deuteranomaly = colorspacious.cspace_convert(self.cmap, cvd_space, "sRGB1")
        deuteranomaly[deuteranomaly<0] = 0
        deuteranomaly[deuteranomaly>1] = 1

        axes[4].imshow(deuteranomaly, origin = 'lower')
        axes[4].set_title('50% Deuteranomaly')

        if 'fig' in locals():
            fig.patch.set_facecolor('white')
            return fig, axes



def get_random_data(seed = 0):
    '''
    returns some random(seeded) data for examples
    '''

    s = 100
    t = 200
    img = np.zeros((s,t,2))
    x_r = np.arange(s)
    y_r = np.arange(t)

    np.random.seed(seed)
    for i in range(100):
        x = int(np.random.random()*s)
        y = int(np.random.random()*t)
        z = 1*np.random.random()
        sigma = 0.02*(np.random.random()+0.5)
        ang = 2*np.pi*np.random.random()

        img[:,:,0] += z*np.sin(ang)*np.exp(-sigma*((x-x_r[:,np.newaxis])**2+(y-y_r[np.newaxis,:])**2))
        img[:,:,1] += z*np.cos(ang)*np.exp(-sigma*((x-x_r[:,np.newaxis])**2+(y-y_r[np.newaxis,:])**2))
    return img
