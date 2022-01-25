import colorspacious
import numpy as np

colorspace = 'CAM02-LCD'

def mask_rgb(rgb, a, b, mask):
    '''
    function that masks an rgb colormap with np.nan according to the string mask

    Args:
        rgb: (l,l,3) matrix
        a,b: values of a and b. if mask = 'circle' anyting with sqrt(a**2+b**2)>1 will be np.nan
        mask: string:
            'circle' -> masks everything outside a circle defined as where sqrt(a**2+b**2)>1
            'no-mask' -> do nothing
            'unavailable'-> masks invalid rgb values (i.e. <0 or >1)
    '''

    if mask == 'unavailable':
        rgb[rgb[:,:,:]<0] = np.nan
        rgb[rgb[:,:,:]>1] = np.nan
        mask = np.isnan(np.sum(rgb[:,:,:], axis = -1))
        rgb[:,:,0][mask] = np.nan
        rgb[:,:,1][mask] = np.nan
        rgb[:,:,2][mask] = np.nan
    elif mask == 'no_mask':
        None
    elif mask == 'circle':
        l = rgb.shape[1]
        a_1 = np.linspace(a[0],a[1],l)
        b_1 = np.linspace(b[0],b[1],l)
        ab = np.sqrt(a_1[:,np.newaxis]**2+b_1[np.newaxis,:]**2)
        mask = ab > 1
        rgb[:,:,0][mask] = np.nan
        rgb[:,:,1][mask] = np.nan
        rgb[:,:,2][mask] = np.nan
    else:
        raise ValueError("mask must be 'no_mask', 'unavailable' or 'circle'")

def set_ab_rot(Jab, ar, br, rot):
    '''
    sets the [:,:,1] and [:,:,2] axes of a Jab colormap to ar and br
    then rotates the ab color plane according to the angle rot

    Args:
        Jab: (l,l,3) colormap
        ar: 1d array, typically made by np.linspace()
        br: 1d array, typically made by np.linspace()
        rot: angle in degrees
    returns:
        None (but Jab changed in-place)
    '''

    if rot==0:
        Jab[:,:,1] = ar[:,np.newaxis]
        Jab[:,:,2] = br[np.newaxis,:]
    else:
        ab = np.sqrt(ar[:,np.newaxis]**2+br[np.newaxis,:]**2)
        Jab[:,:,1] = ar[:,np.newaxis]
        Jab[:,:,2] = br[np.newaxis,:]
        phi = np.arctan2(Jab[:,:,1],Jab[:,:,2])+rot*np.pi/180
        Jab[:,:,2] = ab*np.cos(phi)
        Jab[:,:,1] = ab*np.sin(phi)

def get_const_J(J = 95, a = (-1,1), b = (-1,1), r = 33.0, l=256, mask = 'no_mask', rot = 0):
    '''
    Generates am rgb colormap of (l,l,3) that attempts to keep a constant lightness in the CAM02-LCD colorspace
    The colormap is based on the a-b plane of the Jab colorspace for a constant J.

    Args:
        J: float (lighness), default 95, range approximately 1->128,
        a: tuple of 2 floats, default (-1,1). The limit along the a-axis will be (a[0]*r,a[1]*r)
        b: tuple of 2 floats, default (-1,1). The limit along the b-axis will be (b[0]*r,b[1]*r)
        r: float, default 33.0. The saturation where a or b is 1. (named 'r' for radius in the a-b plane)
        l: int, default 256. Size of the colormap.
        mask: string, default 'no_mask'.
                If 'circle' makes a circular mask, and everything outside will be np.nan
                If 'unavailable' makes a colors that "should" have rgb<0 or rgb>1 when transformed to sRGB will be np.nan
        rot: rotation of the hues on the a-b plane, in degrees

    returns:
        a (l,l,3) numpy array of rgb values
    '''


    Jab = np.zeros((l,l,3))
    Jab[:,:,0] = J
    ar = np.linspace(r*a[0], r*a[1],l)
    br = np.linspace(r*b[0], r*b[1],l)
    set_ab_rot(Jab, ar, br, rot)

    rgb = colorspacious.cspace_convert(Jab, colorspace, "sRGB1")

    mask_rgb(rgb, a, b, mask)

    rgb[rgb[:,:,:]<0] = 0
    rgb[rgb[:,:,:]>1] = 1

    return rgb

def get_var_J(J = [95,128.5], a = (-1,1), b = (-1,1), r = 33.0, l=256, mask = 'no_mask', rot = 0, limit_sat = None):
    '''
    Generates am rgb colormap of (l,l,3) that attempts to keep a constant lightness in the CAM02-LCD colorspace
    The colormap is based on the a-b plane of the Jab colorspace for a constant J.

    Args:
        J: (lighness) tuple of 2 floats, default [95,128.5] defining the range of lightness for the colormap, default 95,
            max range of J approximately 1 to 128.5
        a: tuple of 2 floats, default (-1,1). The limit along the a-axis will be (a[0]*r,a[1]*r)
        b: tuple of 2 floats, default (-1,1). The limit along the b-axis will be (b[0]*r,b[1]*r)
        r: float, default 33.0. The saturation where a or b is 1. (named 'r' for radius in the a-b plane)
        l: int, default 256. Size of the colormap.
        mask: string, default 'no_mask'.
                If 'circle' makes a circular mask, and everything outside will be np.nan
                If 'unavailable' makes a colors that "should" have rgb<0 or rgb>1 when transformed to sRGB will be np.nan
        rot: rotation of the hues on the a-b plane, in degrees

    returns:
        a (l,l,3) numpy array of rgb values
    '''

    Jab = np.zeros((l,l,3))
    ar = np.linspace(r*a[0], r*a[1],l)
    br = np.linspace(r*b[0], r*b[1],l)
    set_ab_rot(Jab, ar, br, rot)
    ab = np.sqrt(ar[:,np.newaxis]**2+br[np.newaxis,:]**2)

    a_1 = np.linspace(a[0], a[1],l)
    b_1 = np.linspace(b[0], b[1],l)

    Jab[:,:,0] = J[0] + (J[1]-J[0])*(1-np.sqrt(a_1[:,np.newaxis]**2+b_1[np.newaxis,:]**2))
    Jab[Jab[:,:,0]<1,0] = 1
    Jab[Jab[:,:,0]>128,0] = 128
    if not (limit_sat is None):
        apply_radial_sat_limit(Jab, limit_sat = limit_sat)
    rgb = colorspacious.cspace_convert(Jab, colorspace, "sRGB1")

    mask_rgb(rgb, a, b, mask)


    rgb[rgb[:,:,:]<0] = 0
    rgb[rgb[:,:,:]>1] = 1

    #print(r)
    return rgb



def parse_name_postfix(cmap, a, b):
    '''
    if a cmap name has a postfix that details the quadrant/side, this will translate that to ranges in a and/or b.
    example: parse_name_postfix('cone tr', a, b) return a and b so that they span the top right quadrant
    inputs a and b so that both can be returned even if only one is changed

    Args:
        cmap: string, potentially with a postfix detailing quadrant/side following a space, i.e. 'cone tr'.
            The postfix translates:
            'b' -> bottom,
            't' -> top,
            'l' -> left,
            'r' -> right.
            Any combination (b,t)+(l,r) is possible to select quadrants
        a: current limits for a, not checked but should be a tuple of length 2
        b: current limits for b, not checked but should be a tuple of length 2

    returns:
        tuple (cmap, a, b):
        cmap (stripped of the postfix)
        a, tuple of length 2
        b, tuple of length 2
    '''
    # check if the cmap name has additional info regarding quadrant/side
    if len(cmap.split(' '))>1:
        param = cmap.split(' ')[1]
        if 'b' in param: a = (0,1)
        if 't' in param: a = (-1,0)
        if 'r' in param: b = (0,1)
        if 'l' in param: b = (-1,0)
        cmap = cmap.split(' ')[0]
    return cmap, a, b

def get_cmap(name, l = None, rot = None, J = None, sat = None, limit_sat = None, a= None, b = None):
    '''
    getter function for named colormaps
    the 'alt' colormaps are rotated 45 degrees
    flat colormaps are: ---------------- 'flat', 'disk'
    colormaps with a bright center: ---- 'peak', 'cone'
    colormaps with a dark center: ------ 'abyss', 'funnel'
    alternate with a bright center: ---- 'hsv', 'fourCorners', 'fourEdges', 'teuling0w' to 'teuling3w'
    colormaps with lighness on y axis: - 'barrel', 'cut', 'blues', 'reds', 'greens', 'yellows'
    teuling colormaps: --------- 'teuling0f', 'teuling1f', 'teuling2f', 'teuling3f', 'teuling0w', 'teuling1w', 'teuling2w', 'teuling3w'

    any matplotlib colormap can also be converted to a colormap with lighness on y-axis

    Args:
        name: string
            For radial colormaps the name may have a postfix separated by a space, i.e. 'cone tr'
            the postfix must be some combination of (t,b) and/or (l,r) which defines the quadrant/side of the colormap to include
            t-> top, b-> bottom, r-> right, l-> left, and 'tr'-> top right, etc.
        l: int, the size of the colormap will be (l,l), defaults to 256 if None
        rot: float, rotation of the colormap (where applicable)
        J: array-like of length 2 (float,float), determins min and max luminocity where applicable
        sat: float, maximum saturation where applicable
        limit_sat: string, 'individual' or 'shared'. How saturation is limited for relevant colormaps when colors outside sRGB are required
                    'individual': each combination J, hue in the colormap has an individual limit to saturation
                    'shared': for each J, all hues share a limit, the maximum where all hues can be represented
        a: range along a-axis, array-like [min,max]
                    Used to move the center of the colormap where applicable.
                    Defaults to (-1,1) which is then multiplied internally with sat
        b: range along b-axis, see a.

    returns:
        a (l,l,3) numpy array of rgb values
    '''

    if l is None: l = 256
    if rot is None: rot = 0
    if sat is None: sat = 33.0
    if a is None: a = (-1,1)
    if b is None: b = (-1,1)
    name, a, b = parse_name_postfix(name, a, b) # if there is a _ in the name, set a and b according to the postfix, and remove the postfix
    if name == 'flat':
        if J is None: J = [95]
        return get_const_J( J = J[0], a = a, b = b, r = sat, l = l, rot = rot)
    elif name == 'disk':
        if J is None: J = [95]
        return get_const_J(J = J[0], a = a, b = b, r = sat, l = l, rot = rot, mask = 'circle')
    elif name == 'peak':
        if J is None: J = [95,128.5]
        return get_var_J(J = J, a = a, b = b, r = sat, l = l, rot = rot, limit_sat = limit_sat)
    elif name == 'cone':
        if J is None: J = [95,128.5]
        return get_var_J(J = J, a = a, b = b, r = sat, l = l, rot = rot, mask = 'circle', limit_sat = limit_sat)
    elif name == 'abyss':
        if J is None: J = [95,1]
        return get_var_J(J = J, a = a, b = b, r = sat, l = l, rot = rot, limit_sat = limit_sat)
    elif name == 'funnel':
        if J is None: J = [95,1]
        return get_var_J(J = J, a = a, b = b, r = sat, l = l, rot = rot, mask = 'circle', limit_sat = limit_sat)
    elif name == 'hsv':
        hsv = np.ones((l,l,3))
        ar = np.linspace(a[0],a[1],l)[:,np.newaxis]*np.ones((l,l))
        br = np.linspace(b[0],b[1],l)[np.newaxis,:]*np.ones((l,l))
        phi = np.arctan2(ar,br)+rot*np.pi/180
        hsv[:,:,0] = phi/np.pi*0.5+0.5
        hsv[:,:,1] = np.sqrt(ar**2+br**2)/np.sqrt(2)
        hsv[:,:,2] = 1
        RGB = matplotlib.colors.hsv_to_rgb(hsv)
        return RGB
    elif name == 'fourEdges':
        return four_edges(l=l, a=a, b=b, rot = rot+90)
    elif name == 'fourCorners':
        return four_edges(l=l, a=(-0.85,0.85), b=(-0.85,0.85), rot = 45)
    # these are shared by all that follow
    if J is None: J = [15,120]
    if limit_sat is None: limit_sat = 'shared'
    # rest of colormaps
    if name == 'barrel':
        return barrel(sat = sat, phi = [-180,180], J =J, l = l, limit_sat = limit_sat)
    elif name == 'cut':
        return cut(a = a, sat = sat, rot = rot, J = J, l = l, limit_sat = limit_sat)
    elif name == 'blues':
        return cut(a = [0,1], sat = sat, rot = 180, J = J, l = l, limit_sat = limit_sat)
    elif name == 'reds':
        return cut(a = [0,1], sat = sat, rot = 90, J = J, l = l, limit_sat = limit_sat)
    elif name == 'greens':
        return cut(a = [0,1], sat = sat, rot = -90, J = J, l = l, limit_sat = limit_sat)
    elif name == 'yellows':
        return cut(a = [0,1], sat = sat, rot = 0, J = J, l = l, limit_sat = limit_sat)
    elif name == 'teuling0f':
        return teuling(l = l, a = 0.32, order = [0,1,2])
    elif name == 'teuling1f':
        return teuling(l = l, a = 0.72, order = [1,0,2])
    elif name == 'teuling2f':
        return teuling(l = l, a = 0.32, order = [1,0,2])
    elif name == 'teuling3f':
        return teuling(l = l, a = 0.32, order = [1,0,2], green_multiplier = 0.75)
    elif name == 'teuling0w':
        return teuling(l = l, a = 0.32, order = [0,1,2], white_center = True)
    elif name == 'teuling1w':
        return teuling(l = l, a = 0.72, order = [1,0,2], white_center = True)
    elif name == 'teuling2w':
        return teuling(l = l, a = 0.32, order = [1,0,2], white_center = True)
    elif name == 'teuling3w':
        return teuling(l = l, a = 0.32, order = [1,0,2], green_multiplier = 0.75, white_center = True)
    elif name == 'orangeBlue':
        return bilinear(l)
    elif name == 'greenPurple':
        return bilinear(l, c0 = [0.5,1,0], c1 = [0.5,0,1])
    elif name == 'greenTealBlue':
        return bilinear(l, c0 = [0,1,0], c1 = [0,0,1])
    elif name == 'redPurpleBlue':
        return bilinear(l, c0 = [1,0,0], c1 = [0,0,1])
    elif name in mpl_cmaps:
        return get_2dcmap_from_mpl(name, J = J, l = l, limit_sat = limit_sat)
    else:
        raise ValueError(f'colormap {name} not known')

def get_sat_limts():
    '''
        returns the a 2d matrix of approximate limits to sat (radius in a-b space) in terms of phi and J
    '''

    if not 'limit' in globals():
        global limit, limit_ax_0_J, limit_ax_1_phi

        phi = np.linspace(-np.pi, np.pi, 256+1)
        J = np.linspace(1,130,128)
        sat = np.linspace(0,70,256)

        J_phi_sat = np.empty((len(J),len(phi),len(sat),3))
        J_phi_sat[:,:,:,0] = J[:,np.newaxis,np.newaxis]
        J_phi_sat[:,:,:,1] = phi[np.newaxis,:,np.newaxis]
        J_phi_sat[:,:,:,2] = sat[np.newaxis,np.newaxis,:]

        Jab = np.empty(J_phi_sat.shape)
        Jab[:,:,:,0] = J_phi_sat[:,:,:,0]
        Jab[:,:,:,1] = J_phi_sat[:,:,:,2]*np.sin(J_phi_sat[:,:,:,1])
        Jab[:,:,:,2] = J_phi_sat[:,:,:,2]*np.cos(J_phi_sat[:,:,:,1])
        rgb = colorspacious.cspace_convert(Jab, colorspace, "sRGB1")
        rgb[rgb>1] = np.nan
        rgb[rgb<0] = np.nan

        flat_rgb = np.sum(rgb, axis = -1)
        flat_rgb[:,:,0] = 0

        # there are some strange regsions in the limits-overview because there are 'jumps' as we go through phi
        # therefore limit the derivative in phi
        for i, _ in enumerate(sat[:-1]):
            flat_rgb[:,0,i]  +=  flat_rgb[:,-1,i]
            flat_rgb[:,-1,i] +=  flat_rgb[:,0,i]
            flat_rgb[:,1:,i+1]  +=  flat_rgb[:,:-1,i]
            flat_rgb[:,:-1,i+1] +=  flat_rgb[:,1:,i]

        flat_rgb[:,0,-1]  +=  flat_rgb[:,-1,-1]
        flat_rgb[:,-1,-1] +=  flat_rgb[:,0,-1]

        valid = np.invert(np.isnan(flat_rgb)) + np.linspace(0,0.9,len(sat))[np.newaxis,np.newaxis,:]
        valid_argmax = np.argmax(valid, axis = -1)

        limit = sat[valid_argmax]
        limit_ax_0_J = J
        limit_ax_1_phi = phi
    return limit, limit_ax_0_J, limit_ax_1_phi

import scipy.interpolate
import matplotlib.cm

def apply_sat_limit(Jab, limit_sat = 'shared'):
    '''
    apply a saturation limit to Jab in order to ensure valid saturation when the limit of the RGB colorspace is reached

    Args:
        Jab: np array of shape (n,m,3) encoded in the colorspace
        limit_sat: 'shared' or 'individual'
            if 'shared', all hues share same limit to saturation (the minimum where all saturation values present in the colormap can be represented)
            if 'individual', different hues have different sauration limits
    returns:
        None (Jab is modified in-place)
    '''

    #limit = sat[valid_argmax]
    #limit_ax_0_J = J
    #limit_ax_1_phi = phi
    limit, limit_ax_0_J, limit_ax_1_phi = get_sat_limts()
    inerpolator = scipy.interpolate.RectBivariateSpline(limit_ax_0_J, limit_ax_1_phi, limit)

    phi = np.arctan2(Jab[:,:,1],Jab[:,:,2])
    sat = np.sqrt(Jab[:,:,1]**2 + Jab[:,:,2]**2)

    max_sat = inerpolator( Jab[:,:,0], phi, grid = False)
    if limit_sat == 'shared':
        max_sat[:,:] = np.min(max_sat, axis=1)[:,np.newaxis]
    mask = sat>max_sat
    #sat[mask] = max_sat[mask]
    change = (max_sat[mask]+0.000000001)/(sat[mask]+0.000000001)
    Jab[mask,1] *= change
    Jab[mask,2] *= change

def apply_radial_sat_limit(Jab, limit_sat = 'shared'):
    '''
    apply a radial saturation limit to Jab in order to make the saturation radial when
    the limit of the RGB colorspace is reached
    the behaviour if limit_sat == 'shared' is different from apply_sat_limit()
    in this function all possible hues are always included, but for apply_sat_limit() only present hues are considered

    Args:
        Jab: np array of shape (n,m,3) encoded in the colorspace
        limit_sat: 'shared' or 'individual'
            if 'shared', all hues share same limit to saturation (the minimum where all are present)
            if 'individual', different hues have different sauration limits
    returns:
        None (Jab is modified in-place)
    '''

    limit, limit_ax_0_J, limit_ax_1_phi = get_sat_limts()
    if limit_sat == 'shared':
        limit_shared = np.min(limit, axis=1)
        inerpolator = scipy.interpolate.interp1d(limit_ax_0_J, limit_shared)
        max_sat = inerpolator( Jab[:,:,0])
    else:
        inerpolator = scipy.interpolate.RectBivariateSpline(limit_ax_0_J, limit_ax_1_phi, limit)
        phi = np.arctan2(Jab[:,:,1],Jab[:,:,2])
        max_sat = inerpolator( Jab[:,:,0], phi, grid = False)
    sat = np.sqrt(Jab[:,:,1]**2 + Jab[:,:,2]**2)
    mask = sat>max_sat
    #sat[mask] = max_sat[mask]
    change = (max_sat[mask]+0.000000001)/(sat[mask]+0.000000001)
    Jab[mask,1] *= change
    Jab[mask,2] *= change

def get_2dcmap_from_mpl(string, J = [15,120], l = 256, limit_sat = 'shared'):
    '''
    Generates a 2d colormap from a 1d colormap found in matplotlib

    Args:
        string: name of the matplotlib colormap
        J: limits to lighness on the y-axis, array like of length 2, default [15,120]
        l: desired size (l,l,3) of the colormap
        limit_sat: string, how to limit the saturation to say within the limits of the RGB colorspace
            'shared': all hues share same limits
            'individual': different hues have different limits
    returns:
        a (l,l,3) numpy array of rgb values
    '''

    cmap = matplotlib.cm.get_cmap(string)
    # make 2d cmap in Jab colorspace
    rgb = np.zeros((l,l,3))
    rgb[:,:,:] = cmap(np.linspace(0,1,l))[np.newaxis,:,:3]
    Jab = colorspacious.cspace_convert(rgb, "sRGB1", colorspace)
    J = np.linspace(J[0], J[1], l)
    Jab[:,:,0] = J[:,np.newaxis]
    # Jab now has colors that cannot be represented in rgb
    # limit the 'saturation' defined as radius in a-b space for a given J according to get_max_ab(J):
    apply_sat_limit(Jab, limit_sat = limit_sat)
    # convert the now limited Jab colorspace to rgb
    rgb = colorspacious.cspace_convert(Jab, colorspace,"sRGB1")
    rgb[rgb<0] = 0
    rgb[rgb>1] = 1
    return rgb

mpl_cmaps = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r',
         'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r',
         'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r',
         'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd',
         'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r',
         'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r',
         'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot',
         'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis',
         'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r',
         'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow',
         'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray',
         'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r',
         'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring',
         'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r',
         'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']

def barrel(sat =  33, phi = [-180,180], J = [15,120], l = 256, limit_sat = 'shared'):
    '''
    Generates a 2d colormap that cycles different hues on the x-axis and has lighness on the y-axis

    Args:
        sat: float, default 33. Desired saturation
        phi: range for the hues on the x-axis in degrees, array like of length 2, default [-180,180]
        J: limits to lighness on the y-axis, array like of length 2, default [15,120]
        l: desired size (l,l,3) of the colormap
        limit_sat: string, how to limit the saturation to say within the limits of the RGB colorspace
            'shared': all hues share same limits
            'individual': different hues have different limits
    returns:
        a (l,l,3) numpy array of rgb values
    '''

    Jab = np.empty((l,l,3))
    J = np.linspace(J[0], J[1], l)
    Jab[:,:,0] = J[:,np.newaxis]
    phi = np.array(phi)/180*np.pi
    phi_linspace = np.linspace(phi[0],phi[1],l)
    Jab[:,:,1] = np.sin(phi_linspace[np.newaxis,:])*sat
    Jab[:,:,2] = np.cos(phi_linspace[np.newaxis,:])*sat
    apply_sat_limit(Jab, limit_sat = limit_sat)
    rgb = colorspacious.cspace_convert(Jab, colorspace,"sRGB1")
    rgb[rgb<0] = 0
    rgb[rgb>1] = 1
    return rgb


def cut(a = (-1,1), sat = 33, rot = 0, J = [15,120], l = 256, limit_sat = 'shared'):
    '''
    Generates a 2d colormap that is bilinear with saturation along the x-axis and lighness on the y-axis
    effectivly the cross-section of the Jab colorspace at some angle rot

    Args:
        sat: float, default 33. Desired saturation
        rot: the hue, or rotation in the a-b plane to make the cut
        J: limits to lighness on the y-axis, array like of length 2, default [15,120]
        l: desired size (l,l,3) of the colormap
        limit_sat: string, how to limit the saturation to say within the limits of the RGB colorspace
            'shared': all hues share same limits
            'individual': different hues have different limits
    returns:
        a (l,l,3) numpy array of rgb values
    '''

    rot = rot*np.pi/180
    Jab = np.empty((l,l,3))
    J = np.linspace(J[0], J[1], l)
    Jab[:,:,0] = J[:,np.newaxis]
    ar = sat*np.linspace(a[0],a[1],l)
    Jab[:,:,1] = np.sin(rot)*ar[np.newaxis,:]
    Jab[:,:,2] = np.cos(rot)*ar[np.newaxis,:]
    apply_sat_limit(Jab, limit_sat = limit_sat)
    rgb = colorspacious.cspace_convert(Jab, colorspace,"sRGB1")
    rgb[rgb<0] = 0
    rgb[rgb>1] = 1
    return rgb

def four_edges(a = (-1,1), b = (-1,1), r_exp = -0.5, f0_exp = 1.7, f1_exp = 1.7, yellow_exp = 0.75, l = 256, rot = 0):
    '''
    Generates a 2d colormap with four colors (r,g,b,y) on the edges (or corners if rot = 45)
    based on the sRGB colorspace, with exponentials

    Args:
        a, b: limits the range of the y,x plane to use for the colormap. Both default to (-1,1), changing this will zoom in/out on different parts of the colormap
        r_exp: float, radial exponent to the lighness, defaults -0.5. Increasing this makes a big white spot in the middle
        f0_exp, f1_exp: floats, exponent for the color along the a and b axes, both default 1.7
        yellow_exp: additioal exponent for yellow, defaults 0.75
        l: size of the colormap, defaults to 256
        rot: rotation of the colormap, in degrees
    returns:
        a (l,l,3) numpy array of rgb values
    '''

    twoD_cmap = np.ones((l,l,3))
    fa = np.linspace(a[0],a[1],l)[:,np.newaxis]*twoD_cmap[:,:,0]
    fb = np.linspace(b[0],b[1],l)[np.newaxis,:]*twoD_cmap[:,:,0]
    if rot == 0:
        f0 = fa
        f1 = fb
    else:
        phi = np.arctan2(fa,fb)+rot*np.pi/180
        r = np.sqrt(fa**2+fb**2)
        f0 = r*np.sin(phi)
        f1 = r*np.cos(phi)

    fac = 1.2*((np.abs(f0)+np.abs(f1))/np.sqrt(f0**2+f1**2))**-0.3
    f0 *= fac
    f1 *= fac
    r = np.sqrt(f0**2+f1**2)
    fr = r**r_exp/np.sqrt(2)**r_exp
    f0 *= fr
    f1 *= fr
    # yellow -> decrease blue
    twoD_cmap[f0<0,2] -= 0.5 * ((-f0[f0<0])**f0_exp)**yellow_exp
    # blue -> decrease red, green
    twoD_cmap[f0>0,0] -= 0.5 * (f0[f0>0])**f0_exp
    twoD_cmap[f0>0,1] -= 0.5 * (f0[f0>0])**f0_exp
    # red -> decrease blue, green
    twoD_cmap[f1<0,1] -= 0.5 * (-f1[f1<0])**f1_exp
    twoD_cmap[f1<0,2] -= 0.5 * (-f1[f1<0])**f1_exp
    # green -> decrease red, blue
    twoD_cmap[f1>0,0] -= 0.5 * (f1[f1>0])**f1_exp
    twoD_cmap[f1>0,2] -= 0.5 * (f1[f1>0])**f1_exp

    twoD_cmap[twoD_cmap>1] = 1
    twoD_cmap[twoD_cmap<0] = 0
    return twoD_cmap

def teuling(l = 256, a = 0.32, order = [1,0,2], white_center = False, green_multiplier = 1.0):
    '''
    Generates a 2d colormap based on:

    Teuling, A. J., R. Stöckli, and Sonia I. Seneviratne. "Bivariate colour maps for visualizing climate data." International journal of climatology 31.9 (2011): 1408-1412.

    Args:
        l: size of the colormap, defaults to 256
        a: float between 0 and 1, determines how the second and third color scales along the x and y axis
        order: The order in which colors are applied. Should be an array of lenght three with the values 0, 1 and 2, corresponding to red, green and blue.
            The color at the first index scales equally (0.5) in both x and y, the two others scale according to a and (1-a) along x and y, and y and x
        white_center: bool, default False. If true, the center will be colored white
        green_multiplier: float, default 1.0. The green part component of the colormap multipled by this value.
            Can help reduce the luminocity of the green corner to produec a more homogeneous-looking map
    returns:
        a (l,l,3) numpy array of rgb values
    '''

    rgb = np.zeros((l,l,3))
    lspace = np.linspace(0,1.0,l)

    rgb[:,:,order[0]] += 0.5*lspace[:,np.newaxis]+0.5*lspace[np.newaxis,:]

    rgb[:,:,order[1]] += a*lspace[:,np.newaxis]+(1-a)*lspace[np.newaxis,::-1]

    rgb[:,:,order[2]] += (1-a)*lspace[::-1,np.newaxis]+a*lspace[np.newaxis,:]

    rgb[:,:,1] *= green_multiplier

    if white_center:
        lim = 1/np.sqrt(2)
        lspace_second = np.linspace(-lim,lim,l)**2
        r = np.sqrt(lspace_second[:,np.newaxis]+ lspace_second[np.newaxis,:])
        rgb[:,:,:] += 0.5*(1-r)[:,:,np.newaxis]
    rgb[rgb<0] = 0
    rgb[rgb>1] = 1
    return rgb


def bilinear(l = 256, c0 = [1,0.5,0], c1 = [0,0.5,1]):
    '''
    Returns an l by l colormap that interpolates linearly between 4 colors;
    black, c0, c1 and c0+c1.

    Args:
        l: size of the colormap, defaults to 256
        c0: [r,g,b] array-like defining the color at the top left corner, defaults to [1,0.5,0] (orange)
        c1: [r,g,b] array-like defining the color at the bottom right corner, defaults to [0,0.5,1]] (light blue)
    returns:
        a (l,l,3) numpy array of rgb values
    '''
    rgb = np.zeros((l,l,3))
    rgb[:,:,:] = np.linspace(0,1,l)[:,np.newaxis,np.newaxis]*np.array(c0)[np.newaxis,np.newaxis,:]
    rgb[:,:,:] += np.linspace(0,1,l)[np.newaxis,:,np.newaxis]*np.array(c1)[np.newaxis,np.newaxis,:]
    rgb[rgb<0] = 0
    rgb[rgb>1] = 1
    return rgb
