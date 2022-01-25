# colorstamps

Python package for 2D colormaps. 
Most included colormaps are based on the 'CAM02-LCD' colorspace as defined in the package colorspacious (https://pypi.org/project/colorspacious/)

Documentation is hosted at https://colorstamps.readthedocs.io/en/latest/ but a summary follows below:


Installation:
pip install colorstamps

Example:
```
import matplotlib.pyplot as plt
import colorstamps

img = colorstamps.helpers.get_random_data() # numpy array of shape (100,200,2) with 2d data to plot    
rgb, stamp = colorstamps.apply_stamp(img[:,:,0], img[:,:,1], 'peak',
                                   vmin_0 = -1.2, vmax_0 = 1.2,
                                   vmin_1 = -1, vmax_1 = 1,
                                 )

fig, axes = plt.subplots(1,2,figsize=(10,3), dpi = 100)    
axes[0].imshow(rgb)

# show colormap as overlay
overlaid_ax = stamp.overlay_ax(axes[0], lower_left_corner = [0.7,0.85], width = 0.2)
overlaid_ax.set_ylabel(r'$\phi$')
overlaid_ax.set_xlabel(r'$\omega$')

# also show colormap as in separate ax to illustrate functionality
stamp.show_in_ax(axes[1])
axes[1].set_ylabel(r'$\phi$')
axes[1].set_xlabel(r'$\omega$')
```


![](docs/source/images/example0.png?raw=true)

Examples of included colormaps are shown below:

![](docs/source/images/colormaps.png?raw=true)
### Customization
All the colormaps are called by colorstamps.stamps.get_cmap() and the following keyword can be used with either colorstamps.apply_stamp() or colorstamps.stamps.get_cmap() to customize the colormaps:
```
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
Additionally, for radial colormaps the name may have a postfix separated by a space, i.e. 'cone tr'
  the postfix must be some combination of (t,b) and/or (l,r) which defines the quadrant/side of the colormap to include
  t-> top, b-> bottom, r-> right, l-> left, and 'tr'-> top right, etc.
```

### Evaluation of colormaps

The package also includes a method for evaluating colormaps, by evaluating the distribution of lighness, hue, saturation, 
and how it can be percieved by those with partial colorblindness
```
stamp = colorstamps.Stamp('hsv')
fig, ax = stamp.eval()
stamp = colorstamps.Stamp('peak')
fig, ax = stamp.eval()
```
![](docs/source/images/eval_hsv.png?raw=true)
![](docs/source/images/eval_peak.png?raw=true)

The package also supports different methods for clipping data to the colormap, using the 'clip' keyword in colorstamps.apply_stamp()

![](docs/source/images/point_outside_colormap.png?raw=true)

### Additional colormaps

Additional colormaps available by converting 1d colormaps in matplotlib to 2d colormaps by varying the lightness along the y-axis:

![](docs/source/images/mpl_colormaps.png?raw=true)

### User-provided colormaps
Custom colormaps may be integrated by providing a numpy array of shape (l,l,3) detailing a 2d colormap instead of a name when calling  colorstamps.apply_stamp()
```
my_cmap = np.zeros((256,256,3))
my_cmap[:,:,0] = np.linspace(0,1,256)[:,np.newaxis]
my_cmap[:,:,2] = np.linspace(0,1,256)[np.newaxis,:]
my_cmap[:,:,1] = 0.5*(my_cmap[:,:,0]+my_cmap[:,:,2])

img = colorstamps.helpers.get_random_data() # numpy array of shape (100,200,2) with 2d data to plot    
fig, ax = plt.subplots(1,1,figsize=(5,3), dpi = 100)    
rgb, stamp = colorstamps.apply_stamp(img[:,:,0], img[:,:,1], my_cmap)

ax.imshow(rgb)
overlaid_ax = stamp.overlay_ax(ax, lower_left_corner = [0.66,0.85], width = 0.2)
overlaid_ax.set_ylabel(r'$\phi$')
overlaid_ax.set_xlabel(r'$\omega$')
```
![](docs/source/images/custom_cmap.png?raw=true)

# support
Any contributions that would provide additional colormaps for are welcome, as are contributions to increase the functionality.
